from __future__ import division
import pyrticle.units as units
import pylinear.array as num
import pylinear.computation as comp




def uniform_on_unit_sphere(dim):
    from random import gauss

    # cf.
    # http://www-alg.ist.hokudai.ac.jp/~jan/randsphere.pdf
    # Algorith due to Knuth

    while True:
        pt = num.array([gauss(0,1) for i in range(dim)])
        n2 = comp.norm_2(pt)
        return pt/n2




def make_kv_distributed_particle(radii, emittances, vz, embed_dim=None):
    """Return (position, velocity) for a random particle
    according to a Kapchinskij-Vladimirskij distribution.
    """
    assert len(radii) == len(emittances)

    x = uniform_on_unit_sphere(len(radii) + len(emittances))
    pos = [xi*ri for xi, ri in zip(x[:len(radii)], radii)]
    momenta = [x_i/r_i*eps_i 
            for x_i, r_i, eps_i in 
            zip(x[len(radii):], radii, emittances)]

    one = sum(x_i**2/r_i**2 for x_i, r_i in zip(pos, radii)) + \
            sum(p_i**2*r_i**2/eps_i**2 
            for p_i, r_i, epsi in zip(momenta, radii, emittances))
    assert abs(one-1) < 1e-15

    if embed_dim is not None:
        assert embed_dim == int(embed_dim)

        while len(pos) < embed_dim:
            pos.append(0*units.M)
        while len(momenta) < embed_dim:
            momenta.append(0)

    z = num.array([0,0,1])

    return (num.array([x_i for x_i in pos]),
            z*vz + num.array([p_i*vz for p_i in momenta]))

    

def add_kv_xy_particles(nparticles, cloud, discr, 
        charge, mass, radii, emittances, beta, z_length, z_pos):
    from random import uniform
    from math import sqrt

    positions = []
    velocities = []

    bbox_min, bbox_max = discr.mesh.bounding_box
    center = (bbox_min+bbox_max)/2
    center[2] = 0
    size = bbox_max-bbox_min

    vz = beta*units.VACUUM_LIGHT_SPEED
    z = num.array([0,0,1])

    for i in range(nparticles):
        pos, v = make_kv_distributed_particle(
                radii, emittances, vz=beta*units.VACUUM_LIGHT_SPEED,
                embed_dim=cloud.mesh_info.dimensions)

        my_beta = comp.norm_2(v)/units.VACUUM_LIGHT_SPEED
        assert abs(beta - my_beta)/beta < 1e-4

        positions.append(center+pos+z*(z_pos+uniform(-z_length, z_length)/2))
        velocities.append(v)

    cloud.add_particles(positions, velocities, charge, mass)




class ODEDefinedFunction:
    def __init__(self, t0, y0, dt):
        self.t = [t0]
        self.y = [y0]
        self.dt = dt

        from hedge.timestep import RK4TimeStepper
        self.forward_stepper = RK4TimeStepper()
        self.backward_stepper = RK4TimeStepper()

    def __call__(self, t):
        def copy_if_necessary(x):
            try:
                return x[:]
            except TypeError:
                return x

        if t < self.t[0]:
            steps = int((self.t[0]-t)/self.dt)+2
            t_list = [self.t[0]]
            y_list = [self.y[0]]
            for n in range(steps):
                y_list.append(copy_if_necessary(self.backward_stepper(
                    y_list[-1], t_list[-1], -self.dt, self.rhs)))
                t_list.append(t_list[-1]-self.dt)

            self.t = t_list[:0:-1] + self.t
            self.y = y_list[:0:-1] + self.y
        elif t >= self.t[-1]:
            steps = int((t-self.t[-1])/self.dt)+2
            t_list = [self.t[-1]]
            y_list = [self.y[-1]]
            for n in range(steps):
                y_list.append(copy_if_necessary(self.forward_stepper(
                    y_list[-1], t_list[-1], self.dt, self.rhs)))
                t_list.append(t_list[-1]+self.dt)

            self.t = self.t + t_list[1:]
            self.y = self.y + y_list[1:]

        from bisect import bisect_right
        below_idx = bisect_right(self.t, t)-1
        assert below_idx >= 0
        above_idx = below_idx + 1

        assert above_idx < len(self.t)
        assert self.t[below_idx] <= t <= self.t[above_idx]

        # FIXME linear interpolation, bad
        slope = ((self.y[above_idx]-self.y[below_idx]) 
                /
                (self.t[above_idx]-self.t[below_idx]))

        return self.y[below_idx] + (t-self.t[below_idx]) * slope

    def rhs(self, t, y):
        raise NotImplementedError




class ChargelessKVRadiusPredictor:
    def __init__(self, a0, eps):
        self.a0 = a0
        self.eps = eps

    def __call__(self, s):
        from math import sqrt
        return sqrt(self.a0**2+(self.eps/self.a0)**2 * s**2)




class KVRadiusPredictor(ODEDefinedFunction):
    """Implement equation (1.74) in Alex Chao's book.

    See equation (1.65) for the definition of xi.
    """
    def __init__(self, a0, eps, eB_2E=0, xi=0):
        ODEDefinedFunction.__init__(self, 0, num.array([a0, 0]), 
                dt=1e-3*(a0**3/eps**2)**2)
        self.eps = eps
        self.xi = xi
        self.eB_2E = eB_2E

    def rhs(self, t, y):
        a = y[0]
        aprime = y[1]
        return num.array([
            aprime, 
            - self.eB_2E**2 * a
            + self.eps**2/a**3
            + self.xi/(2*a)
            ])

    def __call__(self, t):
        return ODEDefinedFunction.__call__(self, t)[0]




class BeamRadiusLoggerBase:
    def __init__(self, dimensions):
        self.dimensions = dimensions

        self.s_collector = []
        self.r_collector = []

        self.nparticles = 0

    def generate_plot(self, title, theory=None, outfile="beam-rad.eps"):
        from Gnuplot import Gnuplot, Data

        gp = Gnuplot()
        gp("set terminal postscript eps")
        gp("set output \"%s\"" % outfile)
        gp.title(title)
        gp.xlabel("s [m]")
        gp.ylabel("Beam Radius [m]")

        data = [Data(
            self.s_collector, 
            self.r_collector,
            title="simulation: %d particles" % self.nparticles,
            with_="lines")]

        if theory is not None:
            data.append(
                    Data(self.s_collector, 
                        [theory(s) for s in self.s_collector],
                        title="theoretical value",
                        with_="lines"))

        gp.plot(*data)

    def write_data(self, filename):
        outf = open(filename, "w")
        for s, r in zip(self.s_collector, self.r_collector):
            outf.write("%g\t%g\n" % (s, r))
        outf.close()

    def relative_error(self, theory):
        true_r = [theory(s) for s in self.s_collector]
        return max(abs(r-r0)/r0 for r, r0 in zip(self.r_collector, true_r))



class MaxBeamRadiusLogger(BeamRadiusLoggerBase):
    def update(self, t, positions, velocities):
        from math import sqrt,pi
        from pytools import argmax

        dim = self.dimensions
        nparticles = len(positions) // dim
        self.nparticles = max(nparticles, self.nparticles)
        pn = argmax(positions[i*dim+0]**2 +positions[i*dim+1]**2
                for i in xrange(nparticles))
        r = comp.norm_2(positions[pn*dim+0:pn*dim+2])
        vz = velocities[pn*dim+2]
        s = t*vz

        self.s_collector.append(s)
        self.r_collector.append(r)




class RMSBeamRadiusLogger(BeamRadiusLoggerBase):
    def update(self, t, positions, velocities):
        from math import sqrt,pi
        from pytools import average

        dim = self.dimensions
        nparticles = len(positions) // dim
        self.nparticles = max(nparticles, self.nparticles)
        r = sqrt(
                average(positions[i*dim+0]**2 +positions[i*dim+1]**2
                for i in xrange(nparticles))
                )
        vz = average(velocities[(dim-1)::dim])
        s = t*vz

        self.s_collector.append(s)
        self.r_collector.append(r)



if __name__ == "__main__":
    class Sin(ODEDefinedFunction):
        def __init__(self):
            ODEDefinedFunction.__init__(self, 0, 
                    num.array([0,1]), 1/7*1e-2)

        def rhs(self, t, y):
            return num.array([y[1], -y[0]])

        def __call__(self, t):
            return ODEDefinedFunction.__call__(self, t)[0]

    from math import pi
    s = Sin()
    assert abs(s(-pi)) < 2e-3
    assert abs(s(pi)) < 2e-3
    assert abs(s(-pi/2)+1) < 2e-3

    kv_env_exact = ChargelessKVRadiusPredictor(1, 1)
    kv_env_num = KVRadiusPredictor(1, 1)

    from hedge.tools import plot_1d
    steps = 20
    final_s = 4
    ds = final_s/steps
    for i in range(steps):
        s = ds*i
        
        a_exact = kv_env_exact(s)
        a_num = kv_env_num(s)
        assert abs(a_exact-a_num)/a_exact < 1e-3
