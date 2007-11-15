from __future__ import division
import pylinear.array as num
import pylinear.computation as comp




def uniform_on_unit_sphere(dim):
    from random import gauss

    # cf.
    # http://www-alg.ist.hokudai.ac.jp/~jan/randsphere.pdf
    # Algorith due to Knuth

    pt = num.array([gauss(0,1) for i in range(dim)])
    n2 = comp.norm_2(pt)
    return pt/n2





    

class KVZIntervalBeam:
    def __init__(self, units, nparticles, p_charge, p_mass,
            radii, emittances, beta, z_length, z_pos):
        """Construct a beam that is KV-distributed in (x,y)
        and uniform over an interval along z.

        @par radii: total (100%) radii
        @par emittances: total (100%) emittances
        """

        assert len(radii) == len(emittances)

        self.units = units

        self.nparticles = nparticles

        self.p_charge = p_charge
        self.p_mass = p_mass

        self.radii = radii
        self.emittances = emittances

        self.beta = beta
        self.gamma = (1-beta)**(-0.5)

        self.z_length = z_length
        self.z_pos = z_pos

    @property
    def rms_radii(self):
        return [r/2 for r in self.radii]

    @property
    def rms_emittances(self):
        return [eps/4 for eps in self.emittances]

    def make_particle(self, vz, dim_x, dim_p):
        """Return (position, velocity) for a random particle
        according to a Kapchinskij-Vladimirskij distribution.
        """

        s = uniform_on_unit_sphere(len(self.radii) + len(self.emittances))
        x = [x_i*r_i for x_i, r_i in zip(s[:len(self.radii)], self.radii)]
        # xp like xprime
        xp = [s_i/r_i*eps_i 
                for s_i, r_i, eps_i in 
                zip(s[len(self.radii):], self.radii, self.emittances)]

        one = sum(x_i**2/r_i**2 for x_i, r_i in zip(x, self.radii)) + \
                sum(xp_i**2*r_i**2/eps_i**2 
                for xp_i, r_i, epsi in zip(xp, self.radii, self.emittances))
        assert abs(one-1) < 1e-15

        while len(x) < dim_x:
            x.append(0)
        while len(xp) < dim_p:
            xp.append(0)

        z = num.array([0,0,1])

        return (num.array(x),
                z*vz + num.array([xp_i*vz for xp_i in xp]))

    def add_to(self, cloud, discr):
        from random import uniform
        from math import sqrt

        positions = []
        velocities = []

        bbox_min, bbox_max = discr.mesh.bounding_box
        center = (bbox_min+bbox_max)/2
        center[2] = 0
        size = bbox_max-bbox_min

        vz = self.beta*self.units.VACUUM_LIGHT_SPEED
        z = num.array([0,0,1])

        for i in range(self.nparticles):
            pos, v = self.make_particle(
                    vz=self.beta*self.units.VACUUM_LIGHT_SPEED,
                    dim_x=cloud.dimensions_pos, dim_p=cloud.dimensions_velocity)

            my_beta = comp.norm_2(v)/self.units.VACUUM_LIGHT_SPEED
            assert abs(self.beta - my_beta)/self.beta < 1e-4

            positions.append(center
                    +pos
                    +z*(self.z_pos+uniform(-self.z_length, self.z_length)/2))
            velocities.append(v)

        cloud.add_particles(positions, velocities, 
                self.p_charge, self.p_mass)

    def get_space_charge_parameter(self):
        from math import pi

        Q = self.p_charge / self.units.EL_CHARGE

        # see (1.3) in Alex Wu Chao and 
        # http://en.wikipedia.org/wiki/Classical_electron_radius
        r0 = 1/(4*pi*self.units.EPSILON0)*( 
                (self.units.EL_CHARGE**2)
                /
                (self.units.EL_MASS*self.units.VACUUM_LIGHT_SPEED**2))

        lambda_ = self.nparticles/self.z_length
        A = self.p_mass / self.units.EL_MASS

        xi = ((4 * Q**2 * r0 * lambda_)
                /
                (A * self.beta**2 * self.gamma**2))
        return xi




def calculate_rms_beam_size(cloud, axis):
    from pytools import average
    from math import sqrt

    xdim = cloud.dimensions_pos

    return sqrt(average(x**2 for x in cloud.positions[axis::xdim]))




def calculate_rms_emittance(cloud, axis):
    from pytools import average
    from math import sqrt

    xdim = cloud.dimensions_pos
    vdim = cloud.dimensions_velocity

    xprime = [px/pz for px, pz in zip(cloud.momenta[axis::vdim], cloud.momenta[2::vdim])]
    mean_x_squared = average(x**2 for x in cloud.positions[axis::xdim])
    mean_p_squared = average(xp**2 for xp in xprime)
    squared_mean_xp = average(
            x*xp 
            for x, xp in zip(cloud.positions[axis::xdim], xprime)
            )**2

    return sqrt(mean_x_squared*mean_p_squared - squared_mean_xp)




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

    See equation (1.65) for the definition of M{xi}.
    M{Q} is the number of electrons in the beam
    """
    def __init__(self, a0, eps, eB_2E=0, xi=0, dt=1e-4):
        ODEDefinedFunction.__init__(self, 0, num.array([a0, 0]), 
                dt=dt*(a0**4/eps**2)**2)
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

    def generate_plot(self, title, sim_label, theories, outfile="beam-rad"):
        theory_data = [
            (name, [theory(s) for s in self.s_collector])
            for name, theory in theories]
        
        try:
            from Gnuplot import Gnuplot, Data

            gp = Gnuplot()
            gp("set terminal postscript eps")
            gp("set output \"%s.eps\"" % outfile)
            gp.title(title)
            gp.xlabel("s [m]")
            gp.ylabel("Beam Radius [m]")

            data = [Data(
                self.s_collector, 
                self.r_collector,
                title=sim_label+" [%d particles]" % self.nparticles,
                with_="lines")]

            for name, data in theory_data:
                data.append(
                        Data(self.s_collector, data, title=name,
                            with_="lines"))

            gp.plot(*data)
        except ImportError:
            pass 

        write_data_file(
            "%s-sim.dat" % outfile,
            self.s_collector,
            self.r_collector,
            sim_label)

        for i, (name, data) in enumerate(theory_data):
            write_data_file(
                "%s-theory-%d.dat" % (outfile, i),
                self.s_collector,
                data,
                name)

    def relative_error(self, theory):
        true_r = [theory(s) for s in self.s_collector]
        return max(abs(r-r0)/r0 for r, r0 in zip(self.r_collector, true_r))




def write_data_file(filename, x, y, comment=None):
    outf = open(filename, "w")

    if comment is not None:
        outf.write("# %s\n" % comment)

    for x_i, y_i in zip(x, y):
        outf.write("%g\t%g\n" % (x_i, y_i))
    outf.close()




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
    def __init__(self, dimensions, axis):
        BeamRadiusLoggerBase.__init__(self, dimensions)
        self.axis = axis

    def update(self, t, positions, velocities):
        from math import sqrt,pi
        from pytools import average

        dim = self.dimensions
        nparticles = len(positions) // dim
        self.nparticles = max(nparticles, self.nparticles)
        
        r = sqrt(average(
            positions[i*dim+self.axis]**2
            for i in xrange(nparticles)
            ))

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

    kv_env_exact = ChargelessKVRadiusPredictor(2.5e-3, 5e-6)
    kv_env_num = KVRadiusPredictor(2.5e-3, 5e-6)

    from hedge.tools import plot_1d
    steps = 50
    for i in range(steps):
        s = kv_env_num.dt/7*i
        
        a_exact = kv_env_exact(s)
        a_num = kv_env_num(s)
        assert abs(a_exact-a_num)/a_exact < 1e-3
