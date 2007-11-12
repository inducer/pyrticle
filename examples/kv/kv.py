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




class KVRadiusPredictor:
    def __init__(self, a0, eps):
        self.a0 = a0
        self.eps = eps

    def __call__(self, s):
        from math import sqrt
        return sqrt(self.a0**2+(self.eps/self.a0)**2 * s**2)




class BeamRadiusLoggerBase:
    def __init__(self, dimensions, a0, eps):
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
