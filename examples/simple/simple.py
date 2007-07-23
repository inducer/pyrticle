from __future__ import division
import pylinear.array as num
import pylinear.computation as comp




def uniform_on_unit_sphere(dim):
    from random import uniform

    while True:
        pt = num.array([uniform(-1,1) for i in range(dim)])
        n2 = comp.norm_2(pt)

        if n2 <= 1:
            break

    return pt/n2




def make_kv_distributed_particle(radii, emittances, vz, c, embed_dim=None):
    """Return (position, velocity) for a random particle
    according to a Kapchinskij-Vladimirskij distribution.
    """
    from math import sqrt

    assert len(radii) == len(emittances)

    x = uniform_on_unit_sphere(len(radii) + len(emittances))
    pos = [xi*ri for xi, ri in zip(x[:len(radii)], radii)]
    momenta = [xi/ri*epsi 
            for xi, ri, epsi in 
            zip(x[len(radii):], radii, emittances)]

    one = sum(xi**2/ri**2 for xi, ri in zip(pos, radii)) + \
            sum(pi**2*ri**2/epsi**2 
            for pi, ri, epsi in zip(momenta, radii, emittances))
    assert abs(one-1) < 1e-15

    if embed_dim is not None:
        assert embed_dim == int(embed_dim)

        while len(pos) < embed_dim:
            pos.append(0)
        while len(momenta) < embed_dim:
            momenta.append(0)

    beta = vz/c
    gamma = 1/sqrt(1-beta**2)

    return num.array(pos), num.array(momenta)*c/gamma

    

def add_kv_xy_particles(nparticles, cloud, discr, 
        charge, mass, emittances, vz, c):
    from random import uniform

    positions = []
    velocities = []

    bbox_min, bbox_max = discr.mesh.bounding_box
    center = (bbox_min+bbox_max)/2
    size = bbox_max-bbox_min

    z = num.array([0,0,1])

    for i in range(nparticles):
        pos, v = make_kv_distributed_particle(
                size[:2]/3, emittances, vz, c,
                embed_dim=cloud.mesh_info.dimensions)

        positions.append(center+pos+z*uniform(-size[2]/2, size[2]/2))
        velocities.append(v+z*vz)

    cloud.add_particles(positions, velocities, charge, mass)





class RTLogger:
    def __init__(self, dimensions):
        self.outf = open("particle-r-t.dat", "w")
        self.dimensions = dimensions

    def __call__(self, t, positions):
        dim = self.dimensions
        nparticles = len(positions) // dim
        for i in xrange(nparticles):
            pstart = i*dim
            pend = (i+1)*dim
            r = comp.norm_2(positions[pstart:pend])
            self.outf.write("%g\t%g\n" % (t,r))
        self.outf.flush()

            


def main():
    from hedge.element import TetrahedralElement
    from hedge.timestep import RK4TimeStepper
    from hedge.mesh import \
            make_box_mesh, \
            make_cylinder_mesh
    from hedge.discretization import \
            Discretization, \
            bind_flux, \
            bind_nabla, \
            bind_mass_matrix, \
            bind_inverse_mass_matrix, \
            pair_with_boundary
    from hedge.visualization import SiloVisualizer
    from hedge.silo import SiloFile
    from hedge.silo import DB_VARTYPE_VECTOR
    from hedge.tools import dot, cross
    from math import sqrt, pi
    from pytools.arithmetic_container import \
            ArithmeticList, concatenate_fields
    from hedge.operators import MaxwellOperator
    from pyrticle.cloud import ParticleCloud
    from random import seed
    seed(0)

    epsilon0 = 8.8541878176e-12 # C**2 / (N m**2)
    mu0 = 4*pi*1e-7 # N/A**2.
    epsilon = 1*epsilon0
    mu = 1*mu0
    c = 1/sqrt(epsilon*mu)

    el_mass = 9.10938215e-31 # kg
    el_charge = 1.602176487e-19 # C

    mesh = make_cylinder_mesh(radius=1, height=2)
    #mesh = make_box_mesh([1,1,2], max_volume=0.01)

    discr = Discretization(mesh, TetrahedralElement(7))
    vis = SiloVisualizer(discr)

    dt = discr.dt_factor(1/sqrt(mu*epsilon))
    final_time = 2/c 
    nsteps = int(final_time/dt)+1
    dt = final_time/nsteps

    print "#elements=%d, dt=%g, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    mass = bind_mass_matrix(discr)

    def l2_norm(field):
        return sqrt(dot(field, mass*field))

    maxwell = MaxwellOperator(discr, epsilon, mu, upwind_alpha=0)

    nparticles = 2000
    cloud = ParticleCloud(discr, epsilon, mu, verbose_vis=False)

    add_kv_xy_particles(nparticles, cloud, discr, 
            charge=el_charge, mass=el_mass,
            emittances=[0.5,0.5], vz=0.5*c, c=c)

    fields = concatenate_fields(
            [discr.volume_zeros() for i in range(2*discr.dimensions)],
            [cloud])

    def rhs(t, y):
        e = y[:3]
        h = y[3:6]

        if False:
            maxwell_rhs = maxwell.rhs(t, y[0:6])
            rho, j = cloud.reconstruct_densities()
        else:
            maxwell_rhs = ArithmeticList(6*[discr.volume_zeros()])
            rho = discr.volume_zeros()
            j = ArithmeticList(3*[discr.volume_zeros()])

        rhs_e = maxwell_rhs[:3]
        rhs_h = maxwell_rhs[3:6]
        return concatenate_fields(
                rhs_e + 1/epsilon*j,
                rhs_h,
                [cloud.rhs(t, e, h)]
                )

    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    rt_logger = RTLogger(cloud.mesh_info.dimensions)

    for step in xrange(nsteps):
        print "timestep %d, t=%g l2[e]=%g l2[h]=%g secs=%f particles=%d" % (
                step, t, l2_norm(fields[0:3]), l2_norm(fields[3:6]),
                time()-last_tstep, len(cloud))
        last_tstep = time()

        if True:
            silo = SiloFile("pic-%04d.silo" % step)

            mesh_scalars, mesh_vectors = \
                    cloud.add_to_silo(silo, time=t, step=step)
            vis.add_to_silo(silo,
                    scalars=mesh_scalars,
                    vectors=[("e", fields[0:3]), 
                        ("h", fields[3:6]), ]
                    + mesh_vectors
                    ,
                    expressions=[
                        ],
                    write_coarse_mesh=True,
                    time=t, step=step)
            silo.close()

            rt_logger(t, cloud.positions)

        fields = stepper(fields, t, dt, rhs)
        cloud.upkeep()


        t += dt




if __name__ == "__main__":
    #import cProfile as profile
    #profile.run("main()", "pic.prof")
    main()
