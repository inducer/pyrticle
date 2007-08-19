from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pyrticle.units as units
import cProfile as profile




def uniform_on_unit_sphere(dim):
    from random import gauss

    # cf.
    # http://www-alg.ist.hokudai.ac.jp/~jan/randsphere.pdf
    # Algorith due to Knuth

    while True:
        pt = num.array([gauss(0,1) for i in range(dim)])
        n2 = comp.norm_2(pt)
        return pt/n2




def make_kv_distributed_particle(radii, emittances, gamma, embed_dim=None):
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

    return (num.array([x_i for x_i in pos]),
            num.array([p_i*units.C0/gamma for p_i in momenta]))

    

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

    z = num.array([0,0,1])

    gamma = (1-beta**2)**-0.5

    for i in range(nparticles):
        pos, v = make_kv_distributed_particle(
                radii, emittances, gamma=gamma,
                embed_dim=cloud.mesh_info.dimensions)

        vz = sqrt(beta**2*units.C0**2 - comp.norm_2_squared(v))
        v += z*vz

        my_beta = comp.norm_2(v)/units.C0
        assert abs(beta - my_beta)/beta < 1e-15

        positions.append(center+pos+z*(z_pos+uniform(-z_length, z_length)/2))
        velocities.append(v)

    cloud.add_particles(positions, velocities, charge, mass)





class RTLogger:
    def __init__(self, dimensions, a0, eps):
        self.outf = open("particle-r-t.dat", "w")
        self.outf_theory = open("particle-r-t-theory.dat", "w")
        self.dimensions = dimensions

        self.a0 = a0
        self.eps = eps

    def __call__(self, t, positions, velocities):
        from math import sqrt,pi
        from pytools import argmax

        dim = self.dimensions
        nparticles = len(positions) // dim
        pn = argmax(positions[i*dim+0]**2 +positions[i*dim+1]**2
                for i in xrange(nparticles))
        r = comp.norm_2(positions[pn*dim+0:pn*dim+2])
        vz = velocities[pn*dim+2]
        s = t*vz
        self.outf.write("%g\t%g\n" % (s,r))
        self.outf.flush()

        a = sqrt(self.a0**2+(self.eps/self.a0)**2 * s**2)
        self.outf_theory.write("%g\t%g\n" % (s,a))
        self.outf_theory.flush()

            


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

    # discretization setup ----------------------------------------------------
    mesh = make_cylinder_mesh(radius=25*units.MM, height=100*units.MM, periodic=True)
    #mesh = make_box_mesh([1,1,2], max_volume=0.01)

    discr = Discretization(mesh, TetrahedralElement(3))
    vis = SiloVisualizer(discr)

    dt = discr.dt_factor(units.C0) / 2
    final_time = 0.2*units.M/units.C0
    nsteps = int(final_time/dt)+1
    dt = final_time/nsteps

    print "#elements=%d, dt=%s, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    mass = bind_mass_matrix(discr)

    def l2_norm(field):
        return sqrt(dot(field, mass*field))

    maxwell = MaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)

    # particles setup ---------------------------------------------------------
    nparticles = 1000

    cloud = ParticleCloud(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            verbose_vis=False)

    cloud_charge = 1e-9 * units.C
    particle_charge = cloud_charge/nparticles
    electrons_per_particle = cloud_charge/nparticles/units.EL_CHARGE
    print "e-/particle = ", electrons_per_particle 

    emittance = 5*pi * units.MM * units.MRAD
    initial_radius = 2.5*units.MM

    #el_energy = 1e6 * units.EV
    #el_energy = units.EL_REST_ENERGY*1.00001
    #el_lorentz_gamma = el_energy/units.EL_REST_ENERGY
    el_lorentz_gamma = 100000
    beta = (1-1/el_lorentz_gamma**2)**0.5
    print "v = %g%% c" % (beta*100)

    add_kv_xy_particles(nparticles, cloud, discr, 
            charge=0, 
            mass=electrons_per_particle*units.EL_MASS,
            radii=[2.5*units.MM, 2.5*units.MM],
            z_length=5*units.MM,
            z_pos=10*units.MM,
            emittances=[emittance, emittance], 
            beta=beta)

    full_charge_at_time = final_time*0.2

    # timestepping ------------------------------------------------------------
    fields = concatenate_fields(
            [discr.volume_zeros() for i in range(2*discr.dimensions)],
            [cloud])

    def rhs(t, y):
        e = y[:3]
        h = y[3:6]

        if False:
            maxwell_rhs = maxwell.rhs(t, y[0:6])
            rho, j = cloud.reconstruct_densities()
            cloud_rhs = [cloud.rhs(t, e, h)]
        else:
            maxwell_rhs = ArithmeticList(6*[discr.volume_zeros()])
            rho = discr.volume_zeros()
            j = ArithmeticList(3*[discr.volume_zeros()])
            cloud_rhs = [ArithmeticList([
                cloud.velocities, 
                0*cloud.velocities, 
                ])]

        rhs_e = maxwell_rhs[:3]
        rhs_h = maxwell_rhs[3:6]
        return concatenate_fields(
                rhs_e + 1/units.EPSILON0*j,
                rhs_h,
                cloud_rhs,
                )

    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    rt_logger = RTLogger(cloud.mesh_info.dimensions,
            initial_radius, emittance)

    for step in xrange(nsteps):
        if False:
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

        rt_logger(t, cloud.positions, cloud.velocities)

        if False:
            myfields = [fields]
            fields = profile.runctx("myfields[0] = stepper(fields, t, dt, rhs)", 
                    globals(), locals(), "pic-%04d.prof" % step)
            fields = myfields[0]
        else:
            fields = stepper(fields, t, dt, rhs)

        cloud.upkeep()

        print "timestep %d, t=%g l2[e]=%g l2[h]=%g secs=%f particles=%d" % (
                step, t, l2_norm(fields[0:3]), l2_norm(fields[3:6]),
                time()-last_tstep, len(cloud))
        if False:
            print "searches: same=%d, normal=%d, vertex=%d, global=%d, periodic=%d" % (
                    cloud.same_searches.pop(),
                    cloud.normal_searches.pop(),
                    cloud.vertex_searches.pop(),
                    cloud.global_searches.pop(),
                    cloud.periodic_hits.pop(),
                    )
            print "shape-adds: neighbor=%d vertex=%d" % (
                    cloud.neighbor_shape_adds.pop(),
                    cloud.vertex_shape_adds.pop(),
                    )

        last_tstep = time()

        t += dt

        if t < full_charge_at_time:
            charge_now = t/full_charge_at_time*particle_charge
            cloud.charges = charge_now * \
                    num.ones((len(cloud.containing_elements),))
             




if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
