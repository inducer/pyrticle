from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pylinear.operator as op
import cProfile as profile




def add_gauss_particles(nparticles, cloud, discr, charge, mass, 
        mean_x, mean_p, sigma_x, sigma_p):
    from random import gauss
    from pyrticle.cloud import v_from_p

    cloud.add_particles(
            positions=[
                num.array([gauss(m, s) for m, s in zip(mean_x, sigma_x)]) 
                for i in range(nparticles)
                ],
            velocities=[v_from_p(
                num.array([gauss(m, s) for m, s in zip(mean_p, sigma_p)]),
                mass, cloud.units.VACUUM_LIGHT_SPEED) 
                for i in range(nparticles)
                ],
            charges=charge, masses=mass)



def main():
    from hedge.element import TriangularElement
    from hedge.timestep import RK4TimeStepper
    from hedge.mesh import \
            make_square_mesh, \
            make_regular_square_mesh, \
            make_regular_rect_mesh
    from hedge.discretization import \
            Discretization, \
            pair_with_boundary
    from hedge.visualization import VtkVisualizer, SiloVisualizer
    from hedge.tools import dot
    from math import sqrt, pi
    from pytools.arithmetic_container import \
            ArithmeticList, join_fields
    from hedge.operators import TEMaxwellOperator, DivergenceOperator
    from pyrticle.cloud import ParticleCloud
    from random import seed
    #seed(0)

    from pyrticle.units import SI
    units = SI()

    # discretization setup ----------------------------------------------------
    full_mesh = make_regular_rect_mesh(
            a=(-0.5, -0.5),
            b=(1.5, 0.5),
            #n=(30, 10), 
            n=(10, 5), 
            periodicity=(True, False))

    #full_mesh = make_regular_square_mesh(n=2, periodicity=(True, False))

    from hedge.parallel import guess_parallelization_context

    pcon = guess_parallelization_context()

    if pcon.is_head_rank:
        mesh = pcon.distribute_mesh(full_mesh)
    else:
        mesh = pcon.receive_mesh()

    discr = pcon.make_discretization(mesh, TriangularElement(5))
    vis = SiloVisualizer(discr)
    #vis = VtkVisualizer(discr, "pic")

    max_op = TEMaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)
    div_op = DivergenceOperator(discr)

    dt = discr.dt_factor(max_op.c) / 2
    final_time = 2*units.M/max_op.c
    nsteps = int(final_time/dt)+1
    dt = final_time/nsteps

    print "#elements=%d, dt=%s, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    def l2_norm(field):
        return sqrt(dot(field, discr.mass_operator*field))
    def l2_error(field, true):
        return l2_norm(field-true)/l2_norm(true)

    # particles setup ---------------------------------------------------------
    nparticles = 1000

    cloud = ParticleCloud(discr, units, dimensions_pos=2, dimensions_velocity=2,
            verbose_vis=True)

    cloud_charge = -1e-9 * units.C
    electrons_per_particle = abs(cloud_charge/nparticles/units.EL_CHARGE)
    print "e-/particle = ", electrons_per_particle 

    avg_x_vel = 0.99*units.VACUUM_LIGHT_SPEED
    mean_v = num.array([avg_x_vel, 0])
    mean_beta = mean_v/units.VACUUM_LIGHT_SPEED
    gamma = units.gamma(mean_v)
    pmass = electrons_per_particle*units.EL_MASS
    mean_p = gamma*pmass*mean_v

    print "beta=%g, gamma=%g" % (comp.norm_2(mean_beta), gamma)

    add_gauss_particles(nparticles, cloud, discr, 
            charge=cloud_charge/nparticles, 
            mass=pmass,
            mean_x=num.zeros((2,)),
            mean_p=mean_p,
            sigma_x=0.1*num.ones((2,)),
            sigma_p=units.gamma(mean_v)*pmass*num.ones((2,))*avg_x_vel*0.1,
            )

    # intial condition --------------------------------------------------------
    from pyrticle.cloud import compute_initial_condition
    fields = compute_initial_condition(pcon, discr, cloud, mean_beta, max_op,
            debug=True)

    # timestepping ------------------------------------------------------------
    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    for step in xrange(nsteps):
        cloud.upkeep()
        fields = stepper(fields, t, dt, fields.rhs)

        if True:
            visf = vis.make_file("pic-%04d" % step)

            mesh_scalars, mesh_vectors = \
                    cloud.add_to_vis(vis, visf, time=t, step=step)
            vis.add_data(visf, [
                        ("divD", max_op.epsilon*div_op(fields.e)),
                        ("e", fields.e), 
                        ("h", fields.h), 
                        ] + mesh_scalars + mesh_vectors,
                    time=t, step=step)
            visf.close()

        print "timestep %d, t=%g l2[e]=%g l2[h]=%g secs=%f particles=%d" % (
                step, t, l2_norm(fields.e), l2_norm(fields.h),
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

    vis.close()




if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
