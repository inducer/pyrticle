from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pylinear.operator as op
import cProfile as profile




def add_gauss_particles(nparticles, cloud, discr, charge, mass, 
        mean_x, mean_p, sigma_x, sigma_p):
    from random import gauss

    cloud.add_particles(
            positions=[
                num.array([gauss(m, s) for m, s in zip(mean_x, sigma_x)]) 
                for i in range(nparticles)
                ],
            velocities=[cloud.units.v_from_p(mass, 
                num.array([gauss(m, s) for m, s in zip(mean_p, sigma_p)])) 
                for i in range(nparticles)
                ],
            charges=charge, masses=mass)



def main():
    from hedge.element import TriangularElement
    from hedge.timestep import RK4TimeStepper
    from hedge.mesh import \
            make_square_mesh, \
            make_regular_square_mesh, \
            make_regular_rect_mesh, \
            make_rect_mesh
    from hedge.discretization import \
            Discretization, \
            pair_with_boundary
    from hedge.visualization import VtkVisualizer, SiloVisualizer
    from hedge.tools import dot
    from math import sqrt, pi
    from pytools.arithmetic_container import \
            ArithmeticList, join_fields
    from hedge.operators import TEMaxwellOperator, DivergenceOperator

    from random import seed
    seed(0)

    from pyrticle.units import SI
    units = SI()

    # discretization setup ----------------------------------------------------
    tube_length = 2
    full_mesh = make_rect_mesh(
            a=(-0.5, -0.5),
            b=(-0.5+tube_length, 0.5),
            periodicity=(True, False),
            subdivisions=(10,5),
            max_area=0.02)

    #full_mesh = make_regular_square_mesh(n=2, periodicity=(True, False))

    from hedge.parallel import guess_parallelization_context

    pcon = guess_parallelization_context()

    if pcon.is_head_rank:
        mesh = pcon.distribute_mesh(full_mesh)
    else:
        mesh = pcon.receive_mesh()

    discr = pcon.make_discretization(mesh, TriangularElement(7))
    vis = SiloVisualizer(discr)
    #vis = VtkVisualizer(discr, "pic")

    max_op = TEMaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)
    div_op = DivergenceOperator(discr)

    dt = discr.dt_factor(max_op.max_eigenvalue())
    #final_time = 15*units.M/max_op.c
    final_time = 50*units.M/max_op.c
    nsteps = int(final_time/dt)+1
    dt = final_time/nsteps

    print "#elements=%d, dt=%s, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    def l2_norm(field):
        return sqrt(dot(field, discr.mass_operator*field))
    def l2_error(field, true):
        return l2_norm(field-true)/l2_norm(true)

    # particles setup ---------------------------------------------------------
    nparticles = 100

    from pyrticle.cloud import ParticleCloud
    from pyrticle.reconstruction import \
            ShapeFunctionReconstructor, \
            NormalizedShapeFunctionReconstructor, \
            AdvectiveReconstructor
    from pyrticle.pusher import MonomialParticlePusher
    cloud = ParticleCloud(discr, units, 
            #AdvectiveReconstructor(
                #activation_threshold=1e-5,
                #kill_threshold=1e-3,
                #upwind_alpha=1),
            ShapeFunctionReconstructor(),
            #NormalizedShapeFunctionReconstructor(),
            MonomialParticlePusher(),
            dimensions_pos=2, dimensions_velocity=2,
            verbose_vis=True)

    cloud_charge = -1e-9 * units.C
    electrons_per_particle = abs(cloud_charge/nparticles/units.EL_CHARGE)
    print "e-/particle = ", electrons_per_particle 

    avg_x_vel = 0.90*units.VACUUM_LIGHT_SPEED
    mean_v = num.array([avg_x_vel, 0])
    #mean_v = num.array([avg_x_vel*0.5, avg_x_vel*0.8])
    mean_beta = mean_v/units.VACUUM_LIGHT_SPEED
    gamma = units.gamma(mean_v)
    pmass = electrons_per_particle*units.EL_MASS
    mean_p = gamma*pmass*mean_v

    sigma_v = num.array([avg_x_vel*1e-3, avg_x_vel*1e-6])
    print "beta=%g, gamma=%g" % (comp.norm_2(mean_beta), gamma)

    add_gauss_particles(nparticles, cloud, discr, 
            charge=cloud_charge/nparticles, 
            mass=pmass,
            mean_x=num.zeros((2,)),
            mean_p=mean_p,
            sigma_x=0.1*num.ones((2,)),
            sigma_p=units.gamma(mean_v)*pmass*sigma_v,
            )

    # intial condition --------------------------------------------------------
    from pyrticle.cloud import compute_initial_condition
    fields = compute_initial_condition(pcon, discr, cloud, mean_beta, max_op,
            debug=True)

    stepper = RK4TimeStepper()

    # diagnostics setup -------------------------------------------------------
    from pytools.log import LogManager, \
            add_simulation_quantities, \
            add_general_quantities, \
            add_run_info, ETA
    from pyrticle.log import add_particle_quantities, add_field_quantities, \
            add_beam_quantities, add_currents
    logmgr = LogManager("2d.dat")
    add_run_info(logmgr)
    add_general_quantities(logmgr)
    add_simulation_quantities(logmgr, dt)
    add_particle_quantities(logmgr, cloud)
    add_field_quantities(logmgr, fields, reconstruct_interval=1)
    add_beam_quantities(logmgr, cloud, axis=1, beam_axis=0)
    add_currents(logmgr, fields, (1,0), tube_length)

    stepper.add_instrumentation(logmgr)
    fields.add_instrumentation(logmgr)
    logmgr.set_constant("beta", comp.norm_2(mean_beta))
    logmgr.set_constant("gamma", gamma)
    logmgr.set_constant("vx", avg_x_vel)
    logmgr.set_constant("Q0", cloud_charge)
    logmgr.set_constant("n_part_0", nparticles)
    logmgr.set_constant("pmass", electrons_per_particle*units.EL_MASS)

    from pytools.log import IntervalTimer
    vis_timer = IntervalTimer("t_vis", "Time spent visualizing")
    logmgr.add_quantity(vis_timer)

    logmgr.add_quantity(ETA(nsteps))

    logmgr.add_watches(["step", "t_sim", "W_field", "t_step", "t_eta", "n_part"])

    # timestepping ------------------------------------------------------------
    t = 0


    for step in xrange(nsteps):
        logmgr.tick()

        cloud.upkeep()
        fields = stepper(fields, t, dt, fields.rhs)

        if step % 20 == 0:
            vis_timer.start()
            visf = vis.make_file("pic-%04d" % step)

            cloud.add_to_vis(vis, visf, time=t, step=step)
            vis.add_data(visf, [
                        ("divD", max_op.epsilon*div_op(fields.e)),
                        ("e", fields.e), 
                        ("h", fields.h), 

                        #("active_elements", 
                            #cloud.pic_algorithm.get_debug_quantity_on_mesh(
                                #"active_elements", cloud.raw_velocities())),
                        ("rho", cloud.reconstruct_rho()),
                        ("j", cloud.reconstruct_j()), 
                        ],
                        time=t, step=step,
                        expressions=[
                            ("rpa", "6e-9*active_elements+rho")
                            ])
            visf.close()
            vis_timer.stop()


        t += dt

    vis.close()

    logmgr.tick()
    logmgr.save()




if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
