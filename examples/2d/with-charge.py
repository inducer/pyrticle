from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pylinear.operator as op
import cProfile as profile
import pytools




class GaussParticleDistribution(pytools.Record):
    def __init__(self, total_charge, total_mass, mean_x, mean_p, sigma_x, sigma_p):
        pytools.Record.__init__(self, locals())

    def add_to(self, cloud, nparticles):
        from random import gauss

        pmass = self.total_mass/nparticles
        cloud.add_particles(
                positions=[
                    num.array([gauss(m, s) for m, s in zip(self.mean_x, self.sigma_x)]) 
                    for i in range(nparticles)
                    ],
                velocities=[cloud.units.v_from_p(pmass, 
                    num.array([gauss(m, s) for m, s in zip(self.mean_p, self.sigma_p)])) 
                    for i in range(nparticles)
                    ],
                charges=self.total_charge/nparticles, 
                masses=pmass)

    def analytic_rho(self, discr):
        from math import exp, pi

        sigma_mat = num.diagonal_matrix(num.power(self.sigma_x, 2))
        inv_sigma_mat = num.diagonal_matrix(num.power(self.sigma_x, -2))

        def distrib(x):
            return 1/((2*pi)**(len(x)/2) * comp.determinant(sigma_mat)**0.5) \
                    * exp(-0.5*(x-self.mean_x)*inv_sigma_mat*(x-self.mean_x))

        rho = self.total_charge * discr.interpolate_volume_function(distrib)

        # check for correctness
        from hedge.discretization import integral
        int_rho = integral(discr, rho)
        rel_err = (int_rho-self.total_charge)/self.total_charge
        assert rel_err < 1e-6

        return rho




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
    from pytools.arithmetic_container import join_fields

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

    from hedge.operators import \
            TEMaxwellOperator, \
            DivergenceOperator, \
            StrongWaveOperator
    from pyrticle.hyperbolic import \
            ECleaningMaxwellOperator, \
            BoneHeadedCleaningMaxwellOperator
    from hedge.mesh import TAG_ALL, TAG_NONE

    max_op = TEMaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)
    max_op = ECleaningMaxwellOperator(max_op, chi=2)
    #wave_op = StrongWaveOperator(c=max_op.c*2, discr=discr,
            #dirichlet_tag=TAG_NONE, radiation_tag=TAG_ALL)

    #max_op = BoneHeadedCleaningMaxwellOperator(max_op, wave_op)
    div_op = DivergenceOperator(discr)

    dt = discr.dt_factor(max_op.max_eigenvalue())/ 10
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
    def make_cloud():
        from pyrticle.cloud import ParticleCloud
        from pyrticle.reconstruction import \
                ShapeFunctionReconstructor, \
                NormalizedShapeFunctionReconstructor, \
                AdvectiveReconstructor
        from pyrticle.pusher import \
                MonomialParticlePusher, \
                AverageParticlePusher

        import sys
        reconstructor_str = sys.argv[1]
        pusher_str = sys.argv[2]

        if reconstructor_str == "advective":
            reconstructor = AdvectiveReconstructor(
                    activation_threshold=1e-5,
                    kill_threshold=1e-3,
                    upwind_alpha=1)
        elif reconstructor_str == "shape":
            reconstructor = ShapeFunctionReconstructor()
        elif reconstructor_str == "normshape":
            reconstructor = NormalizedShapeFunctionReconstructor()
        else:
            raise ValueError, "invalid reconstructor"

        if pusher_str == "monomial":
            pusher = MonomialParticlePusher()
        elif pusher_str == "average":
            pusher = AverageParticlePusher()
        else:
            raise ValueError, "invalid pusher"

        return ParticleCloud(discr, units, reconstructor, pusher,
                dimensions_pos=2, dimensions_velocity=2,
                verbose_vis=True)

    cloud = make_cloud()

    nparticles = 1000

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

    gauss_p = GaussParticleDistribution(
            total_charge=cloud_charge, 
            total_mass=pmass*nparticles,
            mean_x=num.zeros((2,)),
            mean_p=mean_p,
            sigma_x=0.1*num.ones((2,)),
            sigma_p=units.gamma(mean_v)*pmass*sigma_v)
    gauss_p.add_to(cloud, nparticles)
    from pyrticle.cloud import optimize_shape_bandwidth
    optimize_shape_bandwidth(cloud, discr, gauss_p.analytic_rho(discr))

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

    substep = [0]

    for step in xrange(nsteps):
        logmgr.tick()

        #if step % 1 == 0:

        def rhs(t, y):
            fields = y
            vis_timer.start()
            visf = vis.make_file("pic-%04d" % substep[0])

            cloud.add_to_vis(vis, visf, time=t, step=substep[0])
            vis.add_data(visf, [
                        ("divD", max_op.epsilon*div_op(fields.e)),
                        ("e", fields.e), 
                        ("h", fields.h), 
                        ("phi", fields.phi), 

                        #("active_elements", 
                            #cloud.pic_algorithm.get_debug_quantity_on_mesh(
                                #"active_elements", cloud.raw_velocities())),
                        ("rho", cloud.reconstruct_rho()),
                        ("j", cloud.reconstruct_j()), 
                        ],
                        time=t, step=step,
                        expressions=[
                            ])
            visf.close()
            vis_timer.stop()

            substep[0] += 1



            return fields.rhs(t, y)

        cloud.upkeep()
        fields = stepper(fields, t, dt, rhs)


        t += dt

    vis.close()

    logmgr.tick()
    logmgr.save()




if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
