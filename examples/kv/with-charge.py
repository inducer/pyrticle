from __future__ import division
import pylinear.array as num
import pylinear.computation as comp




def main():
    from hedge.element import TetrahedralElement
    from hedge.timestep import RK4TimeStepper
    from hedge.visualization import VtkVisualizer, SiloVisualizer
    from hedge.tools import dot
    from math import sqrt
    from hedge.operators import MaxwellOperator, DivergenceOperator
    from kv import KVZIntervalBeam
    from random import seed
    from pytools.stopwatch import Job

    from pyrticle.units import SI
    units = SI()

    seed(0)

    # user interface ----------------------------------------------------------
    def make_setup():
        from pyrticle.reconstruction import \
                ShapeFunctionReconstructor, \
                NormalizedShapeFunctionReconstructor, \
                AdvectiveReconstructor
        from pyrticle.pusher import \
                MonomialParticlePusher, \
                AverageParticlePusher

        variables = {
                "radial_subdiv": 15,
                "tube_length": 0.1*units.M,
                "tube_radius": 25*units.MM,
                "tube_radius_inner": 2.5*units.MM,
                "max_volume_inner": 10*units.MM**3,
                "max_volume_outer": 100*units.MM**3,

                "element_order": 3,
                "shape_exponent": 2,

                "chi": None,
                "phi_decay": 0,

                "final_time": 0.1*units.M/units.VACUUM_LIGHT_SPEED,

                "pusher": None,
                "reconstructor": None,

                "nparticles": 20000,
                "cloud_charge": -10e-9 * units.C,
                "beam_emittance": 5*units.MM*units.MRAD,
                "beam_radius": 2.5*units.MM,
                "beam_length": 5*units.MM,

                "vis_interval": 100,
                }

        from pyrticle.reconstruction import \
                ShapeFunctionReconstructor, \
                NormalizedShapeFunctionReconstructor, \
                AdvectiveReconstructor
        from pyrticle.pusher import \
                MonomialParticlePusher, \
                AverageParticlePusher

        constants = {
                "num": num,
                "comp": comp,
                "units": units,

                "RecShape": ShapeFunctionReconstructor,
                "RecNormShape": NormalizedShapeFunctionReconstructor,
                "RecAdv": AdvectiveReconstructor,

                "PushMonomial": MonomialParticlePusher,
                "PushAverage": AverageParticlePusher,
                }

        doc = {
                "chi": "relative speed of hyp. cleaning (None for no cleaning)",
                "tube_length": "how long a beam tube [m]",
                "nparticles": "how many particles",
                "beam_radius": "total radius of the beam [m]",
                "beam_length": "Z-wise length of the beam [m]",
                "beam_emittance": "total emittance of the beam [m*rad]",
                "vis_interval": "how often a visualization of the fields is written",
                "max_volume_inner": "max. tet volume in inner mesh [m^3]",
                "max_volume_outer": "max. tet volume in outer mesh [m^3]",
                }

        from pytools import gather_parameters_from_user
        return gather_parameters_from_user(variables, constants, doc)

    setup = make_setup()

    # discretization setup ----------------------------------------------------
    job = Job("mesh")
    #full_mesh = make_cylinder_mesh(radius=25*units.MM, height=25*units.MM, periodic=True,
            #max_volume=1000*units.MM**3, radial_subdivisions=10)
    #full_mesh = make_box_mesh([1,1,2], max_volume=0.01)

    if True:
        from tubemesh import make_cylinder_with_fine_core
        full_mesh = make_cylinder_with_fine_core(
                r=setup.tube_radius, inner_r=setup.tube_radius_inner, 
                min_z=0, max_z=setup.tube_length,
                max_volume_inner=setup.max_volume_inner,
                max_volume_outer=setup.max_volume_outer,
                radial_subdiv=setup.radial_subdiv,
                )
    if False:
        # pillbox cavity
        from tubemesh import make_extrusion_with_fine_core
        full_mesh = make_extrusion_with_fine_core(
                rz=[
                    (1*setup.tube_radius,0),
                    (1*setup.tube_radius,setup.tube_length*0.333),
                    (2*setup.tube_radius,setup.tube_length*0.333),
                    (2*setup.tube_radius,setup.tube_length*0.666),
                    (1*setup.tube_radius,setup.tube_length*0.666),
                    (1*setup.tube_radius,setup.tube_length),
                    ],
                inner_r=setup.tube_radius_inner, 
                max_volume_inner=setup.max_volume_inner,
                max_volume_outer=setup.max_volume_outer,
                radial_subdiv=setup.radial_subdiv,
                )
    job.done()

    from hedge.parallel import guess_parallelization_context

    pcon = guess_parallelization_context()

    if pcon.is_head_rank:
        mesh = pcon.distribute_mesh(full_mesh)
    else:
        mesh = pcon.receive_mesh()

    job = Job("discretization")
    discr = pcon.make_discretization(mesh, TetrahedralElement(setup.element_order))
    job.done()

    vis = SiloVisualizer(discr)
    #vis = VtkVisualizer(discr, "pic")

    max_op = MaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)

    if setup.chi is not None:
        from pyrticle.hyperbolic import ECleaningMaxwellOperator
        max_op = ECleaningMaxwellOperator(max_op, 
                chi=setup.chi, 
                phi_decay=setup.phi_decay)

    div_op = DivergenceOperator(discr)

    dt = discr.dt_factor(max_op.c) / 2
    nsteps = int(setup.final_time/dt)+1
    dt = setup.final_time/nsteps

    print "#elements=%d, dt=%s, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    def l2_norm(field):
        return sqrt(dot(field, discr.mass_operator*field))
    def l2_error(field, true):
        return l2_norm(field-true)/l2_norm(true)

    # particles setup ---------------------------------------------------------
    from pyrticle.cloud import ParticleCloud

    from pyrticle.reconstruction import Reconstructor
    from pyrticle.pusher import Pusher

    assert isinstance(setup.reconstructor, Reconstructor), \
            "must specify valid reconstructor"
    assert isinstance(setup.pusher, Pusher), \
            "must specify valid reconstructor"

    cloud = ParticleCloud(discr, units, 
            setup.reconstructor, setup.pusher,
            dimensions_pos=3, dimensions_velocity=3, 
            verbose_vis=True)

    electrons_per_particle = abs(setup.cloud_charge/setup.nparticles/units.EL_CHARGE)
    print "e-/particle = ", electrons_per_particle 

    el_energy = units.EL_REST_ENERGY*10
    el_lorentz_gamma = el_energy/units.EL_REST_ENERGY
    mean_beta = (1-1/el_lorentz_gamma**2)**0.5
    gamma = 1/sqrt(1-mean_beta**2)
    print "beta = %g, gamma = %g" % (mean_beta, gamma)

    beam = KVZIntervalBeam(units, setup.nparticles, 
            p_charge=setup.cloud_charge/setup.nparticles, 
            p_mass=electrons_per_particle*units.EL_MASS,
            radii=2*[setup.beam_radius],
            emittances=2*[setup.beam_emittance], 
            z_length=setup.beam_length,
            z_pos=10*units.MM,
            beta=mean_beta)
    beam.add_to(cloud, discr)

    from pyrticle.cloud import optimize_shape_bandwidth, guess_shape_bandwidth
    optimize_shape_bandwidth(cloud, discr, beam.analytic_rho(discr), 
            exponent=setup.shape_exponent, plot_l1_errors=True)
    #guess_shape_bandwidth(cloud)

    # initial condition -------------------------------------------------------
    from pyrticle.cloud import compute_initial_condition
    job = Job("initial condition")
    fields = compute_initial_condition(pcon, discr, cloud, 
            mean_beta=num.array([0, 0, mean_beta]), max_op=max_op, debug=True,
            force_zero=False)
    job.done()

    # timestepping setup ------------------------------------------------------
    stepper = RK4TimeStepper()
    t = 0

    # diagnostics setup -------------------------------------------------------
    from pytools.log import LogManager, \
            add_simulation_quantities, \
            add_general_quantities, \
            add_run_info, ETA
    from pyrticle.log import add_particle_quantities, add_field_quantities, \
            add_beam_quantities, add_currents
    logmgr = LogManager("kv.dat")
    add_run_info(logmgr)
    add_general_quantities(logmgr)
    add_simulation_quantities(logmgr, dt)
    add_particle_quantities(logmgr, cloud)
    add_field_quantities(logmgr, fields)
    add_beam_quantities(logmgr, cloud, axis=0, beam_axis=2)
    add_currents(logmgr, fields, (0,0,1), setup.tube_length)

    stepper.add_instrumentation(logmgr)
    fields.add_instrumentation(logmgr)
    logmgr.set_constant("beta", mean_beta)
    logmgr.set_constant("gamma", gamma)
    logmgr.set_constant("vz", mean_beta*units.VACUUM_LIGHT_SPEED)
    logmgr.set_constant("Q0", setup.cloud_charge)
    logmgr.set_constant("n_part_0", setup.nparticles)
    logmgr.set_constant("pmass", electrons_per_particle*units.EL_MASS)

    from pytools.log import IntervalTimer
    vis_timer = IntervalTimer("t_vis", "Time spent visualizing")
    logmgr.add_quantity(vis_timer)

    logmgr.add_quantity(ETA(nsteps))

    logmgr.add_watches(["step", "t_sim", "W_field", "t_step", "t_eta", "n_part"])

    from kv import KVPredictedRadius
    logmgr.add_quantity(KVPredictedRadius(dt, 
        beam_v=mean_beta*units.VACUUM_LIGHT_SPEED,
        predictor=beam.get_rms_predictor(axis=0),
        suffix="x_rms"))
    logmgr.add_quantity(KVPredictedRadius(dt, 
        beam_v=mean_beta*units.VACUUM_LIGHT_SPEED,
        predictor=beam.get_total_predictor(axis=0),
        suffix="x_total"))

    # timestepping ------------------------------------------------------------
    for step in xrange(nsteps):
        logmgr.tick()

        cloud.upkeep()
        fields = stepper(fields, t, dt, fields.rhs)

        if step % setup.vis_interval == 0:
            vis_timer.start()
            visf = vis.make_file("pic-%04d" % step)

            cloud.add_to_vis(vis, visf, time=t, step=step)
            vis.add_data(visf, [
                ("divD", max_op.epsilon*div_op(fields.e)),
                ("e", fields.e), 
                ("h", fields.h), 
                ("j", cloud.reconstruct_j()), 
                ],
                time=t, step=step)
            visf.close()
            vis_timer.stop()

        t += dt

    vis.close()
        
    logmgr.tick()
    logmgr.save()

    _, _, err_table = logmgr.get_expr_dataset("(rx_rms-rx_rms_theory)/rx_rms_theory")
    print "Relative error (rms): %g" % max(err for step, err in err_table)




if __name__ == "__main__":
    #import cProfile as profile
    #profile.run("main()", "pic.prof")
    main()
