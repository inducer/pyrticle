from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pylinear.operator as op
import cProfile as profile




def main():
    from hedge.element import TetrahedralElement
    from hedge.timestep import RK4TimeStepper
    from hedge.mesh import \
            make_box_mesh, \
            make_cylinder_mesh
    from hedge.discretization import \
            Discretization, \
            pair_with_boundary
    from hedge.visualization import VtkVisualizer, SiloVisualizer
    from hedge.tools import dot
    from math import sqrt, pi
    from pytools.arithmetic_container import \
            ArithmeticList, join_fields
    from hedge.operators import MaxwellOperator, DivergenceOperator
    from kv import KVZIntervalBeam
    from random import seed
    from pytools.stopwatch import Job

    from pyrticle.units import SI
    units = SI()

    seed(0)

    # parse command line ------------------------------------------------------
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option(
            "--radial-subdiv", dest="radial_subdiv", default="15",
            help="how many angular subdivisions in the surface mesh")
    parser.add_option(
            "--tube-length", dest="tube_length", default="0.1",
            help="how long a beam tube [m]")
    parser.add_option(
            "--nparticles", dest="nparticles", default="20000",
            help="how many particles")
    parser.add_option(
            "--beam-radius", dest="beam_radius", default="2.5",
            help="radius of the beam [mm]")
    parser.add_option(
            "--emittance", dest="emittance", default="5",
            help="total emittance of the beam [mm*mrad]")
    parser.add_option(
            "--final-time", dest="final_time", default="0.1",
            help="how long to run the computation [m]")
    parser.add_option(
            "--field-dump-interval", dest="field_dump_interval", default="100",
            help="every how many time steps to dump the fields")

    options, args = parser.parse_args()

    radial_subdiv = int(options.radial_subdiv)
    tube_length = float(options.tube_length)*units.M
    nparticles = int(options.nparticles)
    beam_radius = float(options.beam_radius)*units.MM
    emittance = float(options.emittance) * units.MM * units.MRAD
    final_time = float(options.final_time)*units.M/units.VACUUM_LIGHT_SPEED
    field_dump_interval = int(options.field_dump_interval)

    # discretization setup ----------------------------------------------------
    job = Job("mesh")
    #full_mesh = make_cylinder_mesh(radius=25*units.MM, height=25*units.MM, periodic=True,
            #max_volume=1000*units.MM**3, radial_subdivisions=10)
    #full_mesh = make_box_mesh([1,1,2], max_volume=0.01)

    if False:
        from tubemesh import make_cylinder_with_fine_core
        full_mesh = make_cylinder_with_fine_core(
                r=10*beam_radius, inner_r=1*beam_radius, 
                min_z=0, max_z=tube_length,
                max_volume_inner=10*units.MM**3,
                max_volume_outer=100*units.MM**3,
                radial_subdiv=radial_subdiv,
                )
    if True:
        # pillbox cavity
        from tubemesh import make_extrusion_with_fine_core
        full_mesh = make_extrusion_with_fine_core(
                rz=[
                    (10*beam_radius,0),
                    (10*beam_radius,tube_length*0.333),
                    (20*beam_radius,tube_length*0.333),
                    (20*beam_radius,tube_length*0.666),
                    (10*beam_radius,tube_length*0.666),
                    (10*beam_radius,tube_length),
                    ],
                inner_r=1*beam_radius, 
                max_volume_inner=10*units.MM**3,
                max_volume_outer=100*units.MM**3,
                radial_subdiv=radial_subdiv,
                )
    job.done()

    from hedge.parallel import guess_parallelization_context

    pcon = guess_parallelization_context()

    if pcon.is_head_rank:
        mesh = pcon.distribute_mesh(full_mesh)
    else:
        mesh = pcon.receive_mesh()

    job = Job("discretization")
    discr = pcon.make_discretization(mesh, TetrahedralElement(3))
    job.done()

    vis = SiloVisualizer(discr)
    #vis = VtkVisualizer(discr, "pic")

    max_op = MaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)
    div_op = DivergenceOperator(discr)

    dt = discr.dt_factor(max_op.c) / 2
    nsteps = int(final_time/dt)+1
    dt = final_time/nsteps

    print "#elements=%d, dt=%s, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    def l2_norm(field):
        return sqrt(dot(field, discr.mass_operator*field))
    def l2_error(field, true):
        return l2_norm(field-true)/l2_norm(true)

    # particles setup ---------------------------------------------------------
    from pyrticle.cloud import ParticleCloud
    from pyrticle.reconstruction import ShapeFunctionReconstructor
    from pyrticle.pusher import MonomialParticlePusher
    cloud = ParticleCloud(discr, units, 
            ShapeFunctionReconstructor(),
            MonomialParticlePusher(),
            3, 3, verbose_vis=True)

    cloud_charge = -10e-9 * units.C
    electrons_per_particle = abs(cloud_charge/nparticles/units.EL_CHARGE)
    print "e-/particle = ", electrons_per_particle 

    el_energy = units.EL_REST_ENERGY*10
    el_lorentz_gamma = el_energy/units.EL_REST_ENERGY
    #el_lorentz_gamma = 100000
    mean_beta = (1-1/el_lorentz_gamma**2)**0.5
    gamma = 1/sqrt(1-mean_beta**2)
    print "beta = %g, gamma = %g" % (mean_beta, gamma)

    beam = KVZIntervalBeam(units, nparticles, 
            p_charge=cloud_charge/nparticles, 
            p_mass=electrons_per_particle*units.EL_MASS,
            radii=2*[beam_radius],
            emittances=2*[5 * units.MM * units.MRAD], 
            z_length=5*units.MM,
            z_pos=10*units.MM,
            beta=mean_beta)
    beam.add_to(cloud, discr)

    # initial condition -------------------------------------------------------
    from pyrticle.cloud import compute_initial_condition
    job = Job("initial condition")
    fields = compute_initial_condition(pcon, discr, cloud, 
            mean_beta=num.array([0, 0, mean_beta]), max_op=max_op, debug=True)
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
    add_currents(logmgr, fields, (0,0,1), tube_length)

    stepper.add_instrumentation(logmgr)
    fields.add_instrumentation(logmgr)
    logmgr.set_constant("beta", mean_beta)
    logmgr.set_constant("gamma", gamma)
    logmgr.set_constant("vz", units.VACUUM_LIGHT_SPEED)
    logmgr.set_constant("Q0", cloud_charge)
    logmgr.set_constant("n_part_0", nparticles)
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

        if step % field_dump_interval == 0:
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
    #profile.run("main()", "pic.prof")
    main()
