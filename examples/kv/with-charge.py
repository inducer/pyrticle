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
    from pyrticle.cloud import ParticleCloud
    from kv import \
            KVZIntervalBeam, \
            KVRadiusPredictor, \
            MaxBeamRadiusLogger, \
            RMSBeamRadiusLogger
    from tubemesh import make_cylinder_with_fine_core
    from random import seed
    from pytools.stopwatch import Job

    from pyrticle.units import SI
    units = SI()

    seed(0)

    # parse command line ------------------------------------------------------
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option(
            "--radial-subdiv", dest="radial_subdiv", default="10",
            help="how many angular subdivisions in the surface mesh")
    parser.add_option(
            "--tube-length", dest="tube_length", default="0.1",
            help="how long a beam tube [m]")
    parser.add_option(
            "--nparticles", dest="nparticles", default="1000",
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

    full_mesh = make_cylinder_with_fine_core(
            r=10*beam_radius, inner_r=1*beam_radius, 
            min_z=0, max_z=tube_length,
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
    discr = pcon.make_discretization(mesh, TetrahedralElement(1))
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
    cloud = ParticleCloud(discr, units, 3, 3, verbose_vis=True)

    cloud_charge = -1e-9 * units.C
    electrons_per_particle = abs(cloud_charge/nparticles/units.EL_CHARGE)
    print "e-/particle = ", electrons_per_particle 

    el_energy = 5.2e6 * units.EV
    #el_energy = units.EL_REST_ENERGY*1.00001
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
    from kv import calculate_rms_energy_spread

    print "energy spread: %g %%" % (
            calculate_rms_energy_spread(cloud) * 100)

    # initial condition -------------------------------------------------------
    from pyrticle.cloud import compute_initial_condition
    fields = compute_initial_condition(pcon, discr, cloud, 
            mean_beta=num.array([0, 0, mean_beta]), max_op=max_op, debug=True)

    # timestepping ------------------------------------------------------------
    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    rms_r_logger = RMSBeamRadiusLogger(cloud.dimensions_pos, 0)

    rms_theory_with_charge = KVRadiusPredictor(
            beam.rms_radii[rms_r_logger.axis], 
            beam.rms_emittances[rms_r_logger.axis],
            xi=beam.get_rms_space_charge_parameter())

    def write_out_plots():
        rms_r_logger.generate_plot(
                title="Kapchinskij-Vladimirskij Beam Evolution",
                sim_label="RMS, simulated, with space charge", 
                outfile="beam-rad-rms",
                theories=[
                    ("RMS, theoretical, with space charge", rms_theory_with_charge)
                    ])

    from pytools.stopwatch import EtaEstimator
    eta = EtaEstimator(nsteps)

    for step in xrange(nsteps):
        if step % field_dump_interval == 0:
            visf = vis.make_file("pic-%04d" % step)

            mesh_scalars, mesh_vectors = \
                    cloud.add_to_vis(vis, visf, time=t, step=step)
            vis.add_data(visf, [
                ("divD", max_op.epsilon*div_op(fields.e)),
                ("e", fields.e), 
                ("h", fields.h), 
                ]
                + mesh_vectors
                + mesh_scalars,
                time=t, step=step)
            visf.close()

        rms_r_logger.update(t, cloud.positions, cloud.velocities())

        fields = stepper(fields, t, dt, fields.rhs)
        cloud.upkeep()

        print "timestep %d, t=%g l2[e]=%g l2[h]=%g secs=%f eta=%s particles=%d" % (
                step, t, l2_norm(fields.e), l2_norm(fields.h),
                time()-last_tstep, eta.estimate(step), len(cloud))
        last_tstep = time()

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

        t += dt

        if step % 100 == 0:
            write_out_plots()

    vis.close()
        
    write_out_plots()

    print "Relative error: %g" % rms_r_logger.relative_error(rms_theory_with_charge)





if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
