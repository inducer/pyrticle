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

    beam_radius = 2.5*units.MM

    # discretization setup ----------------------------------------------------
    job = Job("mesh")
    #full_mesh = make_cylinder_mesh(radius=25*units.MM, height=25*units.MM, periodic=True,
            #max_volume=1000*units.MM**3, radial_subdivisions=10)
    #full_mesh = make_box_mesh([1,1,2], max_volume=0.01)

    if True:
        full_mesh = make_cylinder_with_fine_core(
                r=10*beam_radius, inner_r=1*beam_radius, 
                min_z=0, max_z=20*beam_radius,
                max_volume_inner=10*units.MM**3,
                max_volume_outer=100*units.MM**3,
                radial_subdiv=10,
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
    final_time = 0.3*units.M/max_op.c
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

    nparticles = 3000
    cloud_charge = 1e-9 * units.C
    electrons_per_particle = cloud_charge/nparticles/units.EL_CHARGE
    print "e-/particle = ", electrons_per_particle 

    emittance = 5 * units.MM * units.MRAD
    initial_radius = 2.5*units.MM

    el_energy = 5.2e6 * units.EV
    #el_energy = units.EL_REST_ENERGY*1.00001
    el_lorentz_gamma = el_energy/units.EL_REST_ENERGY
    #el_lorentz_gamma = 100000
    beta = (1-1/el_lorentz_gamma**2)**0.5
    gamma = 1/sqrt(1-beta**2)
    print "beta = %g, gamma = %g" % (beta, gamma)

    beam = KVZIntervalBeam(units, nparticles, 
            p_charge=cloud_charge/nparticles, 
            p_mass=electrons_per_particle*units.EL_MASS,
            radii=2*[beam_radius],
            emittances=2*[5 * units.MM * units.MRAD], 
            z_length=5*units.MM,
            z_pos=10*units.MM,
            beta=beta)
    beam.add_to(cloud, discr)

    # intial condition --------------------------------------------------------
    def compute_initial_condition():
        from hedge.operators import WeakPoissonOperator
        from hedge.mesh import TAG_ALL, TAG_NONE
        from hedge.data import ConstantGivenFunction, GivenVolumeInterpolant
        from hedge.tools import cross

        # see doc/notes.tm for derivation of IC

        beta_vec = num.array([0,0,beta])

        diff_tensor = num.identity(discr.dimensions)
        diff_tensor[2,2] = 1/gamma**2

        poisson_op = WeakPoissonOperator(discr, 
                diffusion_tensor=ConstantGivenFunction(diff_tensor),
                dirichlet_tag=TAG_ALL,
                neumann_tag=TAG_NONE,
                )

        rho = cloud.reconstruct_rho() 

        from hedge.discretization import ones_on_volume
        print "charge: supposed=%g reconstructed=%g" % (
                cloud_charge,
                ones_on_volume(discr)*(discr.mass_operator*rho),
                )

        from hedge.tools import parallel_cg
        phi = -parallel_cg(pcon, -poisson_op, 
                poisson_op.prepare_rhs(
                    GivenVolumeInterpolant(discr, rho/max_op.epsilon)), 
                debug=True, tol=1e-10)

        etilde = ArithmeticList([1,1,1/gamma])*poisson_op.grad(phi)

        eprime = ArithmeticList([gamma,gamma,1])*etilde

        hprime = (1/max_op.mu)*gamma/max_op.c * cross(beta_vec, etilde)

        rhoprime = gamma*rho
        divDprime_ldg = max_op.epsilon*poisson_op.div(eprime)
        divDprime_ldg2 = max_op.epsilon*poisson_op.div(eprime, gamma*phi)
        divDprime_ldg3 = max_op.epsilon*gamma*\
                (discr.inverse_mass_operator*poisson_op.op(phi))
        divDprime_central = max_op.epsilon*div_op(eprime)

        print "l2 div error ldg: %g" % \
                l2_error(divDprime_ldg, rhoprime)
        print "l2 div error central: %g" % \
                l2_error(divDprime_central, rhoprime)
        print "l2 div error ldg with phi: %g" % \
                l2_error(divDprime_ldg2, rhoprime)
        print "l2 div error ldg with phi 3: %g" % \
                l2_error(divDprime_ldg3, rhoprime)

        if True:
            visf = vis.make_file("ic")
            vis.add_data(visf, [ 
                ("rho", rhoprime), 
                ("divDldg", divDprime_ldg),
                ("divDldg2", divDprime_ldg2),
                ("divDldg3", divDprime_ldg3),
                ("divDcentral", divDprime_central),
                ("phi", phi),
                ("e", eprime), 
                ("h", hprime), 
                ],
                scale_factor=1e30
                )
            cloud.add_to_vis(vis, visf)
            visf.close()

        from pyrticle.cloud import FieldsAndCloud
        return FieldsAndCloud(max_op, eprime, hprime, cloud)

    fields = compute_initial_condition()

    # timestepping ------------------------------------------------------------

    def rhs(t, y):
        rho, j = cloud.reconstruct_densities()

        e, h = max_op.split_fields(y)
        cloud_rhs = cloud.rhs(t, e, h)

        maxwell_rhs = max_op.rhs(t, y[0:eh_components])
        rhs_e, rhs_h = max_op.split_fields(maxwell_rhs)

        return join_fields(
                rhs_e - 1/max_op.epsilon*j,
                rhs_h,
                ).plus([cloud_rhs])

    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    rms_r_logger = RMSBeamRadiusLogger(cloud.dimensions_pos, 0)

    rms_theory_with_charge = KVRadiusPredictor(
            beam.rms_radii[rms_r_logger.axis], 
            beam.rms_emittances[rms_r_logger.axis],
            xi=beam.get_space_charge_parameter())

    def write_out_plots():
        rms_r_logger.generate_plot(
                title="Kapchinskij-Vladimirskij Beam Evolution",
                sim_label="RMS, simulated, with space charge", 
                outfile="beam-rad-rms.eps",
                theories=[
                    ("RMS, theoretical, with space charge", rms_theory_with_charge)
                    ])

    from pytools.stopwatch import EtaEstimator
    eta = EtaEstimator(nsteps)

    for step in xrange(nsteps):
        if True:
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
