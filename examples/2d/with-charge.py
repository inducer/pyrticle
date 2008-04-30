from __future__ import division
import numpy
import numpy.linalg as la
import cProfile as profile
import pytools




class GaussianParticleDistribution(pytools.Record):
    def __init__(self, total_charge, total_mass, mean_x, mean_p, sigma_x, sigma_p):
        pytools.Record.__init__(self, locals())

    def add_to(self, cloud, nparticles):
        from random import gauss

        pmass = self.total_mass/nparticles
        cloud.add_particles(
                positions=[
                    numpy.array([gauss(m, s) for m, s in zip(self.mean_x, self.sigma_x)]) 
                    for i in range(nparticles)
                    ],
                velocities=[cloud.units.v_from_p(pmass, 
                    numpy.array([gauss(m, s) for m, s in zip(self.mean_p, self.sigma_p)])) 
                    for i in range(nparticles)
                    ],
                charges=self.total_charge/nparticles, 
                masses=pmass)

    def analytic_rho(self, discr):
        from math import exp, pi

        sigma_mat = numpy.diag(self.sigma_x**2)
        inv_sigma_mat = numpy.diag(self.sigma_x**(-2))

        from numpy import dot

        def distrib(x):
            x0 = x-self.mean_x
            return 1/((2*pi)**(len(x)/2) * la.det(sigma_mat)**0.5) \
                    * exp(-0.5*dot(x0, dot(inv_sigma_mat, x0)))

        rho = self.total_charge * discr.interpolate_volume_function(distrib)

        # check for correctness
        from hedge.discretization import integral
        int_rho = integral(discr, rho)
        rel_err = (int_rho-self.total_charge)/self.total_charge
        assert rel_err < 1e-2

        return rho




def main():
    from hedge.element import TriangularElement
    from hedge.timestep import RK4TimeStepper
    from hedge.visualization import VtkVisualizer, SiloVisualizer
    from math import sqrt, pi
    from pytools.arithmetic_container import join_fields

    from random import seed
    seed(0)

    from pyrticle.units import SI
    units = SI()

    from pytools.log import LogManager
    logmgr = LogManager("2d.dat", "w")

    # user interface ----------------------------------------------------------
    def make_setup():
        c0 = units.VACUUM_LIGHT_SPEED

        variables = {
                "mesh": None,
                "tube_length": 2,
                "tube_width": 1,
                "tube_periodic": True,
                "tube_max_tri_area": 0.02,

                "element_order": 7,
                "shape_exponent": 2,
                "shape_bandwidth": "optimize",

                "chi": None,
                "phi_decay": 0,
                "phi_filter": None,

                "final_time": 10*units.M/units.VACUUM_LIGHT_SPEED,

                "sigma_x": 0.1*numpy.ones((2,)),
                "mean_x": numpy.zeros((2,)),
                "mean_v": numpy.array([c0*0.9, 0]),
                "sigma_v": numpy.array([c0*0.9*1e-3, c0*0.9*1e-6]),
                "nparticles": 1000,
                "cloud_charge": -1e-9 * units.C,

                "vis_interval": 100,
                "vis_path": ".",
                }
        
        from hedge.mesh import make_rect_mesh

        constants = {
                "units": units,
                "make_rect_mesh": make_rect_mesh,
                }

        doc = {
                "chi": "relative speed of hyp. cleaning (None for no cleaning)",
                "shape_bandwidth": "either 'optimize', 'guess' or a positive real number",
                "phi_filter": "a tuple (min_amp, order) or None, describing the filtering applied to phi in hypclean mode",
                }

        from pyrticle.tools import PICCPyUserInterface
        ui = PICCPyUserInterface(variables, constants, doc)
        return ui.gather()

    setup = make_setup()

    # discretization setup ----------------------------------------------------
    if setup.mesh is None:
        from hedge.mesh import make_rect_mesh
        setup.mesh = make_rect_mesh(
                a=(-0.5, -setup.tube_width/2),
                b=(-0.5+setup.tube_length, setup.tube_width/2),
                periodicity=(setup.tube_periodic, False),
                subdivisions=(10,5),
                max_area=setup.tube_max_tri_area)
    else:
        from hedge.mesh import Mesh
        assert isinstance(setup.mesh, Mesh)

    from hedge.parallel import guess_parallelization_context

    pcon = guess_parallelization_context()

    if pcon.is_head_rank:
        mesh = pcon.distribute_mesh(setup.mesh)
    else:
        mesh = pcon.receive_mesh()

    discr = pcon.make_discretization(mesh, TriangularElement(setup.element_order))
    fine_discr = pcon.make_discretization(mesh, TriangularElement(3*setup.element_order))
    vis = SiloVisualizer(discr)
    fine_vis = SiloVisualizer(fine_discr)
    #vis = VtkVisualizer(discr, "pic")

    from hedge.discretization import Projector
    coarse_to_fine = Projector(discr, fine_discr)

    from hedge.operators import TEMaxwellOperator, DivergenceOperator
    from hedge.mesh import TAG_ALL, TAG_NONE

    max_op = TEMaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)

    em_filters = []
    if setup.chi is not None:
        from pyrticle.hyperbolic import ECleaningMaxwellOperator
        max_op = ECleaningMaxwellOperator(max_op, 
                chi=setup.chi, 
                phi_decay=setup.phi_decay)

        if setup.phi_filter is not None:
            from pyrticle.hyperbolic import PhiFilter
            from hedge.discretization import Filter, ExponentialFilterResponseFunction
            em_filters.append(
                    PhiFilter(max_op, Filter(discr,
                        ExponentialFilterResponseFunction(*setup.phi_filter))))

    div_op = DivergenceOperator(discr)

    dt = discr.dt_factor(max_op.max_eigenvalue())
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
    cloud = ParticleCloud(discr, units, 
            setup.reconstructor, setup.pusher, setup.finder,
            dimensions_pos=2, dimensions_velocity=2,
            verbose_vis=True)

    electrons_per_particle = abs(setup.cloud_charge/setup.nparticles/units.EL_CHARGE)
    print "e-/particle = ", electrons_per_particle 

    mean_beta = setup.mean_v/units.VACUUM_LIGHT_SPEED
    gamma = units.gamma(setup.mean_v)
    pmass = electrons_per_particle*units.EL_MASS
    mean_p = gamma*pmass*setup.mean_v

    print "beta=%g, gamma=%g" % (la.norm(mean_beta), gamma)

    gauss_p = GaussianParticleDistribution(
            total_charge=setup.cloud_charge, 
            total_mass=pmass*setup.nparticles,
            mean_x=setup.mean_x,
            mean_p=mean_p,
            sigma_x=setup.sigma_x,
            sigma_p=gamma*pmass*setup.sigma_v)
    gauss_p.add_to(cloud, setup.nparticles)

    from pyrticle.cloud import set_shape_bandwidth
    set_shape_bandwidth(cloud, setup.shape_bandwidth, setup.shape_exponent,
            gauss_p.analytic_rho(discr))

    # intial condition --------------------------------------------------------
    from pyrticle.cloud import compute_initial_condition
    fields = compute_initial_condition(pcon, discr, cloud, mean_beta, max_op,
            debug=True)
    for f in em_filters:
        fields.add_em_filter(f)

    stepper = RK4TimeStepper()

    # diagnostics setup -------------------------------------------------------
    from pytools.log import \
            add_simulation_quantities, \
            add_general_quantities, \
            add_run_info, ETA
    from pyrticle.log import add_particle_quantities, add_field_quantities, \
            add_beam_quantities, add_currents
    add_run_info(logmgr)
    add_general_quantities(logmgr)
    add_simulation_quantities(logmgr, dt)
    add_particle_quantities(logmgr, cloud)
    add_field_quantities(logmgr, fields, reconstruct_interval=1)
    add_beam_quantities(logmgr, cloud, axis=1, beam_axis=0)
    add_currents(logmgr, fields, (1,0), setup.tube_length)

    stepper.add_instrumentation(logmgr)
    fields.add_instrumentation(logmgr)
    logmgr.set_constant("beta", la.norm(mean_beta))
    logmgr.set_constant("gamma", gamma)
    logmgr.set_constant("mean_v", setup.mean_v)
    logmgr.set_constant("Q0", setup.cloud_charge)
    logmgr.set_constant("n_part_0", setup.nparticles)
    logmgr.set_constant("pmass", electrons_per_particle*units.EL_MASS)
    logmgr.set_constant("chi", setup.chi)
    logmgr.set_constant("shape_radius_setup", setup.shape_bandwidth)
    logmgr.set_constant("shape_radius", cloud.reconstructor.shape_function.radius)
    logmgr.set_constant("shape_exponent", cloud.reconstructor.shape_function.exponent)

    from pytools.log import IntervalTimer
    vis_timer = IntervalTimer("t_vis", "Time spent visualizing")
    logmgr.add_quantity(vis_timer)

    logmgr.add_quantity(ETA(nsteps))

    logmgr.add_watches(["step", "t_sim", "W_field", "t_step", "t_eta", "n_part"])

    # timestepping ------------------------------------------------------------
    t = 0
    rec = cloud.reconstructor

    for step in xrange(nsteps):
        logmgr.tick()

        if step % setup.vis_interval == 0:
            vis_timer.start()
            import os.path
            visf = vis.make_file(os.path.join(
                setup.vis_path, "pic-%04d" % step))

            cloud.add_to_vis(vis, visf, time=t, step=step, beamaxis=0)
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
            try:
                rec.visualize_grid_quantities
            except AttributeError:
                pass
            else:
                rec.visualize_grid_quantities(visf, [
                        ("rho_grid", rec.reconstruct_grid_rho()),
                        ("j_grid", rec.reconstruct_grid_j(cloud.velocities())),
                        ("ones_resid", rec.remap_residual(rec.ones_on_grid())),
                        ("usecount", rec.grid_usecount()),
                        ])

            visf.close()
            vis_timer.stop()

        fields.upkeep()
        fields = stepper(fields, t, dt, fields.rhs)

        t += dt

    vis.close()

    logmgr.tick()
    logmgr.save()




if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
