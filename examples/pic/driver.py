from __future__ import division
import numpy
import numpy.linalg as la
import pytools




class PICCPyUserInterface(pytools.CPyUserInterface):
    def __init__(self, units):
        from pyrticle.reconstruction import \
                ShapeFunctionReconstructor, \
                NormalizedShapeFunctionReconstructor, \
                AdvectiveReconstructor, \
                GridReconstructor, \
                GridFindReconstructor, \
                SingleBrickGenerator, \
                FineCoreBrickGenerator

        from pyrticle.pusher import \
                MonomialParticlePusher, \
                AverageParticlePusher

        from pyrticle.cloud import \
                FaceBasedElementFinder, \
                HeuristicElementFinder

        import pyrticle.geometry
        import pyrticle.distribution

        constants = {
                "numpy": numpy,
                "la": numpy.linalg,

                "units": units,
                "pyrticle": pyrticle,

                "RecShape": ShapeFunctionReconstructor,
                "RecNormShape": NormalizedShapeFunctionReconstructor,
                "RecAdv": AdvectiveReconstructor,
                "RecGrid": GridReconstructor,
                "RecGridFind": GridFindReconstructor,

                "SingleBrickGenerator": SingleBrickGenerator,
                "FineCoreBrickGenerator": FineCoreBrickGenerator,

                "PushMonomial": MonomialParticlePusher,
                "PushAverage": AverageParticlePusher,

                "FindHeuristic": HeuristicElementFinder,
                "FindFaceBased": FaceBasedElementFinder,
                }

        import hedge.data

        variables = {
                "pusher": None,
                "reconstructor": None,
                "finder": FaceBasedElementFinder(),

                "mesh": None,
                "dimensions_pos": None,
                "dimensions_velocity": None,

                "beam_axis": None,
                "beam_diag_axis": None,
                "tube_length": None,

                "element_order": None,
                "maxwell_flux_type": "lf",
                "maxwell_bdry_flux_type": 1,

                "shape_exponent": 2,
                "shape_bandwidth": "optimize",

                "chi": None,
                "phi_decay": 0,
                "phi_filter": None,

                "potential_bc": hedge.data.ConstantGivenFunction(),

                "final_time": None,

                "nparticles": 20000,
                "distribution": None,

                "vis_interval": 100,
                "vis_pattern": "pic-%04d",
                "vis_order": None,
                "output_path": ".",

                "debug": set(["ic", "poisson", "shape_bw"]),
                "dg_debug": set(),

                "watch_vars": ["step", "t_sim", "W_field", "t_step", "t_eta", "n_part"],

                "hook_startup": lambda runner : None,
                "hook_before_step": lambda runner : None,
                "hook_after_step": lambda runner : None,
                "hook_when_done": lambda runner : None,
                "hook_vis_quantities": lambda runner: [
                    ("e", runner.fields.e), 
                    ("h", runner.fields.h), 
                    ("j", runner.cloud.reconstruct_j()), 
                    ],
                "hook_visualize": lambda runner, vis, visf: None,
                }

        doc = {
                "chi": "relative speed of hyp. cleaning (None for no cleaning)",
                "nparticles": "how many particles",
                "vis_interval": "how often a visualization of the fields is written",
                "max_volume_inner": "max. tet volume in inner mesh [m^3]",
                "max_volume_outer": "max. tet volume in outer mesh [m^3]",
                "shape_bandwidth": "either 'optimize', 'guess' or a positive real number",
                "phi_filter": "a tuple (min_amp, order) or None, describing the filtering applied to phi in hypclean mode",
                }

        pytools.CPyUserInterface.__init__(self, variables, constants, doc)
    
    def validate(self, setup):
        pytools.CPyUserInterface.validate(self, setup)

        from pyrticle.reconstruction import Reconstructor
        from pyrticle.pusher import Pusher
        from pyrticle.cloud import ElementFinder
        from pyrticle.distribution import ParticleDistribution

        assert isinstance(setup.reconstructor, Reconstructor), \
                "must specify valid reconstructor"
        assert isinstance(setup.pusher, Pusher), \
                "must specify valid reconstructor"
        assert isinstance(setup.finder, ElementFinder), \
                "must specify valid element finder"
        assert isinstance(setup.distribution, ParticleDistribution), \
                "must specify valid particle distribution"
        assert isinstance(setup.element_order, int), \
                "must specify valid element order"
        assert isinstance(setup.dimensions_pos, int), \
                "must specify valid positional dimension count"
        assert isinstance(setup.dimensions_velocity, int), \
                "must specify valid positional dimension count"




class PICRunner(object):
    def __init__(self):
        from pyrticle.units import SI
        units = SI()
        self.units = SI()

        ui = PICCPyUserInterface(units)
        setup = self.setup = ui.gather()

        from pytools.log import LogManager
        import os.path
        self.logmgr = LogManager(os.path.join(
            setup.output_path, "pic.dat"), "w")

        from hedge.parallel import guess_parallelization_context
        self.pcon = guess_parallelization_context()

        if self.pcon.is_head_rank:
            mesh = self.pcon.distribute_mesh(setup.mesh)
        else:
            mesh = self.pcon.receive_mesh()

        self.discr = discr = \
                self.pcon.make_discretization(mesh, 
                        order=setup.element_order,
                        debug=setup.dg_debug)

        self.logmgr.set_constant("elements_total", len(setup.mesh.elements))
        self.logmgr.set_constant("elements_local", len(mesh.elements))
        self.logmgr.set_constant("element_order", setup.element_order)

        # em operator ---------------------------------------------------------
        maxwell_kwargs = {
                "epsilon": units.EPSILON0, 
                "mu": units.MU0, 
                "flux_type": setup.maxwell_flux_type,
                "bdry_flux_type": setup.maxwell_bdry_flux_type
                }

        if discr.dimensions == 3:
            from hedge.pde import MaxwellOperator
            self.max_op = MaxwellOperator(**maxwell_kwargs)
        elif discr.dimensions == 2:
            from hedge.pde import TEMaxwellOperator
            self.max_op = TEMaxwellOperator(**maxwell_kwargs)
        else:
            raise ValueError, "invalid mesh dimension"

        if setup.chi is not None:
            from pyrticle.hyperbolic import ECleaningMaxwellOperator
            self.max_op = ECleaningMaxwellOperator(self.max_op, 
                    chi=setup.chi, 
                    phi_decay=setup.phi_decay)

            if setup.phi_filter is not None:
                from pyrticle.hyperbolic import PhiFilter
                from hedge.discretization import Filter, ExponentialFilterResponseFunction
                em_filters.append(
                        PhiFilter(max_op, Filter(discr,
                            ExponentialFilterResponseFunction(*setup.phi_filter))))

        # timestepping setup --------------------------------------------------
        goal_dt = discr.dt_factor(self.max_op.max_eigenvalue())
        self.nsteps = int(setup.final_time/goal_dt)+1
        self.dt = setup.final_time/self.nsteps

        from hedge.timestep import RK4TimeStepper
        self.stepper = RK4TimeStepper()

        # particle setup ------------------------------------------------------
        from pyrticle.cloud import ParticleCloud, \
                optimize_shape_bandwidth, \
                guess_shape_bandwidth

        cloud = self.cloud = ParticleCloud(discr, units, 
                setup.reconstructor, setup.pusher, setup.finder,
                dimensions_pos=setup.dimensions_pos, 
                dimensions_velocity=setup.dimensions_velocity, 
                debug=setup.debug)

        cloud.add_particles( 
                setup.distribution.generate_particles(),
                setup.nparticles)

        self.total_charge = setup.nparticles*setup.distribution.mean()[2][0]
        if isinstance(setup.shape_bandwidth, str):
            if setup.shape_bandwidth == "optimize":
                optimize_shape_bandwidth(cloud, 
                        setup.distribution.get_rho_interpolant(
                            discr, self.total_charge),
                        setup.shape_exponent)
            elif setup.shape_bandwidth == "guess":
                guess_shape_bandwidth(cloud, setup.shape_exponent)
            else:
                raise ValueError, "invalid shape bandwidth setting '%s'" % (
                        setup.shape_bandwidth)
        else:
            from pyrticle._internal import PolynomialShapeFunction
            cloud.reconstructor.set_shape_function(
                    PolynomialShapeFunction(
                        float(setup.shape_bandwidth),
                        cloud.mesh_data.dimensions,
                        setup.shape_exponent,
                        ))

        # initial condition ---------------------------------------------------
        if "no_ic" in setup.debug:
            from pyrticle.cloud import FieldsAndCloud
            e, h = self.max_op.split_eh(self.max_op.assemble_fields(discr=discr))
            self.fields = FieldsAndCloud(self.max_op, e, h, cloud)
        else:
            from pyrticle.cloud import compute_initial_condition
            self.fields = compute_initial_condition(self.pcon, discr, cloud, 
                    max_op=self.max_op, 
                    potential_bc=setup.potential_bc, 
                    force_zero=False)

        # instrumentation setup -----------------------------------------------
        self.add_instrumentation(self.logmgr)

    def add_instrumentation(self, logmgr):
        from pytools.log import \
                add_simulation_quantities, \
                add_general_quantities, \
                add_run_info, ETA
        from pyrticle.log import add_particle_quantities, add_field_quantities, \
                add_beam_quantities, add_currents

        setup = self.setup

        add_run_info(logmgr)
        add_general_quantities(logmgr)
        add_simulation_quantities(logmgr, self.dt)
        add_particle_quantities(logmgr, self.cloud)
        add_field_quantities(logmgr, self.fields)
        if setup.beam_axis is not None and setup.beam_diag_axis is not None:
            add_beam_quantities(logmgr, self.cloud, 
                    axis=setup.beam_diag_axis, 
                    beam_axis=setup.beam_axis)
        if setup.tube_length is not None:
            from hedge.tools import unit_vector
            add_currents(logmgr, self.fields, 
                    unit_vector(self.cloud.dimensions_velocity, setup.beam_axis), 
                    setup.tube_length)

        self.stepper.add_instrumentation(logmgr)
        self.fields.add_instrumentation(logmgr)

        mean_beta = self.cloud.mean_beta()
        gamma = self.cloud.units.gamma_from_beta(mean_beta)

        logmgr.set_constant("dt", self.dt)
        logmgr.set_constant("beta", mean_beta)
        logmgr.set_constant("gamma", gamma)
        logmgr.set_constant("v", mean_beta*self.units.VACUUM_LIGHT_SPEED)
        logmgr.set_constant("Q0", self.total_charge)
        logmgr.set_constant("n_part_0", setup.nparticles)
        logmgr.set_constant("pmass", setup.distribution.mean()[3][0])
        logmgr.set_constant("chi", setup.chi)
        logmgr.set_constant("shape_radius_setup", setup.shape_bandwidth)
        logmgr.set_constant("shape_radius", self.cloud.reconstructor.shape_function.radius)
        logmgr.set_constant("shape_exponent", self.cloud.reconstructor.shape_function.exponent)

        from pytools.log import IntervalTimer
        self.vis_timer = IntervalTimer("t_vis", "Time spent visualizing")
        logmgr.add_quantity(self.vis_timer)

        logmgr.add_quantity(ETA(self.nsteps))

        logmgr.add_watches(setup.watch_vars)

    def run(self):
        t = 0
        
        setup = self.setup
        setup.hook_startup(self)

        vis_order = setup.vis_order
        if vis_order is None:
            vis_order = setup.element_order

        if vis_order != setup.element_order:
            vis_discr = self.pcon.make_discretization(self.discr.mesh, 
                            order=vis_order, debug=setup.dg_debug)

            from hedge.discretization import Projector
            vis_proj = Projector(self.discr, vis_discr)
        else:
            vis_discr = self.discr

            def vis_proj(f):
                return f

        from hedge.visualization import SiloVisualizer
        vis = SiloVisualizer(vis_discr)

        for step in xrange(self.nsteps):
            self.logmgr.tick()

            self.fields.upkeep()
            setup.hook_before_step(self)
            self.fields = self.stepper(self.fields, t, self.dt, self.fields.rhs)
            setup.hook_after_step(self)

            if step % setup.vis_interval == 0:
                self.vis_timer.start()
                import os.path
                visf = vis.make_file(os.path.join(
                    setup.output_path, setup.vis_pattern % step))

                self.cloud.add_to_vis(vis, visf, time=t, step=step)
                vis.add_data(visf, 
                        [(name, vis_proj(fld))
                            for name, fld in setup.hook_vis_quantities(self)],
                        time=t, step=step)
                setup.hook_visualize(self, vis, visf)

                visf.close()
                self.vis_timer.stop()

            t += self.dt

        vis.close()
            
        self.logmgr.tick()
        self.logmgr.save()

        setup.hook_when_done(self)



if __name__ == "__main__":
    PICRunner().run()

