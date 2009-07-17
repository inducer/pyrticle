from __future__ import division
import numpy
import numpy.linalg as la
import pytools




class PICCPyUserInterface(pytools.CPyUserInterface):
    def __init__(self, units):
        from pyrticle.deposition.shape import \
                ShapeFunctionDepositor, \
                NormalizedShapeFunctionDepositor
        from pyrticle.deposition.advective import \
                AdvectiveDepositor
        from pyrticle.deposition.grid import \
                GridDepositor
        from pyrticle.deposition.grid_find import \
                GridFindDepositor
        from pyrticle.deposition.grid_base import \
                SingleBrickGenerator, \
                FineCoreBrickGenerator

        from pyrticle.pusher import \
                MonomialParticlePusher, \
                AverageParticlePusher

        import pyrticle.geometry
        import pyrticle.distribution

        constants = {
                "numpy": numpy,
                "la": numpy.linalg,

                "units": units,
                "pyrticle": pyrticle,

                "DepShape": ShapeFunctionDepositor,
                "DepNormShape": NormalizedShapeFunctionDepositor,
                "DepAdv": AdvectiveDepositor,
                "DepGrid": GridDepositor,
                "DepGridFind": GridFindDepositor,

                "SingleBrickGenerator": SingleBrickGenerator,
                "FineCoreBrickGenerator": FineCoreBrickGenerator,

                "PushMonomial": MonomialParticlePusher,
                "PushAverage": AverageParticlePusher,
                }

        import hedge.data

        from hedge.timestep import RK4TimeStepper

        variables = {
                "pusher": None,
                "depositor": None,

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
                "profile_output_filename": None,

                "watch_vars": ["step", "t_sim", "W_field", "t_step", "t_eta", "n_part"],

                "hook_startup": lambda runner: None,
                "hook_after_step": lambda runner, state: None,
                "hook_when_done": lambda runner: None,
                "hook_vis_quantities": lambda observer: [
                    ("e", observer.e), 
                    ("h", observer.h), 
                    ("j", observer.method.deposit_j(observer.state)), 
                    ],
                "hook_visualize": lambda runner, vis, visf, observer: None,

                "timestepper_maker": lambda dt: RK4TimeStepper(),
                "dt_scale": 1,
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

        from pyrticle.deposition import Depositor
        from pyrticle.pusher import Pusher
        from pyrticle.distribution import ParticleDistribution

        assert isinstance(setup.depositor, Depositor), \
                "must specify valid depositor"
        assert isinstance(setup.pusher, Pusher), \
                "must specify valid depositor"
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
        from pyrticle.units import SIUnitsWithNaturalConstants
        self.units = units = SIUnitsWithNaturalConstants()

        ui = PICCPyUserInterface(units)
        setup = self.setup = ui.gather()

        from pytools.log import LogManager
        import os.path
        self.logmgr = LogManager(os.path.join(
            setup.output_path, "pic.dat"), "w")

        from hedge.backends import guess_run_context
        self.rcon = guess_run_context([])

        if self.rcon.is_head_rank:
            mesh = self.rcon.distribute_mesh(setup.mesh)
        else:
            mesh = self.rcon.receive_mesh()

        self.discr = discr = \
                self.rcon.make_discretization(mesh, 
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
            from hedge.models.em import MaxwellOperator
            self.maxwell_op = MaxwellOperator(**maxwell_kwargs)
        elif discr.dimensions == 2:
            from hedge.models.em import TEMaxwellOperator
            self.maxwell_op = TEMaxwellOperator(**maxwell_kwargs)
        else:
            raise ValueError, "invalid mesh dimension"

        if setup.chi is not None:
            from pyrticle.hyperbolic import ECleaningMaxwellOperator
            self.maxwell_op = ECleaningMaxwellOperator(self.maxwell_op, 
                    chi=setup.chi, 
                    phi_decay=setup.phi_decay)

            if setup.phi_filter is not None:
                from pyrticle.hyperbolic import PhiFilter
                from hedge.discretization import Filter, ExponentialFilterResponseFunction
                em_filters.append(
                        PhiFilter(maxwell_op, Filter(discr,
                            ExponentialFilterResponseFunction(*setup.phi_filter))))

        # timestepping setup --------------------------------------------------
        goal_dt = discr.dt_factor(self.maxwell_op.max_eigenvalue())
        self.nsteps = int(setup.final_time/goal_dt)+1
        self.dt = setup.final_time/self.nsteps * setup.dt_scale

        self.stepper = setup.timestepper_maker(self.dt)

        # particle setup ------------------------------------------------------
        from pyrticle.cloud import PicMethod, PicState, \
                optimize_shape_bandwidth, \
                guess_shape_bandwidth

        method = self.method = PicMethod(discr, units, 
                setup.depositor, setup.pusher,
                dimensions_pos=setup.dimensions_pos, 
                dimensions_velocity=setup.dimensions_velocity, 
                debug=setup.debug)

        self.state = method.make_state()
        method.add_particles( 
                self.state,
                setup.distribution.generate_particles(),
                setup.nparticles)

        self.total_charge = setup.nparticles*setup.distribution.mean()[2][0]
        if isinstance(setup.shape_bandwidth, str):
            if setup.shape_bandwidth == "optimize":
                optimize_shape_bandwidth(method, self.state,
                        setup.distribution.get_rho_interpolant(
                            discr, self.total_charge),
                        setup.shape_exponent)
            elif setup.shape_bandwidth == "guess":
                guess_shape_bandwidth(method, self.state, setup.shape_exponent)
            else:
                raise ValueError, "invalid shape bandwidth setting '%s'" % (
                        setup.shape_bandwidth)
        else:
            from pyrticle._internal import PolynomialShapeFunction
            method.depositor.set_shape_function(
                    self.state,
                    PolynomialShapeFunction(
                        float(setup.shape_bandwidth),
                        method.mesh_data.dimensions,
                        setup.shape_exponent,
                        ))

        # initial condition ---------------------------------------------------
        if "no_ic" in setup.debug:
            self.fields = self.maxwell_op.assemble_eh(discr=discr)
        else:
            from pyrticle.cloud import compute_initial_condition
            self.fields = compute_initial_condition(self.rcon, discr, method, self.state,
                    maxwell_op=self.maxwell_op, 
                    potential_bc=setup.potential_bc, 
                    force_zero=False)

        # rhs calculators -----------------------------------------------------
        from pyrticle.cloud import \
                FieldRhsCalculator, \
                FieldToParticleRhsCalculator, \
                ParticleRhsCalculator, \
                ParticleToFieldRhsCalculator
        self.f_rhs_calculator = FieldRhsCalculator(self.method, self.maxwell_op)
        self.p_rhs_calculator = ParticleRhsCalculator(self.method, self.maxwell_op)
        self.f2p_rhs_calculator = FieldToParticleRhsCalculator(self.method, self.maxwell_op)
        self.p2f_rhs_calculator = ParticleToFieldRhsCalculator(self.method, self.maxwell_op)

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

        from pyrticle.log import StateObserver
        self.observer = StateObserver(self.method, self.maxwell_op)
        self.observer.set_fields_and_state(self.fields, self.state)

        add_run_info(logmgr)
        add_general_quantities(logmgr)
        add_simulation_quantities(logmgr, self.dt)
        add_particle_quantities(logmgr, self.observer)
        add_field_quantities(logmgr, self.observer)

        if setup.beam_axis is not None and setup.beam_diag_axis is not None:
            add_beam_quantities(logmgr, self.observer, 
                    axis=setup.beam_diag_axis, 
                    beam_axis=setup.beam_axis)

        if setup.tube_length is not None:
            from hedge.tools import unit_vector
            add_currents(logmgr, self.observer, 
                    unit_vector(self.method.dimensions_velocity, setup.beam_axis), 
                    setup.tube_length)

        self.method.add_instrumentation(logmgr, self.observer)

        self.f_rhs_calculator.add_instrumentation(logmgr)

        if hasattr(self.stepper, "add_instrumentation"):
            self.stepper.add_instrumentation(logmgr)

        mean_beta = self.method.mean_beta(self.state)
        gamma = self.method.units.gamma_from_beta(mean_beta)

        logmgr.set_constant("dt", self.dt)
        logmgr.set_constant("beta", mean_beta)
        logmgr.set_constant("gamma", gamma)
        logmgr.set_constant("v", mean_beta*self.units.VACUUM_LIGHT_SPEED())
        logmgr.set_constant("Q0", self.total_charge)
        logmgr.set_constant("n_part_0", setup.nparticles)
        logmgr.set_constant("pmass", setup.distribution.mean()[3][0])
        logmgr.set_constant("chi", setup.chi)
        logmgr.set_constant("shape_radius_setup", setup.shape_bandwidth)
        logmgr.set_constant("shape_radius", self.method.depositor.shape_function.radius)
        logmgr.set_constant("shape_exponent", self.method.depositor.shape_function.exponent)

        from pytools.log import IntervalTimer
        self.vis_timer = IntervalTimer("t_vis", "Time spent visualizing")
        logmgr.add_quantity(self.vis_timer)

        logmgr.add_quantity(ETA(self.nsteps))

        logmgr.add_watches(setup.watch_vars)

    def inner_run(self): 
        t = 0
        
        setup = self.setup
        setup.hook_startup(self)

        vis_order = setup.vis_order
        if vis_order is None:
            vis_order = setup.element_order

        if vis_order != setup.element_order:
            vis_discr = self.rcon.make_discretization(self.discr.mesh, 
                            order=vis_order, debug=setup.dg_debug)

            from hedge.discretization import Projector
            vis_proj = Projector(self.discr, vis_discr)
        else:
            vis_discr = self.discr

            def vis_proj(f):
                return f

        from hedge.visualization import SiloVisualizer
        vis = SiloVisualizer(vis_discr)

        fields = self.fields
        self.observer.set_fields_and_state(fields, self.state)

        from hedge.tools import make_obj_array
        from pyrticle.cloud import TimesteppablePicState

        def visualize(observer):
            sub_timer = self.vis_timer.start_sub_timer()
            import os.path
            visf = vis.make_file(os.path.join(
                setup.output_path, setup.vis_pattern % step))

            self.method.add_to_vis(vis, visf, observer.state, time=t, step=step)
            vis.add_data(visf, 
                    [(name, vis_proj(fld))
                        for name, fld in setup.hook_vis_quantities(observer)],
                    time=t, step=step)
            setup.hook_visualize(self, vis, visf, observer)

            visf.close()
            sub_timer.stop().submit()

        from hedge.timestep.multirate_ab import TwoRateAdamsBashforthTimeStepper 
        if not isinstance(self.stepper, TwoRateAdamsBashforthTimeStepper): 
            def rhs(t, fields_and_state):
                fields, ts_state = fields_and_state
                state = ts_state.state

                fields_rhs = (
                        self.f_rhs_calculator(t, fields, state)
                        + self.p2f_rhs_calculator(t, fields, state))
                state_rhs = (
                        self.p_rhs_calculator(t, fields, state)
                        + self.f2p_rhs_calculator(t, fields, state))

                return make_obj_array([fields_rhs, state_rhs])
            step_args = (self.dt, rhs)
        else:
            def add_unwrap(rhs):
                def unwrapping_rhs(t, fields, ts_state):
                    if ts_state is not None:
                        state = ts_state.state
                    else:
                        state = None

                    return rhs(t, fields, state)
                return unwrapping_rhs

            step_args = ((
                    add_unwrap(self.f_rhs_calculator),
                    add_unwrap(self.p2f_rhs_calculator),
                    add_unwrap(self.f2p_rhs_calculator),
                    add_unwrap(self.p_rhs_calculator),
                    ),)

        y = make_obj_array([
            fields, 
            TimesteppablePicState(self.method, self.state)
            ])
        del self.state

        try:
            for step in xrange(self.nsteps):
                self.logmgr.tick()

                self.method.upkeep(y[1].state)

                y = self.stepper(y, t, *step_args)

                fields, ts_state = y
                self.observer.set_fields_and_state(fields, ts_state.state)

                setup.hook_after_step(self, self.observer)

                if step % setup.vis_interval == 0:
                    visualize(self.observer)

                t += self.dt
        finally:
            vis.close()
            self.discr.close()
            self.logmgr.save()

        setup.hook_when_done(self)

    def run(self):
        if self.setup.profile_output_filename is not None:
            from cProfile import Profile
            prof = Profile()
            try:
                prof.runcall(self.inner_run)
            finally:
                from lsprofcalltree import KCacheGrind
                kg = KCacheGrind(prof)
                kg.output(open(self.setup.profile_output_filename, "w"))
        else:
            self.inner_run()


if __name__ == "__main__":
    PICRunner().run()

