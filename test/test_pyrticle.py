from __future__ import division




import numpy
import numpy.linalg as la
from pytools.test import mark_test




def test_variance():
    def sg_variance(iterable, entire_pop):
        from pyrticle._internal import StatsGatherer
        sg = StatsGatherer()
        for i in iterable:
            sg.add(i)

        if entire_pop:
            return sg.variance()
        else:
            return sg.variance_sample()


    for entire_pop in [False, True]:
        data = [4, 7, 13, 16]
        from pytools import variance
        orig_variance = variance(data, entire_pop)
        assert abs(sg_variance(data, entire_pop)
                - orig_variance) < 1e-15

        data = [1e9 + x for x in data]
        assert abs(sg_variance(data, entire_pop)
                - orig_variance) < 1e-15




def test_ode_defined_function():
    from pyrticle.tools import ODEDefinedFunction

    class Sin(ODEDefinedFunction):
        def __init__(self):
            ODEDefinedFunction.__init__(self, 0,
                    numpy.array([0,1], dtype=numpy.float64), 1/7*1e-2)

        def rhs(self, t, y):
            return numpy.array([y[1], -y[0]], dtype=numpy.float64)

        def __call__(self, t):
            return ODEDefinedFunction.__call__(self, t)[0]

    from math import pi
    s = Sin()
    assert abs(s(-pi)) < 2e-3
    assert abs(s(pi)) < 2e-3
    assert abs(s(-pi/2)+1) < 2e-3



def test_kv_predictors():
    from pyrticle.distribution import \
            ChargelessKVRadiusPredictor, KVRadiusPredictor
    kv_env_exact = ChargelessKVRadiusPredictor(2.5e-3, 5e-6)
    kv_env_num = KVRadiusPredictor(2.5e-3, 5e-6)

    from hedge.tools import plot_1d
    steps = 50
    for i in range(steps):
        s = kv_env_num.dt/7*i

        a_exact = kv_env_exact(s)
        a_num = kv_env_num(s)
        assert abs(a_exact-a_num)/a_exact < 1e-3




@mark_test.long
def test_kv_with_no_charge():
    from random import seed
    seed(0)

    from pyrticle.units import SIUnitsWithNaturalConstants
    units = SIUnitsWithNaturalConstants()

    # discretization setup ----------------------------------------------------
    from hedge.mesh import make_cylinder_mesh
    from hedge.backends import guess_run_context

    rcon = guess_run_context([])

    tube_length = 100*units.MM
    mesh = make_cylinder_mesh(radius=25*units.MM, height=tube_length, periodic=True)

    discr = rcon.make_discretization(mesh, order=3)

    dt = discr.dt_factor(units.VACUUM_LIGHT_SPEED()) / 2
    final_time = 1*units.M/units.VACUUM_LIGHT_SPEED()
    nsteps = int(final_time/dt)+1
    dt = final_time/nsteps

    # particles setup ---------------------------------------------------------
    from pyrticle.cloud import PicMethod
    from pyrticle.deposition.shape import ShapeFunctionDepositor
    from pyrticle.pusher import MonomialParticlePusher

    method = PicMethod(discr, units,
            ShapeFunctionDepositor(),
            MonomialParticlePusher(),
            3, 3)

    nparticles = 10000
    cloud_charge = 1e-9 * units.C
    electrons_per_particle = cloud_charge/nparticles/units.EL_CHARGE

    el_energy = 5.2e6 * units.EV
    gamma = el_energy/units.EL_REST_ENERGY()
    beta = (1-1/gamma**2)**0.5

    from pyrticle.distribution import KVZIntervalBeam
    beam = KVZIntervalBeam(units,
            total_charge=0,
            p_charge=0,
            p_mass=electrons_per_particle*units.EL_MASS,
            radii=2*[2.5*units.MM],
            emittances=2*[5 * units.MM * units.MRAD],
            z_length=5*units.MM,
            z_pos=10*units.MM,
            beta=beta)

    state = method.make_state()
    method.add_particles(state, beam.generate_particles(), nparticles)

    # diagnostics setup -------------------------------------------------------
    from pytools.log import LogManager
    from pyrticle.log import add_beam_quantities, StateObserver
    observer = StateObserver(method, None)
    logmgr = LogManager(mode="w")
    add_beam_quantities(logmgr, observer, axis=0, beam_axis=2)

    from pyrticle.distribution import KVPredictedRadius
    logmgr.add_quantity(KVPredictedRadius(dt,
        beam_v=beta*units.VACUUM_LIGHT_SPEED(),
        predictor=beam.get_rms_predictor(axis=0),
        suffix="x_rms"))
    logmgr.add_quantity(KVPredictedRadius(dt,
        beam_v=beta*units.VACUUM_LIGHT_SPEED(),
        predictor=beam.get_total_predictor(axis=0),
        suffix="x_total"))

    # timestep loop -----------------------------------------------------------
    vel = method.velocities(state)
    from hedge.tools import join_fields
    def rhs(t, y):
        return join_fields([
            vel,
            0*vel,
            0, # drecon
            ])

    from hedge.timestep import RK4TimeStepper
    stepper = RK4TimeStepper()
    t = 0

    from pyrticle.cloud import TimesteppablePicState
    ts_state = TimesteppablePicState(method, state)

    for step in xrange(nsteps):
        observer.set_fields_and_state(None, ts_state.state)

        logmgr.tick()

        ts_state = stepper(ts_state, t, dt, rhs)
        method.upkeep(ts_state.state)

        t += dt

    logmgr.tick()

    _, _, err_table = logmgr.get_expr_dataset("(rx_rms-rx_rms_theory)/rx_rms_theory")
    rel_max_rms_error = max(err for step, err in err_table)
    assert rel_max_rms_error < 0.01




@mark_test.long
def test_efield_vs_gauss_law():
    from hedge.mesh import \
            make_box_mesh, \
            make_cylinder_mesh
    from math import sqrt, pi
    from pytools.arithmetic_container import \
            ArithmeticList, join_fields
    from random import seed
    from pytools.stopwatch import Job

    from pyrticle.units import SIUnitsWithNaturalConstants
    units = SIUnitsWithNaturalConstants()

    seed(0)

    nparticles = 10000
    beam_radius = 2.5 * units.MM
    emittance = 5 * units.MM * units.MRAD
    final_time = 0.1*units.M/units.VACUUM_LIGHT_SPEED()
    field_dump_interval = 1
    tube_length = 20*units.MM

    # discretization setup ----------------------------------------------------
    from pyrticle.geometry import make_cylinder_with_fine_core
    mesh = make_cylinder_with_fine_core(
            r=10*beam_radius, inner_r=1*beam_radius,
            min_z=0, max_z=tube_length,
            max_volume_inner=10*units.MM**3,
            max_volume_outer=100*units.MM**3,
            radial_subdiv=10,
            )

    from hedge.backends import guess_run_context
    rcon = guess_run_context([])
    discr = rcon.make_discretization(mesh, order=3)

    from hedge.models.em import MaxwellOperator
    max_op = MaxwellOperator(
            epsilon=units.EPSILON0,
            mu=units.MU0,
            flux_type=1)

    from hedge.models.nd_calculus import DivergenceOperator
    div_op = DivergenceOperator(discr.dimensions)

    # particles setup ---------------------------------------------------------
    from pyrticle.cloud import PicMethod
    from pyrticle.deposition.shape import ShapeFunctionDepositor
    from pyrticle.pusher import MonomialParticlePusher

    method = PicMethod(discr, units,
            ShapeFunctionDepositor(),
            MonomialParticlePusher(),
            3, 3)

    # particle ic ---------------------------------------------------------
    cloud_charge = -1e-9 * units.C
    electrons_per_particle = abs(cloud_charge/nparticles/units.EL_CHARGE)

    el_energy = 10*units.EL_REST_ENERGY()
    el_lorentz_gamma = el_energy/units.EL_REST_ENERGY()
    beta = (1-1/el_lorentz_gamma**2)**0.5
    gamma = 1/sqrt(1-beta**2)

    from pyrticle.distribution import KVZIntervalBeam
    beam = KVZIntervalBeam(units, total_charge=cloud_charge,
            p_charge=cloud_charge/nparticles,
            p_mass=electrons_per_particle*units.EL_MASS,
            radii=2*[beam_radius],
            emittances=2*[5 * units.MM * units.MRAD],
            z_length=tube_length,
            z_pos=tube_length/2,
            beta=beta)

    state = method.make_state()

    method.add_particles(state, beam.generate_particles(), nparticles)

    # field ic ----------------------------------------------------------------
    from pyrticle.cloud import guess_shape_bandwidth
    guess_shape_bandwidth(method, state, 2)

    from pyrticle.cloud import compute_initial_condition

    from hedge.data import ConstantGivenFunction
    fields = compute_initial_condition(
            rcon,
            discr, method, state, maxwell_op=max_op,
            potential_bc=ConstantGivenFunction())

    # check against theory ----------------------------------------------------
    q_per_unit_z = cloud_charge/beam.z_length
    class TheoreticalEField:
        shape = (3,)

        def __call__(self, x, el):
            r = la.norm(x[:2])
            if r >= max(beam.radii):
                xy_unit = x/r
                xy_unit[2] = 0
                return xy_unit*((q_per_unit_z)
                        /
                        (2*pi*r*max_op.epsilon))
            else:
                return numpy.zeros((3,))

    def theory_indicator(x, el):
        r = la.norm(x[:2])
        if r >= max(beam.radii):
            return 1
        else:
            return 0

    from hedge.tools import join_fields, to_obj_array
    e_theory = to_obj_array(discr.interpolate_volume_function(TheoreticalEField()))
    theory_ind = discr.interpolate_volume_function(theory_indicator)

    e_field, h_field = max_op.split_eh(fields)
    restricted_e = join_fields(*[e_i * theory_ind for e_i in e_field])

    def l2_error(field, true):
        return discr.norm(field-true)/discr.norm(true)

    outer_l2 = l2_error(restricted_e, e_theory)
    assert outer_l2 < 0.08

    if False:
        visf = vis.make_file("e_comparison")
        mesh_scalars, mesh_vectors = \
                method.add_to_vis(vis, visf)
        vis.add_data(visf, [
            ("e", restricted_e),
            ("e_theory", e_theory),
            ]
            + mesh_vectors
            + mesh_scalars
            )
        visf.close()




@mark_test.long
def test_with_static_fields():
    from pyrticle.units import SIUnitsWithNaturalConstants

    units = SIUnitsWithNaturalConstants()

    from hedge.discretization.local import TetrahedronDiscretization
    from hedge.mesh import \
            make_box_mesh, \
            make_cylinder_mesh
    from hedge.discretization import Discretization

    # discretization setup ----------------------------------------------------
    radius = 1*units.M
    full_mesh = make_cylinder_mesh(radius=radius, height=2*radius, periodic=True,
            radial_subdivisions=30)

    from hedge.backends import guess_run_context

    pcon = guess_run_context([])

    if pcon.is_head_rank:
        mesh = pcon.distribute_mesh(full_mesh)
    else:
        mesh = pcon.receive_mesh()

    discr = pcon.make_discretization(mesh, order=1)

    # particles setup ---------------------------------------------------------
    def get_setup(case):
        c = units.VACUUM_LIGHT_SPEED()
        from static_field import LarmorScrew, EBParallel
        if case == "screw":
            return LarmorScrew(units,
                    mass=units.EL_MASS, charge=units.EL_CHARGE, c=c,
                    vpar=c*0.8, vperp=c*0.1, bz=1e-3,
                    nparticles=4)
        elif case == "epb":
            return EBParallel(units,
                    mass=units.EL_MASS, charge=units.EL_CHARGE, c=c,
                    ez=1e+5, bz=1e-3, radius=0.5*radius, nparticles=1)
        else:
            raise ValueError, "invalid test case"

    from pyrticle.pusher import \
            MonomialParticlePusher, \
            AverageParticlePusher

    from static_field import run_setup

    for pusher in [MonomialParticlePusher, AverageParticlePusher]:
        for case in ["screw", "epb"]:
            casename = "%s-%s" % (case, pusher.__class__.__name__.lower())
            run_setup(units, casename, get_setup(case), discr, pusher)




def test_shape_functions():
    from pyrticle.tools import \
            CInfinityShapeFunction, \
            PolynomialShapeFunction

    from hedge.mesh import \
            make_uniform_1d_mesh, \
            make_rect_mesh, make_box_mesh

    from hedge.backends import guess_run_context
    rcon = guess_run_context([])

    for r in [0.1, 10]:
        for mesh in [
                make_uniform_1d_mesh(-r, r, 10),
                make_rect_mesh(
                    (-r,-r), (r,r),
                    max_area=(r/10)**2),
                make_box_mesh(
                    (-r,-r,-r), (r,r,r),
                    max_volume=(r/10)**3),
                ]:
            discr = rcon.make_discretization(mesh, order=3)
            for sfunc in [
                    PolynomialShapeFunction(r, discr.dimensions, 2),
                    PolynomialShapeFunction(r, discr.dimensions, 4),
                    CInfinityShapeFunction(r, discr.dimensions),
                    ]:
                num_sfunc = discr.interpolate_volume_function(
                        lambda x, el: sfunc(x))
                int_sfunc = discr.integral(num_sfunc)
                assert abs(int_sfunc-1) < 4e-5




if __name__ == "__main__":
    from py.test.cmdline import main
    main([__file__])
