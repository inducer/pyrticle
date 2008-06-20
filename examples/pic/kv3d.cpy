import random as _random
_random.seed(0)

pusher = PushMonomial()
reconstructor = RecGrid(
        FineCoreBrickGenerator(core_axis=2),
        el_tolerance=0.1,
        method="simplex_reduce")

dimensions_pos = 3
dimensions_velocity = 3

beam_axis = 2
beam_diag_axis = 0
tube_length = 100*units.MM

_cloud_charge = -10e-9 * units.C
nparticles = 20000
element_order = 3
final_time = 0.1*units.M/units.VACUUM_LIGHT_SPEED
_electrons_per_particle = abs(_cloud_charge/nparticles/units.EL_CHARGE)

_el_energy = units.EL_REST_ENERGY*10
_gamma = _el_energy/units.EL_REST_ENERGY
_mean_beta = (1-1/_gamma**2)**0.5

def hook_when_done(runner):
    _, _, err_table = runner.logmgr.get_expr_dataset("(rx_rms-rx_rms_theory)/rx_rms_theory")
    print "Relative error (rms): %g" % max(err for step, err in err_table)

def hook_startup(runner):
    from pyrticle.distribution import KVPredictedRadius
    runner.logmgr.add_quantity(KVPredictedRadius(runner.dt, 
        beam_v=distribution.beta*units.VACUUM_LIGHT_SPEED,
        predictor=distribution.get_rms_predictor(axis=0),
        suffix="x_rms"))
    runner.logmgr.add_quantity(KVPredictedRadius(runner.dt, 
        beam_v=distribution.beta*units.VACUUM_LIGHT_SPEED,
        predictor=distribution.get_total_predictor(axis=0),
        suffix="x_total"))

mesh = pyrticle.geometry.make_cylinder_with_fine_core(
        r=25*units.MM, inner_r=2.5*units.MM, 
        min_z=0, max_z=tube_length,
        max_volume_inner=10*units.MM**3,
        max_volume_outer=100*units.MM**3,
        radial_subdiv=4)

distribution = pyrticle.distribution.KVZIntervalBeam(
        units, total_charge=_cloud_charge, 
        p_charge=_cloud_charge/nparticles, 
        p_mass=_electrons_per_particle*units.EL_MASS,
        radii=[2.5*units.MM]*2,
        emittances=[5*units.MM*units.MRAD]*2, 
        z_length=5*units.MM,
        z_pos=10*units.MM,
        beta=_mean_beta)
