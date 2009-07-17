import random as _random
_random.seed(0)

pusher = PushMonomial()
reconstructor = RecShape()

dimensions_pos = 3
dimensions_velocity = 3

beam_axis = 2
beam_diag_axis = 0
tube_length = 100*units.MM

_cloud_charge = -10e-9 * units.C
nparticles = 20000
element_order = 3
final_time = 0.1*units.M/units.VACUUM_LIGHT_SPEED()
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
        beam_v=distribution.beta*units.VACUUM_LIGHT_SPEED(),
        predictor=distribution.get_rms_predictor(axis=0),
        suffix="x_rms"))
    runner.logmgr.add_quantity(KVPredictedRadius(runner.dt, 
        beam_v=distribution.beta*units.VACUUM_LIGHT_SPEED(),
        predictor=distribution.get_total_predictor(axis=0),
        suffix="x_total"))

# -----------------------------------------------------------------------------
# mesh ------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def _get_aps_path():
    root_path = "."

    from os.path import join, isdir
    while not isdir(join(root_path, "examples")):
        root_path = join(root_path, "..")
    return join(root_path, "examples", "aps")



def _make_mesh():
    aps_path = _get_aps_path()

    from os.path import join, isdir

    import sys
    sys.path.append(aps_path)

    from superfish import parse_superfish_format
    rz = parse_superfish_format(join(aps_path, "gun.am"), 
            max_point_dist=0.3)
    from rzmesh import make_mesh_info_with_inner_tube
    mesh_info = make_mesh_info_with_inner_tube(rz, 
            tube_r=0.1, radial_subdiv=10,
            max_inner_volume=4e-4
            )
    from meshpy.tet import build
    generated_mesh = build(mesh_info, verbose=True, volume_constraints=True)
    print len(generated_mesh.elements)

    # generated_mesh is in cm, we use m
    points = numpy.array(generated_mesh.points)*1e-2

    from hedge.mesh import make_conformal_mesh
    return make_conformal_mesh(
            points,
            generated_mesh.elements,
            lambda fvi, el, fn: ["pec"],
            periodicity=None)

mesh = _make_mesh()

# -----------------------------------------------------------------------------
# particles -------------------------------------------------------------------
# -----------------------------------------------------------------------------
distribution = pyrticle.distribution.KVZIntervalBeam(
        units, total_charge=_cloud_charge, 
        p_charge=_cloud_charge/nparticles, 
        p_mass=_electrons_per_particle*units.EL_MASS,
        radii=[2.5*units.MM]*2,
        emittances=[5*units.MM*units.MRAD]*2, 
        z_length=5*units.MM,
        z_pos=10*units.MM,
        beta=_mean_beta)
