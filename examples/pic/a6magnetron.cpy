import random as _random
_random.seed(0)

pusher = PushMonomial()
reconstructor = RecGrid(
        #el_tolerance=0.1,
        method="simplex_reduce",
        jiggle_radius=0.0)
#reconstructor = RecAdv()
#reconstructor = RecShape()
#reconstructor = RecGridFind()

debug.add("shape_bw")

dimensions_pos = 2
dimensions_velocity = 2

beam_axis = 0
beam_diag_axis = 1
tube_length = 2

shape_bandwidth = "optimize"

_cloud_charge = 10e-9 * units.C
nparticles = 1
element_order = 3
final_time = 10*units.M/units.VACUUM_LIGHT_SPEED
_electrons_per_particle = abs(_cloud_charge/nparticles/units.EL_CHARGE)

def _make_a6():
    import sys, os
    sys.path.append(os.path.join(os.getcwd(), "mesh"))

    from magnetron import A6Triangulator as _A6Triangulator
    return _A6Triangulator()

_a6 = _make_a6()

def _make_mesh():
    tri_out = _a6.make_triangulation(1e-6)

    from hedge.mesh import MeshPyFaceMarkerLookup
    fmlookup = MeshPyFaceMarkerLookup(tri_out)

    marker2tags = {
            _a6.anode_marker: ["anode"],
            _a6.cathode_marker: ["cathode"],
            _a6.open_marker: ["open"],
            }

    def boundary_tagger(fvi, el, fn):
        return marker2tags[fmlookup(fvi)]

    from hedge.mesh import make_conformal_mesh

    return make_conformal_mesh(
            tri_out.points,
            tri_out.elements,
            boundary_tagger)

mesh = _make_mesh()

def _make_potential():
    from magnetron import A6Triangulator

    def pot(x, el):
        if la.norm(x) > (_a6.radius_cathode+_a6.radius_anode)/2:
            return 0
        else:
            return -100

    from hedge.data import GivenFunction
    return GivenFunction(pot)

potential_bc = _make_potential()

_c0 = units.VACUUM_LIGHT_SPEED

_mean_v = numpy.array([_c0*0.9,0])
_sigma_v = numpy.array([_c0*0.9*1e-3, _c0*1e-5])

_mean_beta = _mean_v/units.VACUUM_LIGHT_SPEED
_gamma = units.gamma_from_v(_mean_v)
_pmass = _electrons_per_particle*units.EL_MASS
_mean_p = _gamma*_pmass*_mean_v

distribution = pyrticle.distribution.JointParticleDistribution([
    pyrticle.distribution.GaussianPos([0.0,0.0], [0.01, 0.01]),
    pyrticle.distribution.GaussianMomentum(
        #_mean_p, _sigma_v*_gamma*_pmass, 
        0*_mean_p, 1e-10*_sigma_v*_gamma*_pmass, 
        units,
        pyrticle.distribution.DeltaChargeMass(
            _cloud_charge/nparticles,
            _pmass))
    ])

vis_interval = 10

if isinstance(reconstructor, RecGrid):
    def hook_visualize(runner, vis, visf):
        rec = runner.cloud.reconstructor
        rec.visualize_grid_quantities(visf, [
                ("rho_grid", rec.reconstruct_grid_rho()),
                ("j_grid", rec.reconstruct_grid_j(runner.cloud.velocities())),
                ("ones_resid", rec.remap_residual(rec.ones_on_grid())),
                ("rho_resid", rec.remap_residual(rec.reconstruct_grid_rho())),
                ("usecount", rec.grid_usecount()),
                ])

