import random as _random
_random.seed(0)

pusher = PushMonomial()
reconstructor = RecGrid(
        el_tolerance=0.1,
        method="simplex_reduce")
#reconstructor = RecAdv()
#reconstructor = RecShape()
#reconstructor = RecGridFind()

dimensions_pos = 2
dimensions_velocity = 2

beam_axis = 0
beam_diag_axis = 1
tube_length = 2

shape_bandwidth = 0.2

_cloud_charge = -10e-9 * units.C
nparticles = 1
element_order = 3
final_time = 10*units.M/units.VACUUM_LIGHT_SPEED
_electrons_per_particle = abs(_cloud_charge/nparticles/units.EL_CHARGE)

_tube_width = 1
import hedge.mesh as _mesh
mesh = _mesh.make_rect_mesh(
        a=(-0.5, -_tube_width/2),
        b=(-0.5+tube_length, _tube_width/2),
        periodicity=(True, False),
        subdivisions=(10,5),
        max_area=0.02)

_c0 = units.VACUUM_LIGHT_SPEED

_mean_v = numpy.array([_c0*0.9,0])
_sigma_v = numpy.array([_c0*0.9*1e-3, _c0*1e-5])

_mean_beta = _mean_v/units.VACUUM_LIGHT_SPEED
_gamma = units.gamma_from_v(_mean_v)
_pmass = _electrons_per_particle*units.EL_MASS
_mean_p = _gamma*_pmass*_mean_v

distribution = pyrticle.distribution.JointParticleDistribution([
    pyrticle.distribution.GaussianPos([0,0], [0.1, 0.1]),
    pyrticle.distribution.GaussianMomentum(
        _mean_p, _sigma_v*_gamma*_pmass, units,
        pyrticle.distribution.DeltaChargeMass(
            _cloud_charge/nparticles,
            _pmass))
    ])

vis_interval = 10
