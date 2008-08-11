# -----------------------------------------------------------------------------
# pic setup
# -----------------------------------------------------------------------------
pusher = PushMonomial()
reconstructor = RecShape()

dimensions_pos = 2
dimensions_velocity = 2

final_time = 10*units.M/units.VACUUM_LIGHT_SPEED

vis_interval = 10

# -----------------------------------------------------------------------------
# geometry and field discretization
# -----------------------------------------------------------------------------
element_order = 3

_tube_width = 1*units.M
_tube_length = 2*units.M

import hedge.mesh as _mesh
mesh = _mesh.make_rect_mesh(
        [0, -_tube_width/2],
        [_tube_length, _tube_width/2],
        periodicity=(True, False),
        subdivisions=(10,5),
        max_area=0.02)

# -----------------------------------------------------------------------------
# particle setup
# -----------------------------------------------------------------------------
_cloud_charge = 10e-9 * units.C
nparticles = 2000

_electrons_per_particle = abs(_cloud_charge/nparticles/units.EL_CHARGE)
_pmass = _electrons_per_particle*units.EL_MASS

_c0 = units.VACUUM_LIGHT_SPEED

_mean_v = numpy.array([_c0*0.9,0])
_sigma_v = numpy.array([_c0*0.9*1e-3, _c0*1e-5])

_mean_beta = _mean_v/units.VACUUM_LIGHT_SPEED
_gamma = units.gamma_from_v(_mean_v)
_mean_p = _gamma*_pmass*_mean_v

distribution = pyrticle.distribution.JointParticleDistribution([
    pyrticle.distribution.GaussianPos(
        [_tube_length*0.25,0.0], 
        [0.1, 0.1]),
    pyrticle.distribution.GaussianMomentum(
        _mean_p, _sigma_v*_gamma*_pmass, 
        units,
        pyrticle.distribution.DeltaChargeMass(
            _cloud_charge/nparticles,
            _pmass))
    ])
