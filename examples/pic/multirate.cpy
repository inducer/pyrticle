import random as _random
_random.seed(0)

pusher = PushMonomial()
#reconstructor = RecGrid(
        ##el_tolerance=0.1,
        #method="simplex_reduce",
        #jiggle_radius=0.0)
#reconstructor = RecAdv()
reconstructor = RecShape()
#reconstructor = RecGridFind()

debug.add("shape_bw")
#debug.add("no_ic")

dimensions_pos = 2
dimensions_velocity = 2

beam_axis = 0
beam_diag_axis = 1
tube_length = 2

#shape_bandwidth = "optimize"
shape_bandwidth = 0.1

_cloud_charge = 10e-9 * units.C
nparticles = 1
element_order = 8
final_time = 100*units.M/units.VACUUM_LIGHT_SPEED
_electrons_per_particle = abs(_cloud_charge/nparticles/units.EL_CHARGE)

shape_exponent = 2

from hedge.timestep import TwoRateAdamsBashforthTimeStepper as _TwoRateAB
_step_ratio = 100
timestepper_maker = lambda dt: _TwoRateAB(dt, step_ratio=_step_ratio, order=4)
dt_scale = 0.36/2*_step_ratio

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
    pyrticle.distribution.GaussianPos([0.0,0.0], [0.01, 0.01]),
    pyrticle.distribution.GaussianMomentum(
        #_mean_p, _sigma_v*_gamma*_pmass, 
        _mean_p, 1e-10*_sigma_v*_gamma*_pmass, 
        units,
        pyrticle.distribution.DeltaChargeMass(
            _cloud_charge/nparticles,
            _pmass))
    ])

vis_interval = 10
vis_order = 8

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
