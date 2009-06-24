# -----------------------------------------------------------------------------
# Two particles 2D test case:
# -----------------------------------------------------------------------------
# Idea: Two particles shall fly in a rectangular box with periodic BC's.
# -----------------------------------------------------------------------------
# pic setup
# -----------------------------------------------------------------------------
pusher = PushMonomial()
depositor = DepShape()

dimensions_pos = 2
dimensions_velocity = 2

final_time = 10 * units.M/units.VACUUM_LIGHT_SPEED

vis_interval = 10
# -----------------------------------------------------------------------------
# geometry and field discretization
# -----------------------------------------------------------------------------
element_order = 3

#shape_bandwidth = "optimize"
shape_bandwidth = 0.1

import hedge.mesh as _mesh
mesh = _mesh.make_rect_mesh(
        [-1,-1],
        [1, 1],
        periodicity=(True, True),
        subdivisions=(10,5),
        max_area=0.02)

# -----------------------------------------------------------------------------
# timestepper setup
# -----------------------------------------------------------------------------
from hedge.timestep import TwoRateAdamsBashforthTimeStepper as _TwoRateAB
_enable_multirate = False
if _enable_multirate:
    _step_ratio = 10
    timestepper_maker = lambda dt: _TwoRateAB(dt, step_ratio=_step_ratio, order=5,
            #largest_first=True)
            fastest_first=True)
            #slowest_first=True)
    _rk4_stability = 2
    #dt_scale = 0.36/_rk4_stability*_step_ratio # ab4
    dt_scale = 0.18/_rk4_stability*_step_ratio # ab5

# -----------------------------------------------------------------------------
# particle setup
# -----------------------------------------------------------------------------
nparticles = 2

# Create a list of spatial coordinates:
_x = [[-0.5,0],
        [0.5,0.0],
        [0.0,0.0]]
# Create a list of velocities:
_v = [[0,0],
          [0,0],
          [0,0]]
_q = [1e-09]
_m = [5.6856296568531526e-21]


distribution = pyrticle.distribution.ManualParticleDistribution(_x, _v, _q, _m)

if False:
    def hook_startup(runner):
        # Here to set the particles manually without any distribution:
        runner.method.add_particles(runner.state, 
            [
                (numpy.array([1,2,3]), numpy.array([4,5,6]), 1, 1),
                (numpy.array([1,2,3]), numpy.array([4,5,6]), 1, 1),
                (numpy.array([1,2,3]), numpy.array([4,5,6]), 1, 1),
                (numpy.array([1,2,3]), numpy.array([4,5,6]), 1, 1),
                ])

