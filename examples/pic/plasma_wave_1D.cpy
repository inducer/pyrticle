"""1D Plasma Wave test case:
taken from: G.B. Jacobs, J.S. Hesthaven / JCP 214 (2006), High-order nodal 
discontious Galerkinn particle-in-cell method on unstructured grids

Comments: Even thought this is a 1D test case it is performed in a full two-
dimensional solver. This has been done in the referenced paper as well. 
Quotation: [p.111]
"..., the computational domain has a length of 2*pi in the x-direction of the 
plasma wave propagation. In the y-direction the grid has length 1.5 an is meshed 
with approximately two triangles so as to simulate a 1D setting, i.e., a full
two-dimensional solver is used in this case."
"""

from numpy import pi as _pi
from numpy import sin as _sin
import numpy as _nu

# -----------------------------------------------------------------------------
# pic setup
# -----------------------------------------------------------------------------
pusher = PushMonomial()
depositor = DepShape()

dimensions_pos = 2
dimensions_velocity = 2

#final_time = 100 * units.M/units.VACUUM_LIGHT_SPEED
final_time = 0.01

vis_interval = 10
vis_pattern =  "plasma_wave_1D-%04d"
vis_order =  None

# -----------------------------------------------------------------------------
# geometry and field discretization
# -----------------------------------------------------------------------------
element_order = 5
element_order = 4

shape_exponent = 10
shape_bandwidth = 1

#potential_bc = hedge.data.ConstantGivenFunction()

_tube_length = 2 * _pi
_tube_width = 1.5
import hedge.mesh as _mesh

mesh = _mesh.make_rect_mesh(
        [0,-_tube_width/2],
        [2*_pi,_tube_width/2],
        periodicity=(True, True),
        subdivisions=(31,3),
        max_area=0.2)

#mesh = _mesh.make_regular_rect_mesh(
#        [-0.1,-0.7],
#        [_tube_length+0.1,8],
#            n=(50,3),
#            periodicity=(True, False))
            #max_area=0.15)       

# -----------------------------------------------------------------------------
# timestepper setup
# -----------------------------------------------------------------------------
from hedge.timestep import TwoRateAdamsBashforthTimeStepper as _TwoRateAB
_enable_multirate = False
if _enable_multirate:
    _step_ratio = 10
    timestepper_maker = lambda dt: _TwoRateAB(dt, step_ratio=_step_ratio, order=5,
            #largest_first=True)
            slowest_first=True)
    _rk4_stability = 2
    #dt_scale = 0.36/_rk4_stability*_step_ratio # ab4
    dt_scale = 0.18/_rk4_stability*_step_ratio # ab5


# -----------------------------------------------------------------------------
# particle setup
# -----------------------------------------------------------------------------

# get 25 particles in y-direction for a 1D equivalent setup in a 2D environment
_npart_y = 25
# Number of particles in x direction:
_npart_x = 320

nparticles = _npart_y * _npart_x 
_part_charge = [0.001177]
#_part_charge = [1e-09]
#_part_m = [5.6856296568531526e-21]
_part_m = _part_charge

watch_vars = ["step", "t_sim", "W_field", "t_step", "t_eta", "n_part"]


# -----------------------------------------------------------------------------
# particle distribution setup
# -----------------------------------------------------------------------------
_x_equ_dist = _tube_length/_npart_x
_y_equ_dist = _tube_width/_npart_y

def _x_dist_f(_x):
    return _x + 0.001 * _sin(2*_x)

_p_state = []

for _y in range(_npart_y+1):
    for _x in range(_npart_x):
        _p_state.append(([_x_dist_f(_x*_x_equ_dist),_y*_y_equ_dist-_tube_width/2], \
                [0.001,0], \
                _part_charge, \
                _part_m))


distribution = pyrticle.distribution.PlasmaWaveDistribution(_p_state, \
        _part_charge, \
        _part_m)

