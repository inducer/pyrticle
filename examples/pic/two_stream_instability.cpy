"""1D - two stream instability test case:
taken from: G.B. Jacobs, J.S. Hesthaven / JCP 214 (2006), High-order nodal 
discontious Galerkinn particle-in-cell method on unstructured grids

Comments: Even thought this is a 1D test case it is performed in a full two-
dimensional solver. This has been done in the referenced paper as well.
Quotation: [p.111]
"..., the computational domain has a length of 2*pi in the x-direction of the
plasma wave propagation. In the y-direction the grid has length 1.5 an is meshed
with approximately two triangles so as to simulate a 1D setting, i.e., a full
two-dimensional solver is used in this case."

Decription of the test case:
256 particles with R = 0.5 are released according to:

                      x = x_eq + A sin(k x_eq),

with A = 0.0001, k = 2 and a unit velocity.  Another 256 are released with A =
-0.0001 and unit velocity in the opposite direction.

Observed quantities:
- total energy
- kinetic energy
- potential energy

Expected result:
Appearance of the instability in the sudden drop of the kinetic energy.

"""

from numpy import pi as _pi
from numpy import sin as _sin
import numpy as _nu

#profile_output_filename = "wave.log"

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
vis_pattern =  "two_stream_inst-%04d"
vis_order =  None

# -----------------------------------------------------------------------------
# geometry and field discretization
# -----------------------------------------------------------------------------
element_order = 5
element_order = 4

shape_exponent = 2
shape_bandwidth = 1

#potential_bc = hedge.data.ConstantGivenFunction()

_tube_length = 2 * _pi
_tube_width = 1.5
import hedge.mesh as _mesh

mesh = _mesh.make_rect_mesh(
        [0,-_tube_width/2],
        [2*_pi,_tube_width/2],
        periodicity=(True, True),
        subdivisions=(30,10),
        max_area=0.2)

# -----------------------------------------------------------------------------
# timestepper setup
# -----------------------------------------------------------------------------
timestepper_order = 3

from hedge.timestep.multirate_ab.methods import methods as _methods
if False:
    _methods_man = ['f_f_1a', 'f_f_1b',
            's_f_1', 's_f_1_nr',
            's_f_2a', 's_f_2a_nr',
            's_f_2b', 's_f_2b_nr',
            's_f_3a', 's_f_3a_nr',
            's_f_3b', 's_f_3b_nr',
            's_f_4', 's_f_4_nr']

from hedge.timestep.multirate_ab import \
        TwoRateAdamsBashforthTimeStepper as _TwoRateAB

_enable_multirate = False
if _enable_multirate:
    def _dt_getter(discr, op, order):
       from hedge.timestep import AdamsBashforthTimeStepper
       # Take the stable step size of single-rate AB
       return discr.dt_factor(op.max_eigenvalue(),
               AdamsBashforthTimeStepper, 
               order)

    dt_getter = _dt_getter

    _step_ratio = 100
    dt_scale = 0.8
    timestepper_maker = lambda dt:_TwoRateAB(
            'f_f_1a',
            dt,
            _step_ratio,
            order=timestepper_order)

# -----------------------------------------------------------------------------
# particle setup
# -----------------------------------------------------------------------------

# get 25 particles in y-direction for a 1D equivalent setup in a 2D environment
_npart_y = 25
# Number of particles in x direction:
_npart_x = 256

nparticles = _npart_y * _npart_x *2
#_part_charge = [0.001177]
_part_charge = [1e-09]
_part_m = [5.6856296568531526e-21]
#_part_m = _part_charge

watch_vars = ["step", "t_sim", "W_field", "t_step", "t_eta", "n_part"]


# -----------------------------------------------------------------------------
# particle distribution setup
# -----------------------------------------------------------------------------
_x_equ_dist = _tube_length/_npart_x
_y_equ_dist = _tube_width/_npart_y

def _x_dist_f_1(_x):
    return _x + 0.0001 * _sin(2*_x)

def _x_dist_f_2(_x):
    return _x - 0.0001 * _sin(2*_x)

_p_state = []
_unit_v = 1

# Release 256 particles in one direction with unit speed.
for _y in range(_npart_y+1):
    for _x in range(_npart_x):
        _p_state.append(([_x_dist_f_1(_x*_x_equ_dist),_y*_y_equ_dist-_tube_width/2], \
                [_unit_v,0], \
                _part_charge, \
                _part_m))

# Release 256 particles in the opposite direction.
for _y in range(_npart_y+1):
    for _x in range(_npart_x):
        _p_state.append(([_x_dist_f_2(_x*_x_equ_dist),_y*_y_equ_dist-_tube_width/2], \
                [_unit_v*(-1),0], \
                _part_charge, \
                _part_m))


distribution = pyrticle.distribution.PlasmaWaveDistribution(_p_state, \
        _part_charge, \
        _part_m)

