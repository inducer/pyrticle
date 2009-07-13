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

Decription of the test case:
320 particles are equidistantly distributed with a superimposed one-dimensional
sine deviation:

                      x = x_eq + A sin(k x_eq),

where the amplitude of the deviation is A = 0.001 and the wavenumber k = 2. The
cloud influence area is R = 0.5. Physical parameters for the particle are q =
0.001177 and q/m = 1.0.

Observed quantities:
- total energy
- kinetic energy
- potential energy

Expected result:
Preserve energy equally well.
"""

from numpy import pi as _pi
from numpy import sin as _sin
import numpy as _nu

#profile_output_filename = "wave.log"

debug = ["vis_files", "ic", "poisson"]

# -----------------------------------------------------------------------------
# pic setup
# -----------------------------------------------------------------------------
pusher = PushMonomial()
depositor = DepGridFind()

dimensions_pos = 2
dimensions_velocity = 2

#final_time = 0.1*units.M/units.VACUUM_LIGHT_SPEED
#final_time = 100 * units.M/units.VACUUM_LIGHT_SPEED
final_time = 0.001

vis_interval = 10
vis_pattern =  "plasma_wave_1D-%04d"
vis_order =  None

def hook_vis_quantities(observer):
    return [
                    ("e", observer.e), 
                    ("h", observer.h), 
                    ("j", observer.method.deposit_j(observer.state)), 
                    ("rho", observer.method.deposit_rho(observer.state)), 
                    ]

# -----------------------------------------------------------------------------
# geometry and field discretization
# -----------------------------------------------------------------------------
element_order = 5

shape_exponent = 2
shape_bandwidth = 0.5

#potential_bc = hedge.data.ConstantGivenFunction()

_tube_length = 2 * _pi
_tube_width = 1.5
import hedge.mesh as _mesh

mesh = _mesh.make_rect_mesh(
        [0,-_tube_width/2],
        [2*_pi,_tube_width/2],
        periodicity=(True, True),
        subdivisions=(30,10),
        max_area=0.05)


# -----------------------------------------------------------------------------
# timestepper setup
# -----------------------------------------------------------------------------
timestepper_order = 5

def _dt_getter(discr, op, order):
    from hedge.timestep import AdamsBashforthTimeStepper
    return discr.dt_factor(op.max_eigenvalue(),
            AdamsBashforthTimeStepper, 
            order)

dt_getter = _dt_getter

from hedge.timestep import TwoRateAdamsBashforthTimeStepper as _TwoRateAB
_enable_multirate = True
if _enable_multirate:
    _step_ratio = 10 * 0.8
    dt_scale = _step_ratio
    timestepper_maker = lambda dt: _TwoRateAB(dt,
            step_ratio=_step_ratio,
            order=timestepper_order,
            slowest_first=True,
            fastest_first=False,
            substepping=False)


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

