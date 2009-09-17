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
from numpy import sqrt as _sqrt
import numpy as _nu

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
final_time = 10000 * 1.60348185077e-13

vis_interval = 10
vis_pattern =  "plasma_wave_1D-%04d"
vis_order =  None
log_file = "plasma_wave_1D_TRAB4_10.dat"

def hook_vis_quantities(observer):
    try:
        observer.div_operator
    except AttributeError:
        from hedge.models.nd_calculus import DivergenceOperator

        observer.div_operator = DivergenceOperator(dimensions_pos).bind(observer.discr)

    return [
                    ("e", observer.e), 
                    ("h", observer.h), 
                    ("j", observer.method.deposit_j(observer.state)), 
                    ("rho", observer.method.deposit_rho(observer.state)), 
                    ("div_d", units.EPSILON0*observer.div_operator(observer.e)), 
                    ]

# -----------------------------------------------------------------------------
# geometry and field discretization
# -----------------------------------------------------------------------------
_mesh_scale = 0.01
element_order = 5

_max_area = 0.02 * _mesh_scale**2

shape_exponent = 2
shape_bandwidth = 2*_sqrt(_max_area)

ic_tol = 1e-7

#potential_bc = hedge.data.ConstantGivenFunction()

_tube_length = 2 * _pi * _mesh_scale
_tube_width = 1.5 * _mesh_scale
import hedge.mesh as _mesh

mesh = _mesh.make_rect_mesh(
        [0,-_tube_width/2],
        [_tube_length,_tube_width/2],
        periodicity=(True, True),
        subdivisions=(30,10),
        max_area=_max_area)


# -----------------------------------------------------------------------------
# timestepper setup
# -----------------------------------------------------------------------------

from hedge.timestep.multirate_ab.methods import methods as _methods
from hedge.timestep.multirate_ab import \
        TwoRateAdamsBashforthTimeStepper as _TwoRateAB

_enable_multirate = True

if _enable_multirate:
    timestepper_order = 4
    
    def _dt_getter(discr, op, order):
       from hedge.timestep import AdamsBashforthTimeStepper
       # Take the stable step size of single-rate AB
       return discr.dt_factor(op.max_eigenvalue(),
               AdamsBashforthTimeStepper, 
               order)

    dt_getter = _dt_getter
    
    _step_ratio = 10

    dt_scale = 0.6 * _step_ratio

    timestepper_maker = lambda dt:_TwoRateAB(
            'slowest_first_1',
            dt,
            _step_ratio,
            order=timestepper_order)

# -----------------------------------------------------------------------------
# particle setup
# -----------------------------------------------------------------------------

# get 25 particles in y-direction for a 1D equivalent setup in a 2D environment
_npart_y = 25
# Number of particles in x direction:
_npart_x = 200

nparticles = _npart_y * _npart_x

# Get the particle charge in order to obtain a period T of _steps_per_wave * _dt
_steps_per_wave = 1000
_c_charge_mass  = units.EL_CHARGE/units.EL_MASS
_number_density = nparticles / (_tube_width * _tube_length)
_dt = 1.60348185077e-13
_part_charge = (2*_pi/(_steps_per_wave * _dt))**2 * units.EPSILON0 \
        / (_number_density * _c_charge_mass)

#_part_charge = 0.000001
#_part_charge = 1.457e-6
_part_m = _part_charge * units.EL_MASS/units.EL_CHARGE

_omega = _sqrt(_number_density * _part_charge**2/(_part_m * units.EPSILON0 ))

_frequency = _omega/(2*_pi)
        
print "omega", _omega
print "f", _frequency
print "T", 1/_frequency
print "q_part", _part_charge

def rho_static_getter(discr):
    return -_nu.array([_number_density*_part_charge])


# -----------------------------------------------------------------------------
# particle distribution setup
# -----------------------------------------------------------------------------
_x_equ_dist = _tube_length/_npart_x
_y_equ_dist = _tube_width/_npart_y

def _x_dist_f(_x):
    return _x + 0.001 * _sin(2*_x/_mesh_scale) * _mesh_scale

_p_state = []

for _y in range(_npart_y+1):
    for _x in range(_npart_x):
        _p_state.append(([_x_dist_f(_x*_x_equ_dist),_y*_y_equ_dist-_tube_width/2], \
                [0,0], \
                [_part_charge], \
                [_part_m]))


distribution = pyrticle.distribution.PlasmaWaveDistribution(_p_state, \
        [_part_charge], \
        [_part_m])

def hook_startup(runner):
    runner.logmgr.set_constant("haumichblau", 5)
