# -*- coding: utf8 -*-

"""2D - Weibel instability
taken from: R.L. Morse, C.W. Nielson, Numerical simulation of the Weibel
instability in one and two dimensions, Phys. Fluids 14 (4) (1971) 830– 840.

Decription of the test case: The Weibel instability simulations are performed
on a unit square with periodic boundary conditions. We consider a quasi-neutral
plasma with a thermal velocity ratio of 5 of the velocity in x, u_the = 0.25
and y, v_the = 0.05 direction. The plasma frequency is 15 times the length of
the square, i.e., ω_pe = 15 resulting in q/m = -(15 π) with the electron charge
density set to ρ = -1.  With these settings, magnetic waves develop with a
dominant frequency in the y-direction. The wave number decreases in time as
the thermal velocities approach the equilibrium state.  A study with the FDTD
method indicates that a 256 · 256 grid with 36 particles per cell yields a
reasonably converged solution. The results of this simulation are used in the
remainder of this section for comparison.  The high-order simulations are
performed on a grid with 200 triangles using a ﬁfth order scheme. We track
Np · Np particles in this domain for two time units.


Observed quantities:
- potential energy

Expected result:
With these settings, magnetic waves develop with a dominant frequency in the
y-direction. The wave number decreases in time as the thermal velocities
approach the equilibrium state.
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
vis_pattern =  "weibel-%04d"
vis_order =  None

# -----------------------------------------------------------------------------
# geometry and field discretization
# -----------------------------------------------------------------------------
element_order = 5

shape_exponent = 2
shape_bandwidth = 1

_x_len = 1
_y_len = 1

_area_total_domain = _x_len * _x_len
_number_of_cells = 200
_max_area = _area_total_domain/_number_of_cells

import hedge.mesh as _mesh
mesh = _mesh.make_rect_mesh(
        [0,0],
        [_x_len,_y_len],
        periodicity=(True, True),
        subdivisions=(30,30),
        max_area=_max_area)

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
# plasma frequency:
_omega_plasma = -(15*_pi)
# electron charge density:
_p_chrg = 1e-09
_part_charge = [1e-09]
# N_p:
_N_p = 300
# Number of particles:
_npart = _N_p * _N_p
nparticles = _npart
# particle mass:
_part_m = [_p_chrg/_omega_plasma]
# thermal velocities:
_v_x_thermal = 0.25
_v_y_thermal = 0.05

# -----------------------------------------------------------------------------
# particle distribution setup
# -----------------------------------------------------------------------------
_x_equ_dist = _x_len/_N_p
_y_equ_dist = _y_len/_N_p

_p_state = []

_a = 0
for _i in range(_N_p+1):
    _a += 1
    for _x in range(_N_p):
        _a += 1
        _y = _i * _y_equ_dist
        _p_state.append(([_x*_x_equ_dist,_y],
                [_v_x_thermal,_v_y_thermal],
                _part_charge,
                _part_m))
#print _a, _npart
#raw_input()

distribution = pyrticle.distribution.PlasmaWaveDistribution(_p_state, \
        _part_charge, \
        _part_m)


# -----------------------------------------------------------------------------
# Log setup:
# -----------------------------------------------------------------------------
watch_vars = ["step", "t_sim", "W_field", "t_step", "t_eta", "n_part"]


