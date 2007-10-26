from math import pi, sqrt




if True:
    RAD = 1
    M = 1
    N = 1
    A = 1
    S = 1
    KG = 1
    C = A*S
    J = N*M
    V = J/C
    F = C/V
else:
    from unum.units import *




MM = 1e-3 * M
MRAD = 1e-3 * RAD
M_S = M/S

EPSILON0 = 8.8541878176e-12 * F / M
MU0 = 4*pi*1e-7 * N / A**2
VACUUM_LIGHT_SPEED = 1/(EPSILON0*MU0)**0.5

EL_MASS = 9.10938215e-31 * KG
EL_CHARGE = 1.602176487e-19 * C
EL_REST_ENERGY = EL_MASS*VACUUM_LIGHT_SPEED**2

EV = EL_CHARGE * V
