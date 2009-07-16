"""Systems of physical units of measurement"""

__copyright__ = "Copyright (C) 2007, 2008 Andreas Kloeckner"

__license__ = """
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see U{http://www.gnu.org/licenses/}.
"""




from math import pi, sqrt
import numpy
import numpy.linalg as la




class SIUnits:
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
    T = V*S/M**2

    MM = 1e-3 * M
    MRAD = 1e-3 * RAD
    M_S = M/S

    @classmethod
    def VACUUM_LIGHT_SPEED(cls):
        return 1/(cls.EPSILON0*cls.MU0)**0.5

    @classmethod
    def EL_REST_ENERGY(cls):
        return cls.EL_MASS*cls.VACUUM_LIGHT_SPEED()**2

    @classmethod
    def EV(cls):
        return cls.EL_CHARGE * V

    @classmethod
    def gamma_from_beta(cls, beta):
        return (1-numpy.dot(beta, beta))**(-0.5)

    @classmethod
    def gamma_from_v(cls, v):
        value = 1-numpy.dot(v, v)/cls.VACUUM_LIGHT_SPEED()**2
        if value <= 0:
            raise RuntimeError, "particle velocity >= speed of light"
        return value**(-0.5)

    @classmethod
    def v_from_p(cls, m, p):
        c = cls.VACUUM_LIGHT_SPEED()
        v =  c*p*(numpy.dot(p, p)+c*c*m*m)**(-0.5)
        assert la.norm(v) < c
        return v



class SIUnitsWithNaturalConstants(SIUnits):
    EPSILON0 = 8.8541878176e-12 * SIUnits.F / SIUnits.M
    MU0 = 4*pi*1e-7 * SIUnits.N / SIUnits.A**2

    EL_MASS = 9.10938215e-31 * SIUnits.KG
    EL_CHARGE = 1.602176487e-19 * SIUnits.C
    PROTON_MASS = 1.672621637e-27 * SIUnits.KG

    EV = EL_CHARGE * SIUnits.V





class SIUnitsWithUnityConstants(SIUnits):
    EPSILON0 = 1 * SIUnits.F / SIUnits.M
    MU0 = 1 * SIUnits.N / SIUnits.A**2

    EL_MASS = 1 * SIUnits.KG
    EL_CHARGE = 1 * SIUnits.C
    PROTON_MASS = 1 * SIUnits.KG
