"""Kapchinskij-Vladimirskij beam physics"""

from __future__ import division

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




import numpy
import numpy.linalg as la
from pytools.log import SimulationLogQuantity




def uniform_on_unit_sphere(dim):
    from random import gauss

    # cf.
    # http://www-alg.ist.hokudai.ac.jp/~jan/randsphere.pdf
    # Algorith due to Knuth

    pt = numpy.array([gauss(0,1) for i in range(dim)])
    n2 = la.norm(pt)
    return pt/n2





    

class KVZIntervalBeam:
    def __init__(self, units, nparticles, p_charge, p_mass,
            radii, emittances, beta, z_length, z_pos):
        """Construct a beam that is KV-distributed in (x,y)
        and uniform over an interval along z.

        @par radii: total (100%) radii
        @par emittances: total (100%) emittances
        """

        assert len(radii) == len(emittances)

        self.units = units

        self.nparticles = nparticles

        self.p_charge = p_charge
        self.p_mass = p_mass

        self.radii = radii
        self.emittances = emittances

        self.beta = beta
        self.gamma = (1-beta**2)**(-0.5)

        self.z_length = z_length
        self.z_pos = z_pos

    @property
    def rms_radii(self):
        return [r/2 for r in self.radii]

    @property
    def rms_emittances(self):
        return [eps/4 for eps in self.emittances]

    def make_particle(self, vz, dim_x, dim_p):
        """Return (position, velocity) for a random particle
        according to a Kapchinskij-Vladimirskij distribution.
        """

        s = uniform_on_unit_sphere(len(self.radii) + len(self.emittances))
        x = [x_i*r_i for x_i, r_i in zip(s[:len(self.radii)], self.radii)]
        # xp like xprime
        xp = [s_i/r_i*eps_i 
                for s_i, r_i, eps_i in 
                zip(s[len(self.radii):], self.radii, self.emittances)]

        one = sum(x_i**2/r_i**2 for x_i, r_i in zip(x, self.radii)) + \
                sum(xp_i**2*r_i**2/eps_i**2 
                for xp_i, r_i, epsi in zip(xp, self.radii, self.emittances))
        assert abs(one-1) < 1e-15

        while len(x) < dim_x:
            x.append(0)
        while len(xp) < dim_p:
            xp.append(0)

        z = numpy.array([0,0,1])

        return (numpy.array(x),
                z*vz + numpy.array([xp_i*vz for xp_i in xp]))

    def add_to(self, cloud, discr):
        from random import uniform
        from math import sqrt

        positions = []
        velocities = []

        bbox_min, bbox_max = discr.mesh.bounding_box()
        center = (bbox_min+bbox_max)/2
        center[2] = 0
        size = bbox_max-bbox_min

        vz = self.beta*self.units.VACUUM_LIGHT_SPEED
        z = numpy.array([0,0,1])

        for i in range(self.nparticles):
            pos, v = self.make_particle(
                    vz=self.beta*self.units.VACUUM_LIGHT_SPEED,
                    dim_x=cloud.dimensions_pos, dim_p=cloud.dimensions_velocity)

            my_beta = la.norm(v)/self.units.VACUUM_LIGHT_SPEED
            
            if abs(self.beta) > 1e-4:
                assert abs(self.beta - my_beta)/self.beta < 1e-4

            positions.append(center
                    +pos
                    +z*(self.z_pos+uniform(-self.z_length, self.z_length)/2))
            velocities.append(v)

        cloud.add_particles(positions, velocities, 
                self.p_charge, self.p_mass)

    def analytic_rho(self, discr):
        from pytools import product

        z_min = self.z_pos - self.z_length/2
        z_max = self.z_pos + self.z_length/2

        def distrib(x):
            if (z_min <= x[2] <= z_max and
                    sum((xi/ri)**2 for xi, ri in zip(x, self.radii)) <= 1):
                return 1
            else:
                return 0

        from math import pi
        from pyrticle._internal import gamma
        n = len(self.radii)
        distr_vol = (2*pi)**(n/2) / (gamma(n/2)*n) * product(self.radii) * self.z_length

        total_charge = self.p_charge * self.nparticles 
        unscaled_rho = discr.interpolate_volume_function(distrib)
        rho = total_charge * unscaled_rho / distr_vol

        # check for correctness
        from hedge.discretization import integral
        int_rho = integral(discr, rho)
        vol = integral(discr, unscaled_rho)
        rel_err = (int_rho-total_charge)/total_charge
        assert rel_err < 0.2

        return rho

    def get_total_space_charge_parameter(self):
        from math import pi

        # see http://en.wikipedia.org/wiki/Classical_electron_radius
        r0 = 1/(4*pi*self.units.EPSILON0)*( 
                (self.units.EL_CHARGE**2)
                /
                (self.units.EL_MASS*self.units.VACUUM_LIGHT_SPEED**2))

        total_charge = abs(self.p_charge*self.nparticles)

        lambda_ = total_charge/(self.z_length*self.units.EL_CHARGE)

        #print "total_charge", total_charge
        #print "z_length", self.z_length
        #print "lambda", lambda_

        #print  "beta", self.beta
        #print  "gamma", self.gamma
        #print "xi", xi

        # factor of 2 here is uncertain
        # from S.Y.Lee, Accelerator Physics, p. 68
        # 2nd ed. 
        # (2.140), space charge term (2.136)

        xi = 4*((lambda_ * r0) / (self.beta**2 * self.gamma**3))

        return xi

    def get_rms_space_charge_parameter(self):
        # by rms scaling analysis on the KV ODE
        return self.get_total_space_charge_parameter()/4

    def get_chargeless_rms_predictor(self, axis):
        return ChargelessKVRadiusPredictor(
                self.rms_radii[axis], self.rms_emittances[axis])

    def get_rms_predictor(self, axis):
        return KVRadiusPredictor(
                self.rms_radii[axis], self.rms_emittances[axis],
                xi=self.get_rms_space_charge_parameter())

    def get_total_predictor(self, axis):
        return KVRadiusPredictor(
                self.radii[axis], self.emittances[axis],
                xi=self.get_total_space_charge_parameter())




class ODEDefinedFunction:
    def __init__(self, t0, y0, dt):
        self.t = [t0]
        self.y = [y0]
        self.dt = dt

        from hedge.timestep import RK4TimeStepper
        self.forward_stepper = RK4TimeStepper()
        self.backward_stepper = RK4TimeStepper()

    def __call__(self, t):
        def copy_if_necessary(x):
            try:
                return x[:]
            except TypeError:
                return x

        if t < self.t[0]:
            steps = int((self.t[0]-t)/self.dt)+1
            t_list = [self.t[0]]
            y_list = [self.y[0]]
            for n in range(steps):
                y_list.append(self.backward_stepper(
                    copy_if_necessary(y_list[-1]), 
                    t_list[-1], -self.dt, self.rhs))
                t_list.append(t_list[-1]-self.dt)

            self.t = t_list[:0:-1] + self.t
            self.y = y_list[:0:-1] + self.y
        elif t >= self.t[-1]:
            steps = int((t-self.t[-1])/self.dt)+1
            t_list = [self.t[-1]]
            y_list = [self.y[-1]]
            for n in range(steps):
                y_list.append(self.forward_stepper(
                    copy_if_necessary(y_list[-1]), 
                    t_list[-1], self.dt, self.rhs))
                t_list.append(t_list[-1]+self.dt)

            self.t = self.t + t_list[1:]
            self.y = self.y + y_list[1:]

        from bisect import bisect_right
        below_idx = bisect_right(self.t, t)-1
        assert below_idx >= 0
        above_idx = below_idx + 1

        assert above_idx < len(self.t)
        assert self.t[below_idx] <= t <= self.t[above_idx]

        # FIXME linear interpolation, bad
        slope = ((self.y[above_idx]-self.y[below_idx]) 
                /
                (self.t[above_idx]-self.t[below_idx]))

        return self.y[below_idx] + (t-self.t[below_idx]) * slope

    def rhs(self, t, y):
        raise NotImplementedError





class ChargelessKVRadiusPredictor:
    def __init__(self, a0, eps):
        self.a0 = a0
        self.eps = eps

    def __call__(self, s):
        from math import sqrt
        return sqrt(self.a0**2+(self.eps/self.a0)**2 * s**2)




class KVRadiusPredictor(ODEDefinedFunction):
    """Implement equation (1.74) in Alex Chao's book.

    See equation (1.65) for the definition of M{xi}.
    M{Q} is the number of electrons in the beam
    """
    def __init__(self, a0, eps, eB_2E=0, xi=0, dt=1e-4):
        ODEDefinedFunction.__init__(self, 0, numpy.array([a0, 0]), 
                dt=dt*(a0**4/eps**2)**2)
        self.eps = eps
        self.xi = xi
        self.eB_2E = eB_2E

    def rhs(self, t, y):
        a = y[0]
        aprime = y[1]
        return numpy.array([
            aprime, 
            - self.eB_2E**2 * a
            + self.eps**2/a**3
            + self.xi/(2*a)
            ])

    def __call__(self, t):
        return ODEDefinedFunction.__call__(self, t)[0]




class KVPredictedRadius(SimulationLogQuantity):
    def __init__(self, dt, beam_v, predictor, suffix, name=None):
        if name is None:
            name = "r%s_theory" % suffix

        SimulationLogQuantity.__init__(self, dt, name, "m", 
                "Theoretical RMS Beam Radius")

        self.beam_v = beam_v
        self.predictor = predictor
        self.t = 0

    def __call__(self):
        s = self.beam_v * self.t
        self.t += self.dt
        return self.predictor(s)




if __name__ == "__main__":
    class Sin(ODEDefinedFunction):
        def __init__(self):
            ODEDefinedFunction.__init__(self, 0, 
                    numpy.array([0,1]), 1/7*1e-2)

        def rhs(self, t, y):
            return numpy.array([y[1], -y[0]])

        def __call__(self, t):
            return ODEDefinedFunction.__call__(self, t)[0]

    from math import pi
    s = Sin()
    assert abs(s(-pi)) < 2e-3
    assert abs(s(pi)) < 2e-3
    assert abs(s(-pi/2)+1) < 2e-3

    kv_env_exact = ChargelessKVRadiusPredictor(2.5e-3, 5e-6)
    kv_env_num = KVRadiusPredictor(2.5e-3, 5e-6)

    from hedge.tools import plot_1d
    steps = 50
    for i in range(steps):
        s = kv_env_num.dt/7*i
        
        a_exact = kv_env_exact(s)
        a_num = kv_env_num(s)
        assert abs(a_exact-a_num)/a_exact < 1e-3
