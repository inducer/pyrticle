"""Data logging"""

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




from __future__ import division
from pytools.log import LogQuantity
import pylinear.array as num
import pylinear.computation as comp




# Particle quantities ---------------------------------------------------------
class ParticleCount(LogQuantity):
    def __init__(self, cloud, name="n_part"):
        LogQuantity.__init__(self, name, "1", "Particle Count")
        self.cloud = cloud

    def __call__(self):
        return len(self.cloud)




class ParticleMomentum(LogQuantity):
    def __init__(self, cloud, name="p_part"):
        LogQuantity.__init__(self, name, "N*s", "Particle Momentum")
        self.cloud = cloud

    def __call__(self):
        total_momentum = num.array(
                [num.sum(p_i) for p_i in 
                    self.cloud.momenta.get_alist_of_components()])
        return comp.norm_2(total_momentum)




class KineticEnergy(LogQuantity):
    def __init__(self, cloud, name="W_part"):
        LogQuantity.__init__(self, name, "J", "Kinetic Energy")
        self.cloud = cloud

    def __call__(self):
        total_kin_energy = num.sum(
                self.cloud.pic_algorithm.kinetic_energies())
        return total_kin_energy




class ParticleCharge(LogQuantity):
    def __init__(self, cloud, name="Q_part"):
        LogQuantity.__init__(self, name, "C", "Total particle charge")
        self.cloud = cloud

    def __call__(self):
        return self.cloud.pic_algorithm.total_charge()




def add_particle_quantities(mgr, cloud):
    mgr.add_quantity(ParticleCount(cloud))
    mgr.add_quantity(ParticleMomentum(cloud))
    mgr.add_quantity(KineticEnergy(cloud))
    mgr.add_quantity(ParticleCharge(cloud))




# EM quantities ---------------------------------------------------------------
class FieldEnergy(LogQuantity):
    def __init__(self, fields, name="W_field"):
        LogQuantity.__init__(self, name, "J", "Field Energy")
        self.fields = fields

    def __call__(self):
        max_op = self.fields.maxwell_op

        e = self.fields.e
        h = self.fields.h
        d = max_op.epsilon * e
        b = max_op.mu * h

        from hedge.tools import dot
        energy_density = 1/2*(
                dot(e, d, num.multiply) 
                + dot(h, b, num.multiply))

        from hedge.discretization import integral
        return integral(max_op.discr, energy_density)




class FieldMomentum(LogQuantity):
    def __init__(self, fields, name="p_field"):
        LogQuantity.__init__(self, name, "N*s", "Field Momentum")
        self.fields = fields

    def __call__(self):
        max_op = self.fields.maxwell_op

        mu0 = self.fields.cloud.units.MU0
        c = self.fields.cloud.units.VACUUM_LIGHT_SPEED

        e = self.fields.e
        h = self.fields.h

        poynting_s = max_op.h_cross(1/mu0*e, h, 
                three_mult=lambda lc, x, y: lc*num.multiply(x,y))

        norm_poynting_s = num.sqrt(
                sum(num.power(si, 2) for si in poynting_s))
        momentum_density = norm_poynting_s/c**2

        from hedge.discretization import integral
        return integral(max_op.discr, momentum_density)




class DivergenceError(LogQuantity):
    def __init__(self, f_and_c, name="div_err"):
        LogQuantity.__init__(self, name, "C", "Divergence Error")
        self.f_and_c = f_and_c
        self.discr = self.f_and_c.maxwell_op.discr
        from hedge.operators import DivergenceOperator
        self.div_op = DivergenceOperator(self.discr)

    def __call__(self):
        rho = self.f_and_c.cloud.reconstruct_rho()
        max_op = self.f_and_c.maxwell_op
        d = max_op.epsilon * self.f_and_c.e
        div_d = self.div_op(d)
        
        from hedge.discretization import integral
        return integral(self.discr, num.absolute(div_d-rho))




class ReconstructedCharge(LogQuantity):
    def __init__(self, f_and_c, name="Q_rec"):
        LogQuantity.__init__(self, name, "C", "Total Charge")
        self.f_and_c = f_and_c

    def __call__(self):
        from hedge.discretization import integral
        return integral(self.f_and_c.maxwell_op.discr, 
                self.f_and_c.cloud.reconstruct_rho())




class NetCurrent(LogQuantity):
    def __init__(self, f_and_c, direction, name="I"):
        LogQuantity.__init__(self, name, "A", 
                "Net Current along %s" % str(list(direction)))
        self.f_and_c = f_and_c
        self.direction = direction

    def __call__(self):
        from hedge.discretization import integral
        from hedge.tools import dot
        return integral(self.f_and_c.maxwell_op.discr, 
                dot(self.direction,
                    self.f_and_c.cloud.reconstruct_j()))




def add_field_quantities(mgr, f_and_c):
    mgr.add_quantity(FieldEnergy(f_and_c))
    mgr.add_quantity(FieldMomentum(f_and_c))
    mgr.add_quantity(DivergenceError(f_and_c))
    mgr.add_quantity(ReconstructedCharge(f_and_c))




# Beam quantities -------------------------------------------------------------
def _axis_name(axis):
    if axis == 0: return "X"
    elif axis == 1: return "Y"
    elif axis == 2: return "Z"
    else: return "axis %d" % axis




class RMSBeamRadius(LogQuantity):
    def __init__(self, cloud, axis, name=None):
        if name is None:
            name = "r%s" % _axis_name(axis)

        LogQuantity.__init__(self, name, "m", 
                "RMS Beam Radius along %s" % _axis_name(axis))
        self.cloud = cloud
        self.axis = axis

    def __call__(self):
        from pyrticle.beam import calculate_rms_beam_size
        return calculate_rms_beam_size(self.cloud, self.axis)
    


class RMSBeamEmittance(LogQuantity):
    def __init__(self, cloud, axis, beam_axis, name=None):
        if name is None:
            name = "eps%s" % _axis_name(axis)

        LogQuantity.__init__(self, name, "m*rad", 
                "RMS Beam Emittance along %s" % _axis_name(axis))
        self.cloud = cloud
        self.axis = axis
        self.beam_axis = beam_axis

    def __call__(self):
        from pyrticle.beam import calculate_rms_emittance
        return calculate_rms_emittance(self.cloud, self.axis, self.beam_axis)
    



class RMSBeamEnergySpread(LogQuantity):
    def __init__(self, cloud, name="Espread"):
        LogQuantity.__init__(self, name, "m*rad", 
                "RMS Beam Energy Spread")
        self.cloud = cloud

    def __call__(self):
        from pyrticle.beam import calculate_rms_energy_spread
        return calculate_rms_energy_spread(self.cloud)




def add_beam_quantities(mgr, f_and_c, axis=0, beam_axis=2):
    current_dir = num.unit_vector(f_and_c.cloud.dimensions_velocity, beam_axis)

    mgr.add_quantity(RMSBeamRadius(f_and_c.cloud, axis))
    mgr.add_quantity(RMSBeamEmittance(f_and_c.cloud, axis, beam_axis))
    mgr.add_quantity(RMSBeamEnergySpread(f_and_c.cloud))
    mgr.add_quantity(NetCurrent(f_and_c, current_dir))




if __name__ == "__main__":
    for i in _join_by_first_of_tuple([[(0,1), (1,2)],[(0,4), (1,5)]]):
        print i
