"""Data logging"""

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




from pytools.log import LogQuantity, MultiLogQuantity
from hedge.log import axis_name
import numpy




# generic StatsGatherer logging -----------------------------------------------
class StatsGathererLogQuantity(MultiLogQuantity):
    def __init__(self, gatherer, basename, unit, description):
        MultiLogQuantity.__init__(self, 
                names=[
                    "%s_mean" % basename,
                    "%s_stddev" % basename,
                    "%s_min" % basename,
                    "%s_max" % basename,
                    ],
                    units=4 * [unit], 
                    descriptions=[
                        "Mean of %s" % description,
                        "Standard deviation of %s" % description,
                        "Minimum of %s" % description,
                        "Maximum of %s" % description,
                        ])

        self.gatherer = gatherer

    def __call__(self):
        if self.gatherer.count():
            return [self.gatherer.mean(), 
                    self.gatherer.standard_deviation(),
                    self.gatherer.minimum(),
                    self.gatherer.maximum()
                    ]
        else:
            return [None, None, None, None]




# Particle quantities ---------------------------------------------------------
class ParticleCount(LogQuantity):
    def __init__(self, cloud, name="n_part"):
        LogQuantity.__init__(self, name, "1", "Particle Count")
        self.cloud = cloud

    def __call__(self):
        return len(self.cloud)




class ParticleMomentum(MultiLogQuantity):
    def __init__(self, cloud, names=None):
        vdim = cloud.dimensions_velocity
        if names is None:
            names = ["p%s_part" % axis_name(i) for i in range(vdim)]

        MultiLogQuantity.__init__(self, names, 
            units=["N*s"] * vdim, 
            descriptions=["Particle Momentum"] * vdim)

        self.cloud = cloud

    def __call__(self):
        from pyrticle._internal import particle_momentum
        return particle_momentum(self.cloud.pic_algorithm)




class KineticEnergy(LogQuantity):
    def __init__(self, cloud, name="W_part"):
        LogQuantity.__init__(self, name, "J", "Kinetic Energy")
        self.cloud = cloud

    def __call__(self):
        from pyrticle._internal import kinetic_energies
        total_kin_energy = numpy.sum(kinetic_energies(
                self.cloud.pic_algorithm))
        return total_kin_energy




class ParticleCharge(LogQuantity):
    def __init__(self, cloud, name="Q_part"):
        LogQuantity.__init__(self, name, "C", "Total particle charge")
        self.cloud = cloud

    def __call__(self):
        from pyrticle._internal import total_charge
        return total_charge(self.cloud.pic_algorithm)




def add_particle_quantities(mgr, cloud):
    mgr.add_quantity(ParticleCount(cloud))
    mgr.add_quantity(ParticleMomentum(cloud))
    mgr.add_quantity(KineticEnergy(cloud))
    mgr.add_quantity(ParticleCharge(cloud))




# EM quantities ---------------------------------------------------------------
class DivergenceEQuantities(MultiLogQuantity):
    def __init__(self, f_and_c, names=["divD", "err_divD_l1"]):
        MultiLogQuantity.__init__(self, names, 
                units=["C", "C"], 
                descriptions=["Central divergence of D", "L1 Divergence Error"])

        self.f_and_c = f_and_c
        self.discr = self.f_and_c.maxwell_op.discr
        from hedge.pde import DivergenceOperator
        self.div_op = DivergenceOperator(self.discr,
                f_and_c.maxwell_op.get_eh_subset()[:3])

    def __call__(self):
        rho = self.f_and_c.cloud.reconstruct_rho()
        max_op = self.f_and_c.maxwell_op
        d = max_op.epsilon * self.f_and_c.e
        div_d = self.div_op(d)
        
        return [
                self.discr.integral(div_d),
                self.discr.integral(numpy.absolute(div_d-rho))
                ]




class ReconstructedCharge(LogQuantity):
    def __init__(self, f_and_c, name="Q_rec"):
        LogQuantity.__init__(self, name, "C", "Total Charge")
        self.f_and_c = f_and_c

    def __call__(self):
        return self.f_and_c.maxwell_op.discr.integral(
                self.f_and_c.cloud.reconstruct_rho())




class FieldCurrent(LogQuantity):
    def __init__(self, f_and_c, direction, tube_length, name="I_field"):
        LogQuantity.__init__(self, name, "A", 
                "Integrate current density along %s" % str(list(direction)))
        self.f_and_c = f_and_c
        self.direction = direction
        self.tube_length = tube_length

    def __call__(self):
        return self.f_and_c.maxwell_op.discr.integral(
                numpy.dot(self.direction,
                    self.f_and_c.cloud.reconstruct_j()))/self.tube_length




def add_field_quantities(mgr, f_and_c, reconstruct_interval=5):
    from hedge.log import \
            EMFieldEnergy, \
            EMFieldMomentum, \
            EMFieldDivergenceB
    mgr.add_quantity(EMFieldEnergy(f_and_c))
    mgr.add_quantity(EMFieldMomentum(f_and_c, f_and_c.cloud.units.VACUUM_LIGHT_SPEED))
    mgr.add_quantity(EMFieldDivergenceB(f_and_c.maxwell_op, f_and_c))
    mgr.add_quantity(DivergenceEQuantities(f_and_c), reconstruct_interval)
    mgr.add_quantity(ReconstructedCharge(f_and_c), reconstruct_interval)




# Beam quantities -------------------------------------------------------------
class RMSBeamRadius(LogQuantity):
    def __init__(self, cloud, axis, name=None):
        if name is None:
            name = "r%s_rms" % axis_name(axis)

        LogQuantity.__init__(self, name, "m", 
                "RMS Beam Radius along %s" % axis_name(axis))
        self.cloud = cloud
        self.axis = axis

    def __call__(self):
        from pyrticle.beam import rms_beam_size
        return rms_beam_size(self.cloud, self.axis)




class RMSBeamEmittance(LogQuantity):
    def __init__(self, cloud, axis, beam_axis, name=None):
        if name is None:
            name = "eps%s_rms" % axis_name(axis)

        LogQuantity.__init__(self, name, "m*rad", 
                "RMS Beam Emittance along %s" % axis_name(axis))
        self.cloud = cloud
        self.axis = axis
        self.beam_axis = beam_axis

    def __call__(self):
        from pyrticle.beam import rms_emittance
        return rms_emittance(self.cloud, self.axis, self.beam_axis)
    



class RMSBeamEnergySpread(LogQuantity):
    def __init__(self, cloud, name="Espread"):
        LogQuantity.__init__(self, name, "J", "RMS Beam Energy Spread")
        self.cloud = cloud

    def __call__(self):
        from pyrticle.beam import rms_energy_spread
        return rms_energy_spread(self.cloud)




class ParticleCurrent(LogQuantity):
    def __init__(self, cloud, direction, tube_length, name="I_part"):
        LogQuantity.__init__(self, name, "A",
                "Particle Current along %s" % str(list(direction)))
        self.cloud = cloud
        self.direction = numpy.asarray(direction)
        self.tube_length = tube_length

    def __call__(self):
        from pyrticle._internal import particle_current
        return numpy.dot(
                particle_current(
                    self.cloud.pic_algorithm, 
                    self.cloud.velocities(),
                    self.tube_length),
                self.direction)




def add_beam_quantities(mgr, cloud, axis=0, beam_axis=2):
    mgr.add_quantity(RMSBeamRadius(cloud, axis))
    mgr.add_quantity(RMSBeamEmittance(cloud, axis, beam_axis))
    mgr.add_quantity(RMSBeamEnergySpread(cloud))




def add_currents(mgr, f_and_c, direction, tube_length):
    mgr.add_quantity(FieldCurrent(f_and_c, direction, tube_length))
    mgr.add_quantity(ParticleCurrent(f_and_c.cloud, direction, tube_length))




if __name__ == "__main__":
    for i in _join_by_first_of_tuple([[(0,1), (1,2)],[(0,4), (1,5)]]):
        print i
