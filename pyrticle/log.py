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




# -----------------------------------------------------------------------------
class StateObserver:
    def __init__(self, method, maxwell_op):
        self.method = method
        self.maxwell_op = maxwell_op

    def set_fields_and_state(self, fields, state):
        self.fields = fields
        self.state = state

    @property
    def e(self):
        e, h = self.maxwell_op.split_eh(self.fields)
        return e

    @property
    def h(self):
        e, h = self.maxwell_op.split_eh(self.fields)
        return h

    @property
    def have_phi(self):
        from pyrticle.hyperbolic import CleaningMaxwellOperator
        return isinstance(self.maxwell_op, CleaningMaxwellOperator)

    @property
    def phi(self):
        if not self.have_phi:
            raise AttributeError("Observer.phi is only available when "
                    "using hyperbolic cleaning")

        e, h, phi = self.maxwell_op.split_ehphi(self.fields)
        return phi

    @property
    def discr(self):
        return self.method.discretization




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
        data = self.gatherer()
        if data.count():
            return [data.mean(),
                    data.standard_deviation(),
                    data.minimum(),
                    data.maximum()
                    ]
        else:
            return [None, None, None, None]




# Particle quantities ---------------------------------------------------------
class ParticleCount(LogQuantity):
    def __init__(self, observer, name="n_part"):
        LogQuantity.__init__(self, name, "1", "Particle Count")
        self.observer = observer

    def __call__(self):
        return len(self.observer.state)




class ParticleMomentum(MultiLogQuantity):
    def __init__(self, observer, names=None):
        vdim = observer.method.dimensions_velocity
        if names is None:
            names = ["p%s_part" % axis_name(i) for i in range(vdim)]

        MultiLogQuantity.__init__(self, names,
            units=["N*s"] * vdim,
            descriptions=["Particle Momentum"] * vdim)

        self.observer = observer

    def __call__(self):
        from pyrticle._internal import particle_momentum
        return particle_momentum(self.observer.state.particle_state)




class KineticEnergy(LogQuantity):
    def __init__(self, observer, name="W_part"):
        LogQuantity.__init__(self, name, "J", "Kinetic Energy")
        self.observer = observer

    def __call__(self):
        from pyrticle._internal import kinetic_energies
        total_kin_energy = numpy.sum(kinetic_energies(
                self.observer.state.particle_state,
                self.observer.method.units.VACUUM_LIGHT_SPEED()))
        return total_kin_energy




class ParticleCharge(LogQuantity):
    def __init__(self, observer, name="Q_part"):
        LogQuantity.__init__(self, name, "C", "Total particle charge")
        self.observer = observer

    def __call__(self):
        from pyrticle._internal import total_charge
        return total_charge(self.observer.state.particle_state)




def add_particle_quantities(mgr, observer):
    mgr.add_quantity(ParticleCount(observer))
    mgr.add_quantity(ParticleMomentum(observer))
    mgr.add_quantity(KineticEnergy(observer))
    mgr.add_quantity(ParticleCharge(observer))




# EM quantities ---------------------------------------------------------------
class DivergenceEQuantities(MultiLogQuantity):
    def __init__(self, observer, names=["divD", "err_divD_l1"]):
        MultiLogQuantity.__init__(self, names,
                units=["C", "C"],
                descriptions=["Central divergence of D", "L1 Divergence Error"])

        self.observer = observer
        self.discr = self.observer.discr
        from hedge.models.nd_calculus import DivergenceOperator
        self.bound_div_op = DivergenceOperator(self.discr.dimensions,
                observer.maxwell_op.get_eh_subset()[:3]).bind(self.discr)

    def __call__(self):
        rho = self.observer.method.deposit_rho(self.observer.state)
        max_op = self.observer.maxwell_op
        d = max_op.epsilon * self.observer.e
        div_d = self.bound_div_op(d)

        return [
                self.discr.integral(div_d),
                self.discr.integral(numpy.absolute(div_d-rho))
                ]




class DepositedCharge(LogQuantity):
    def __init__(self, observer, name="Q_rec"):
        LogQuantity.__init__(self, name, "C", "Total Charge")
        self.observer = observer

    def __call__(self):
        rho = self.observer.method.deposit_rho(self.observer.state)
        return self.observer.discr.integral(rho)




class FieldCurrent(LogQuantity):
    def __init__(self, observer, direction, tube_length, name="I_field"):
        LogQuantity.__init__(self, name, "A",
                "Integrate current density along %s" % str(list(direction)))
        self.observer = observer
        self.direction = direction
        self.tube_length = tube_length

    def __call__(self):
        return self.observer.discr.integral(
                numpy.dot(self.direction,
                    self.observer.method.deposit_j(self.observer.state)))/self.tube_length




def add_field_quantities(mgr, observer, deposit_interval=5):
    from hedge.log import \
            ElectricFieldEnergy, \
            MagneticFieldEnergy, \
            EMFieldMomentum, \
            EMFieldDivergenceB
    mgr.add_quantity(ElectricFieldEnergy(observer))
    mgr.add_quantity(MagneticFieldEnergy(observer))
    mgr.add_quantity(EMFieldMomentum(observer, observer.method.units.VACUUM_LIGHT_SPEED()))
    mgr.add_quantity(EMFieldDivergenceB(observer.maxwell_op, observer))
    mgr.add_quantity(DivergenceEQuantities(observer), deposit_interval)
    mgr.add_quantity(DepositedCharge(observer), deposit_interval)




# Beam quantities -------------------------------------------------------------
class RMSBeamRadius(LogQuantity):
    def __init__(self, observer, axis, name=None):
        if name is None:
            name = "r%s_rms" % axis_name(axis)

        LogQuantity.__init__(self, name, "m",
                "RMS Beam Radius along %s" % axis_name(axis))
        self.observer = observer
        self.axis = axis

    def __call__(self):
        from pyrticle._internal import rms_beam_size
        return rms_beam_size(self.observer.state.particle_state, self.axis)




class RMSBeamEmittance(LogQuantity):
    def __init__(self, observer, axis, beam_axis, name=None):
        if name is None:
            name = "eps%s_rms" % axis_name(axis)

        LogQuantity.__init__(self, name, "m*rad",
                "RMS Beam Emittance along %s" % axis_name(axis))
        self.observer = observer
        self.axis = axis
        self.beam_axis = beam_axis

    def __call__(self):
        from pyrticle._internal import rms_beam_emittance
        return rms_beam_emittance(self.observer.state.particle_state, self.axis, self.beam_axis)




class RMSBeamEnergySpread(LogQuantity):
    def __init__(self, observer, name="Espread"):
        LogQuantity.__init__(self, name, "J", "RMS Beam Energy Spread")
        self.observer = observer

    def __call__(self):
        from pyrticle._internal import rms_energy_spread
        return rms_energy_spread(self.observer.state.particle_state,
                self.observer.method.units.VACUUM_LIGHT_SPEED())




class ParticleCurrent(LogQuantity):
    def __init__(self, observer, direction, tube_length, name="I_part"):
        LogQuantity.__init__(self, name, "A",
                "Particle Current along %s" % str(list(direction)))
        self.observer = observer
        self.direction = numpy.asarray(direction)
        self.tube_length = tube_length

    def __call__(self):
        from pyrticle._internal import particle_current
        return numpy.dot(
                particle_current(
                    self.observer.state.particle_state,
                    self.observer.method.velocities(self.observer.state),
                    self.tube_length),
                self.direction)




def add_beam_quantities(mgr, observer, axis=0, beam_axis=2):
    mgr.add_quantity(RMSBeamRadius(observer, axis))
    mgr.add_quantity(RMSBeamEmittance(observer, axis, beam_axis))
    mgr.add_quantity(RMSBeamEnergySpread(observer))




def add_currents(mgr, observer, direction, tube_length):
    mgr.add_quantity(FieldCurrent(observer, direction, tube_length))
    mgr.add_quantity(ParticleCurrent(observer, direction, tube_length))




if __name__ == "__main__":
    for i in _join_by_first_of_tuple([[(0,1), (1,2)],[(0,4), (1,5)]]):
        print i
