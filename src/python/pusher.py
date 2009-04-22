"""Python interface for different particle pushers"""

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




import pytools.log
import numpy
import numpy.linalg as la
import pyrticle._internal as _internal




class EmptyPusherState:
    pass




class Pusher(object):
    def __init__(self):
        from pytools.log import IntervalTimer, EventCounter

        self.force_timer = IntervalTimer(
                "t_force",
                "Time spent calculating forces")

    def initialize(self, method):
        self.method = method

    def make_state(self, state):
        return EmptyPusherState()

    def advance_state(self, state):
        return state.pusher_state

    def add_instrumentation(self, mgr, observer):
        mgr.set_constant("pusher", self.__class__.__name__)

        mgr.add_quantity(self.force_timer)

    def upkeep(self, state):
        pass

    def _forces(self, state, velocities, *field_args):
        return self.backend.forces(
                ps=state.particle_state,
                velocities=velocities,
                vis_listener=state.vis_listener,
                *field_args
                )

    def forces(self, state, velocities, *field_args):
        self.force_timer.start()
        forces = self._forces(state, velocities, *field_args)
        self.force_timer.stop()
        return forces

    def note_move(self, state, orig, dest, size):
        pass

    def note_change_size(self, state, count):
        pass




# monomial pusher -------------------------------------------------------------
class MonomialParticlePusher(Pusher):
    def initialize(self, method):
        Pusher.initialize(self, method)

        backend_class = getattr(_internal, "MonomialPusher" 
                + method.get_dimensionality_suffix())
        self.backend = backend_class(method.mesh_data)

        # add monomial basis data ---------------------------------------------
        from hedge.polynomial import generic_vandermonde
        from pyrticle._internal import MonomialBasisFunction, lu

        for i, eg in enumerate(method.mesh_data.discr.element_groups):
            ldis = eg.local_discretization

            from pyrticle._internal import LocalMonomialDiscretization
            lmd = LocalMonomialDiscretization()

            lmd.basis.extend([MonomialBasisFunction(*idx)
                for idx in ldis.node_tuples()])

            mon_vdm_t = numpy.asarray(
                    generic_vandermonde(ldis.unit_nodes(), lmd.basis).T,
                    order="C")

            lmd.lu_vandermonde_t, lmd.lu_piv_vandermonde_t = lu(mon_vdm_t)
            self.backend.local_discretizations.append(lmd)

            self.backend.ldis_indices.extend([i]*len(eg.members))





# average pusher instrumentation ----------------------------------------------
class AverageFieldStdDeviation(pytools.log.LogQuantity):
    def __init__(self, observer, field, name, unit):
        if name is None:
            name = "avg_pt_%s_stddev" % field

        pytools.log.LogQuantity.__init__(self, name, unit, 
                "Average standard deviation of the %s-field "
                "recorded by the average pusher" % field)

        self.observer = observer
        self.vis_field = "pt_%s_stddev" % field

    def __call__(self):
        vis_info = self.observer.state.vis_listener.particle_vis_map

        if self.vis_field not in vis_info:
            return None
        sd = vis_info[self.vis_field]
        if not sd:
            return None

        from pyrticle.tools import NumberShiftableVector
        sd = NumberShiftableVector.unwrap(sd)

        return sd.sum()/len(sd)




class AverageEFieldStdDeviation(AverageFieldStdDeviation):
    def __init__(self, observer, name=None):
        AverageFieldStdDeviation.__init__(self, observer,
                "e", name, "V/m")




class AverageBFieldStdDeviation(AverageFieldStdDeviation):
    def __init__(self, observer, name=None):
        AverageFieldStdDeviation.__init__(self, observer,
                "b", name, "T")




# average pusher --------------------------------------------------------------
class AverageParticlePusher(Pusher):
    def initialize(self, method):
        Pusher.initialize(self, method)

        eg, = method.mesh_data.discr.element_groups
        ldis = eg.local_discretization

        backend_class = getattr(_internal, "AveragingPusher" 
                + method.get_dimensionality_suffix())
        self.backend = backend_class(method.mesh_data, ldis.mass_matrix())

    def make_state(self, state):
        return self.backend.PusherState()

    def add_instrumentation(self, mgr, observer):
        Pusher.add_instrumentation(self, mgr, observer)

        mgr.add_quantity(AverageEFieldStdDeviation(observer))
        mgr.add_quantity(AverageBFieldStdDeviation(observer))

        from pyrticle.log import StatsGathererLogQuantity
        mgr.add_quantity(StatsGathererLogQuantity(
            lambda : observer.state.pusher_state.e_normalization_stats,
            "avgpush_enorm", "1", 
            "normalization constants applied to E-field during averaging particle push"))
        mgr.add_quantity(StatsGathererLogQuantity(
            lambda : observer.state.pusher_state.b_normalization_stats,
            "avgpush_bnorm", "1", 
            "normalization constants applied to B-field during averaging particle push"))

    def _forces(self, state, velocities, *field_args):
        state.pusher_state.e_normalization_stats.reset()
        state.pusher_state.b_normalization_stats.reset()
        return self.backend.forces(
                particle_state=state.particle_state,
                pusher_state=state.pusher_state,
                depositor=self.method.depositor.backend,
                depositor_state=state.depositor_state,
                velocities=velocities,
                vis_listener=state.vis_listener,
                *field_args)
