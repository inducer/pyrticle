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




class Pusher(object):
    def __init__(self):
        from pytools.log import IntervalTimer, EventCounter

        self.force_timer = IntervalTimer(
                "t_force",
                "Time spent calculating forces")

    def initialize(self, cloud):
        self.cloud = cloud

    def add_instrumentation(self, mgr):
        mgr.add_quantity(self.force_timer)

    def forces(self, velocities, verbose_vis, *field_args):
        self.force_timer.start()
        forces = self.cloud.pic_algorithm.forces(
                velocities=velocities,
                verbose_vis=verbose_vis,
                *field_args
                )
        self.force_timer.stop()
        return forces




# monomial pusher -------------------------------------------------------------
class MonomialParticlePusher(Pusher):
    name = "Monomial"

    def initialize(self, cloud):
        Pusher.initialize(self, cloud)

        # add monomial basis data ---------------------------------------------
        from hedge.polynomial import generic_vandermonde
        from pyrticle._internal import MonomialBasisFunction, lu

        for i, eg in enumerate(cloud.mesh_data.discr.element_groups):
            ldis = eg.local_discretization

            from pyrticle._internal import LocalMonomialDiscretization
            lmd = LocalMonomialDiscretization()

            lmd.basis.extend([MonomialBasisFunction(*idx)
                for idx in ldis.node_tuples()])

            mon_vdm_t = numpy.asarray(
                    generic_vandermonde(ldis.unit_nodes(), lmd.basis).T,
                    order="C")

            lvt, uvt, perm = lu(mon_vdm_t)
            from pyublas import permutation_matrix
            lmd.l_vandermonde_t = numpy.asarray(lvt, order="C")
            lmd.u_vandermonde_t = numpy.asarray(uvt, order="C")
            lmd.p_vandermonde_t = permutation_matrix(from_indices=perm)
            cloud.pic_algorithm.local_discretizations.append(lmd)

            cloud.pic_algorithm.ldis_indices.extend([i]*len(eg.members))





# average pusher instrumentation ----------------------------------------------
class AverageFieldStdDeviation(pytools.log.LogQuantity):
    def __init__(self, pusher, field, name, unit):
        if name is None:
            name = "avg_pt_%s_stddev" % field

        pytools.log.LogQuantity.__init__(self, name, unit, 
                "Average standard deviation of the %s-field "
                "recorded by the average pusher" % field)

        self.pusher = pusher
        self.vis_field = "pt_%s_stddev" % field

    def __call__(self):
        vis_info = self.pusher.cloud.vis_info
        if self.vis_field not in vis_info:
            return None
        sd = vis_info[self.vis_field]
        if not sd:
            return None

        from pyrticle.tools import NumberShiftableVector
        sd = NumberShiftableVector.unwrap(sd)

        return sd.sum()/len(sd)




class AverageEFieldStdDeviation(AverageFieldStdDeviation):
    def __init__(self, pusher, name=None):
        AverageFieldStdDeviation.__init__(self, pusher,
                "e", name, "V/m")




class AverageBFieldStdDeviation(AverageFieldStdDeviation):
    def __init__(self, pusher, name=None):
        AverageFieldStdDeviation.__init__(self, pusher,
                "b", name, "T")




# average pusher --------------------------------------------------------------
class AverageParticlePusher(Pusher):
    name = "Average"

    def initialize(self, cloud):
        Pusher.initialize(self, cloud)

        eg, = cloud.mesh_data.discr.element_groups
        ldis = eg.local_discretization

        cloud.pic_algorithm.setup_averaging_particle_pusher(
                ldis.mass_matrix())

    def add_instrumentation(self, mgr):
        Pusher.add_instrumentation(self, mgr)

        mgr.add_quantity(AverageEFieldStdDeviation(self))
        mgr.add_quantity(AverageBFieldStdDeviation(self))

        from pyrticle.log import StatsGathererLogQuantity
        mgr.add_quantity(StatsGathererLogQuantity(
            self.cloud.pic_algorithm.e_normalization_stats,
            "avgpush_enorm", "1", 
            "normalization constants applied to E-field during averaging particle push"))
        mgr.add_quantity(StatsGathererLogQuantity(
            self.cloud.pic_algorithm.b_normalization_stats,
            "avgpush_bnorm", "1", 
            "normalization constants applied to B-field during averaging particle push"))

    def forces(self, velocities, verbose_vis, *field_args):
        self.cloud.pic_algorithm.e_normalization_stats.reset()
        self.cloud.pic_algorithm.b_normalization_stats.reset()
        return Pusher.forces(self, velocities, verbose_vis, *field_args)

    

