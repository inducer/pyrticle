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
import pylinear.array as num
import pylinear.computation as comp




class Pusher(object):
    def add_instrumentation(self, mgr):
        pass




class MonomialParticlePusher(Pusher):
    name = "Monomial"

    def initialize(self, cloud):
        # add monomial basis data ---------------------------------------------
        from hedge.polynomial import generic_vandermonde
        from pyrticle._internal import MonomialBasisFunction

        for i, eg in enumerate(cloud.mesh_data.discr.element_groups):
            ldis = eg.local_discretization

            from pyrticle._internal import LocalMonomialDiscretization
            lmd = LocalMonomialDiscretization()

            lmd.basis.extend([MonomialBasisFunction(*idx)
                for idx in ldis.node_tuples()])

            mon_vdm_t = generic_vandermonde(ldis.unit_nodes(), lmd.basis).T

            lmd.l_vandermonde_t, \
                    lmd.u_vandermonde_t, perm, sign = comp.lu(mon_vdm_t)
            lmd.p_vandermonde_t = num.permutation_matrix(from_indices=perm)
            cloud.pic_algorithm.local_discretizations.append(lmd)

            cloud.pic_algorithm.ldis_indices.extend([i]*len(eg.members))





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
            return 0
        sd = vis_info[self.vis_field]
        if not sd:
            return 0

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




class AverageParticlePusher(Pusher):
    name = "Average"

    def initialize(self, cloud):
        self.cloud = cloud

        eg, = cloud.mesh_data.discr.element_groups
        ldis = eg.local_discretization

        cloud.pic_algorithm.setup_averaging_particle_pusher(
                ldis.mass_matrix())

    def add_instrumentation(self, mgr):
        mgr.add_quantity(AverageEFieldStdDeviation(self))
        mgr.add_quantity(AverageBFieldStdDeviation(self))

