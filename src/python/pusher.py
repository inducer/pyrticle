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




import pylinear.array as num
import pylinear.computation as comp




class MonomialParticlePusher:
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





class AverageParticlePusher:
    name = "Average"

    def initialize(self, cloud):
        eg, = cloud.mesh_data.discr.element_groups
        ldis = eg.local_discretization

        cloud.pic_algorithm.setup_averaging_particle_pusher(
                ldis.mass_matrix())
