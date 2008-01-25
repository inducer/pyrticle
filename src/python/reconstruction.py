"""Python interface for reconstructors"""

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




import pylinear.array as num
import pylinear.computation as comp




class ShapeFunctionReconstructor(object):
    name = "Shape"

    def initialize(self, cloud):
        cloud.pic_algorithm.set_radius(0.5*cloud.mesh_data.min_vertex_distance())

    def add_instrumentation(self, mgr):
        pass

    def add_particle_hook(self, pn):
        pass

    def rhs(self):
        return 0

    def add_rhs(self, rhs):
        return 0


class AdvectiveReconstructor(object):
    name = "Advective"

    def initialize(self, cloud):
        self.cloud = cloud
        discr = cloud.mesh_data.discr
        
        assert len(discr.face_groups) == 1
        assert len(discr.element_groups) == 1

        eg = discr.element_groups[0]
        fg, fmm = discr.face_groups[0]
        ldis = eg.local_discretization

        cloud.pic_algorithm.setup_advective_reconstructor(
                len(ldis.face_indices()),
                ldis.node_count(),
                ldis.mass_matrix(),
                ldis.inverse_mass_matrix(),
                fmm,
                fg)

        for i, diffmat in enumerate(ldis.differentiation_matrices()):
            cloud.pic_algorithm.add_local_diff_matrix(i, diffmat)

        self.radius = 0.5*cloud.mesh_data.min_vertex_distance()

    def add_instrumentation(self, mgr):
        pass

    def add_particle_hook(self, pn):
        self.cloud.pic_algorithm.add_advective_particle(self.radius, pn)

    def rhs(self):
        return self.cloud.pic_algorithm.get_advective_particle_rhs(self.cloud.raw_velocities())

    def add_rhs(self, rhs):
        self.cloud.pic_algorithm.apply_advective_particle_rhs(rhs)
