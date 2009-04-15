"""Python interface for depositors"""

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




from pyrticle.deposition import Depositor
import pyrticle._internal as _internal




class ShapeFunctionDepositor(Depositor):
    name = "Shape"

    def initialize(self, method):
        Depositor.initialize(self, method)

        backend_class = getattr(_internal, "InterpolatingDepositor" 
                + method.get_dimensionality_suffix())
        self.backend = backend_class(method.mesh_data)

    def make_state(self, state):
        return self.backend.DepositorState()

    def set_shape_function(self, state, sf):
        Depositor.set_shape_function(self, state, sf)
        self.backend.shape_function = sf

    def note_move(self, state, orig, dest, size):
        pass

    def note_change_size(self, state, count):
        pass





class NormalizedShapeFunctionDepositor(Depositor):
    name = "NormShape"

    def initialize(self, method):
        Depositor.initialize(self, method)

        eg, = method.discretization.element_groups
        ldis = eg.local_discretization

        backend_class = getattr(_internal, "NormalizingInterpolatingDepositor" 
                + method.get_dimensionality_suffix())
        self.backend = backend_class(method.mesh_data, ldis.mass_matrix())

    def add_instrumentation(self, mgr, observer):
        Depositor.add_instrumentation(self, mgr, observer)

        from pyrticle.log import StatsGathererLogQuantity
        mgr.add_quantity(StatsGathererLogQuantity(
            lambda : observer.state.depositor_state.stats.normalization_stats,
            "normshape_norm", "1", 
            "normalization constants applied during deposition"))

        mgr.add_quantity(StatsGathererLogQuantity(
            lambda : observer.state.depositor_state.stats.centroid_distance_stats,
            "normshape_centroid_dist", "m", 
            "distance of shape center from element centroid"))

        mgr.add_quantity(StatsGathererLogQuantity(
            lambda : observer.state.depositor_state.stats.el_per_particle_stats,
            "normshape_el_per_particle", "1", 
            "number of elements per particle"))

    def set_shape_function(self, state, sf):
        Depositor.set_shape_function(self, state, sf)
        self.backend.shape_function = sf

    def note_move(self, state, orig, dest, size):
        pass

    def note_change_size(self, state, count):
        pass
