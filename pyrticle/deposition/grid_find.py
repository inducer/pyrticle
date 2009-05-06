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
from pyrticle.deposition.grid_base import GridVisualizer, SingleBrickGenerator
import pyrticle._internal as _internal
import numpy




# grid find deposition ----------------------------------------------------
class GridFindDepositor(Depositor, GridVisualizer):
    name = "GridFind"
    iterator_type = _internal.BrickIterator

    def __init__(self, brick_generator=None):
        Depositor.__init__(self)
        self.brick_generator = brick_generator

    def initialize(self, method):
        Depositor.initialize(self, method)

        backend_class = getattr(_internal, "GridFindDepositor" 
                + method.get_dimensionality_suffix())
        backend = self.backend = backend_class(method.mesh_data)

        discr = method.discretization

        grid_node_num_to_nodes = {}

        if self.brick_generator is None:
            bbox_min, bbox_max = discr.mesh.bounding_box()
            max_bbox_size = max(bbox_max-bbox_min)
            self.brick_generator = SingleBrickGenerator(
                    mesh_margin=1e-3*max_bbox_size,
                    overresolve=0.2)

        from pyrticle._internal import Brick
        for i, (stepwidths, origin, dims) in enumerate(
                self.brick_generator(discr)):
            backend.bricks.append(
                    Brick(i, backend.grid_node_count(), stepwidths, origin, dims))

        from pyrticle._internal import BoxFloat
        for eg in discr.element_groups:
            ldis = eg.local_discretization

            for el in eg.members:
                el_bbox = BoxFloat(*el.bounding_box(discr.mesh.points))
                el_slice = discr.find_el_range(el.id)

                for brk in backend.bricks:
                    if brk.bounding_box().intersect(el_bbox).is_empty():
                        continue

                    for node_num in range(el_slice.start, el_slice.stop):
                        try:
                            cell_number = brk.which_cell(
                                        discr.nodes[node_num])
                        except ValueError:
                            pass
                        else:
                            grid_node_num_to_nodes.setdefault(
                                    brk.index(cell_number), []).append(node_num)
                        
        from pytools import flatten
        unassigned_nodes = (set(xrange(len(discr))) 
                - set(flatten(grid_node_num_to_nodes.itervalues())))

        if unassigned_nodes:
            raise RuntimeError("dep_grid_find: unassigned mesh nodes found. "
                    "you should specify a mesh_margin when generating "
                    "bricks")

        usecounts = numpy.zeros(
                (backend.grid_node_count(),))
        for gnn in xrange(backend.grid_node_count()):
            grid_nodes = grid_node_num_to_nodes.get(gnn, [])
            usecounts[gnn] = len(grid_nodes)
            backend.node_number_list_starts.append(len(backend.node_number_lists))
            backend.node_number_lists.extend(grid_nodes)
        backend.node_number_list_starts.append(len(backend.node_number_lists))

        if "depositor" in method.debug:
            from hedge.visualization import SiloVisualizer
            vis = SiloVisualizer(discr)
            visf = vis.make_file("grid-find-debug")
            vis.add_data(visf, [])
            self.visualize_grid_quantities(visf, [
                    ("usecounts", usecounts)
                    ])
            visf.close()

            if "interactive" in cloud.debug:
                from matplotlib.pylab import hist, show
                hist(usecounts, bins=20)
                show()

    def set_shape_function(self, state, sf):
        Depositor.set_shape_function(self, state, sf)
        self.backend.shape_function = sf

    def make_state(self, state):
        return self.backend.DepositorState()

    def note_move(self, state, orig, dest, size):
        state.depositor_state.note_move(orig, dest, size)

    def note_change_size(self, state, count):
        state.depositor_state.note_change_size(count)

    def _deposit_densities(self, state, velocities, pslice):
        return self.backend.deposit_densities(
                state.depositor_state,
                state.particle_state,
                velocities, pslice)

    def _deposit_j(self, state, velocities, pslice):
        return self.backend.deposit_j(
                state.depositor_state,
                state.particle_state,
                velocities, pslice)

    def _deposit_rho(self, state, pslice):
        return self.backend.deposit_rho(
            state.depositor_state,
            state.particle_state,
            pslice)

