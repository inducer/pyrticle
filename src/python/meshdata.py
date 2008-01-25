# Pyrticle - Particle in Cell in Python
# Python interface for mesh data
# Copyright (C) 2007 Andreas Kloeckner
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.




import pyrticle._internal as _internal
import pylinear.array as num
import pylinear.computation as comp




MeshData = _internal.MeshData




def _add_mesh_data_methods():
    def fill_from_hedge(self, discr):
        self.discr = discr

        # add periodicity -----------------------------------------------------
        from pyrticle._internal import PeriodicityAxis
        for axis, ((ax_min, ax_max), periodicity_tags) in enumerate(zip(
                zip(*discr.mesh.bounding_box), discr.mesh.periodicity)):
            if periodicity_tags is not None:
                pa = PeriodicityAxis()
                pa.axis = axis
                pa.min = ax_min
                pa.max = ax_max
                pa.width = ax_max-ax_min
                self.periodicities.append(pa)

        # add elements --------------------------------------------------------
        mesh = discr.mesh

        neighbor_map = {}
        for face, (e2, f2) in discr.mesh.both_interfaces():
            neighbor_map[face] = e2.id
        from hedge.mesh import TAG_ALL
        for face in discr.mesh.tag_to_boundary[TAG_ALL]:
            neighbor_map[face] = MeshData.INVALID_ELEMENT

        from pyrticle._internal import ElementInfo

        self.element_info.reserve(len(mesh.elements))
        for i, el in enumerate(mesh.elements):
            ei = ElementInfo()
            ei.id = i
            ei.inverse_map = el.inverse_map
            ei.jacobian = abs(el.map.jacobian)
            ei.start, ei.end = discr.find_el_range(el.id)
            ei.vertices.extend([vi for vi in el.vertex_indices])
            ei.normals.extend(el.face_normals)
            ei.neighbors.extend(
                    list(set(neighbor_map[el,fi] 
                        for fi in xrange(len(el.faces)))))

            self.element_info.append(ei)

        # add vertices --------------------------------------------------------
        vertex_to_element_map = {}

        for el in mesh.elements:
            for vi in el.vertex_indices:
                vertex_to_element_map.setdefault(vi, set()).add(el.id)
                for other_vi in mesh.periodic_opposite_vertices.get(vi, []):
                    vertex_to_element_map.setdefault(other_vi, set()).add(el.id)

        self.vertices.reserve(len(mesh.points))
        self.vertices.extend(mesh.points)

        from pyrticle._internal import UnsignedVector
        self.vertex_adj_elements.reserve(
                2*discr.dimensions*len(mesh.points))
        self.vertex_adj_element_starts.reserve(len(mesh.points))
        self.vertex_adj_element_starts.append(0)

        for vi in xrange(len(mesh.points)):
            self.vertex_adj_elements.extend(vertex_to_element_map[vi])
            self.vertex_adj_element_starts.append(
                    len(self.vertex_adj_elements))

        # add nodes -----------------------------------------------------------
        self.nodes.reserve(len(discr.nodes))
        self.nodes.extend(discr.nodes)

    def find_containing_element(self, point):
        for el in self.discr.mesh.elements:
            if el.contains_point(point):
                return el
        return None

    def min_vertex_distance(self):
        def min_vertex_distance_for_el(el):
            vertices = [self.discr.mesh.points[vi] 
                    for vi in el.vertex_indices]

            return min(min(comp.norm_2(vi-vj)
                    for i, vi in enumerate(vertices)
                    if i != j)
                    for j, vj in enumerate(vertices))

        return min(min_vertex_distance_for_el(el) 
                for el in self.discr.mesh.elements)

    MeshData.fill_from_hedge = fill_from_hedge
    MeshData.min_vertex_distance = min_vertex_distance
_add_mesh_data_methods()
