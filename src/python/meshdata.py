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

        # add periodicity
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

    def find_containing_element(self, point):
        for el in self.discr.mesh.elements:
            if el.contains_point(point):
                return el
        return None

    def add_elements(self, discr):
        ldis_indices = self.add_local_discretizations(discr)

        mesh = discr.mesh

        neighbor_map = {}
        for face, (e2, f2) in discr.mesh.both_interfaces():
            neighbor_map[face] = e2.id
        from hedge.mesh import TAG_ALL
        for face in discr.mesh.tag_to_boundary[TAG_ALL]:
            neighbor_map[face] = MeshData.INVALID_ELEMENT

        for el in mesh.elements:
            (estart, eend), ldis = discr.find_el_data(el.id)

            self.add_element(
                    el.inverse_map,
                    ldis_indices[ldis],
                    estart, eend,
                    [vi for vi in el.vertex_indices],
                    el.face_normals,
                    set(neighbor_map[el,fi] for fi in xrange(len(el.faces)))
                    )

    def add_vertices(self, discr):
        mesh = discr.mesh

        vertex_to_element_map = {}
        for el in mesh.elements:
            for vi in el.vertex_indices:
                vertex_to_element_map.setdefault(vi, set()).add(el.id)
                for other_vi in mesh.periodic_opposite_vertices.get(vi, []):
                    vertex_to_element_map.setdefault(other_vi, set()).add(el.id)

        for vi in xrange(len(mesh.points)):
            self.add_vertex(vi, 
                    mesh.points[vi],
                    vertex_to_element_map[vi])


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
    MeshData.add_elements = add_elements
    MeshData.add_vertices = add_vertices
    MeshData.min_vertex_distance = min_vertex_distance
_add_mesh_data_methods()
