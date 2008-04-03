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
import numpy
import numpy.linalg as la
from pytools import monkeypatch_class




MeshData = _internal.MeshData





class MeshData(_internal.MeshData):
    __metaclass__ = monkeypatch_class

    def fill_from_hedge(self, discr):
        self.discr = discr

        # add periodicity -----------------------------------------------------
        from pyrticle._internal import PeriodicityAxis
        for axis, ((ax_min, ax_max), periodicity_tags) in enumerate(zip(
                zip(*discr.mesh.bounding_box), discr.mesh.periodicity)):
            pa = PeriodicityAxis()
            if periodicity_tags is not None:
                pa.min = ax_min
                pa.max = ax_max
            else:
                pa.min = 0
                pa.max = 0
            self.periodicities.append(pa)

        # add elements --------------------------------------------------------
        mesh = discr.mesh

        neighbor_map = {}
        for face, (e2, f2) in discr.mesh.both_interfaces():
            neighbor_map[face] = e2.id
        from hedge.mesh import TAG_ALL
        for face in discr.mesh.tag_to_boundary[TAG_ALL]:
            neighbor_map[face] = MeshData.INVALID_ELEMENT

        from pyrticle._internal import ElementInfo, FaceInfo

        self.element_info.reserve(len(mesh.elements))
        for i, el in enumerate(mesh.elements):
            ei = ElementInfo()
            ei.id = i
            ei.inverse_map = el.inverse_map
            ei.jacobian = abs(el.map.jacobian)
            ei.start, ei.end = discr.find_el_range(el.id)
            ei.vertices.extend([vi for vi in el.vertex_indices])

            all_face_vertex_indices = el.face_vertices(el.vertex_indices)

            for face_idx, face_normal in enumerate(el.face_normals):
                fi = FaceInfo()
                fi.normal = face_normal
                fi.neighbor = neighbor_map[el, face_idx] 

                this_fvi = all_face_vertex_indices[face_idx]
                try:
                    fi.neighbor_periodicity_axis = mesh.periodic_opposite_faces[this_fvi][1]
                except KeyError:
                    fi.neighbor_periodicity_axis = MeshData.INVALID_AXIS

                this_face_vertices = [mesh.points[fvi] 
                        for fvi in all_face_vertex_indices[face_idx]]

                rhs = [numpy.dot(face_normal, fv) for fv in this_face_vertices]
                from pytools import all_roughly_equal, average
                assert all_roughly_equal(rhs, 1e-14)

                fi.face_plane_eqn_rhs = average(rhs)
                fi.face_centroid = average(this_face_vertices)

                fi.face_radius_from_centroid = max(
                        la.norm(fv - fi.face_centroid) for fv in this_face_vertices)

                ei.faces.append(fi)

            self.element_info.append(ei)

        # add vertices --------------------------------------------------------
        vertex_to_element_map = {} # map vertex to (via_periodic, el_id)

        for el in mesh.elements:
            for vi in el.vertex_indices:
                vertex_to_element_map.setdefault(vi, set()).add((
                    _internal.MeshData.INVALID_AXIS, el.id))
                for other_vi, per_axis in mesh.periodic_opposite_vertices.get(vi, []):
                    vertex_to_element_map.setdefault(other_vi, set()).add((per_axis, el.id))

        self.vertices.reserve(len(mesh.points))
        self.vertices.extend(mesh.points)

        from pyrticle._internal import UnsignedVector
        self.vertex_adj_elements.reserve(
                2*discr.dimensions*len(mesh.points))
        self.vertex_adj_element_starts.reserve(len(mesh.points))
        self.vertex_adj_element_starts.append(0)

        for vi in xrange(len(mesh.points)):
            adj_periodicity_axes, adj_ids = zip(*vertex_to_element_map[vi])
            self.vertex_adj_elements.extend(adj_ids)
            self.vertex_adj_periodicity_axes.extend(adj_periodicity_axes)
            self.vertex_adj_element_starts.append(
                    len(self.vertex_adj_elements))

        # add nodes -----------------------------------------------------------
        self.nodes.reserve(len(discr.nodes))
        self.nodes.extend(discr.nodes)

    def min_vertex_distance_for_el(self, el):
        vertices = [self.discr.mesh.points[vi] 
                for vi in el.vertex_indices]

        return min(min(la.norm(vi-vj)
                for i, vi in enumerate(vertices)
                if i != j)
                for j, vj in enumerate(vertices))

    def advisable_particle_radius(self):
        vertex_distances = [self.min_vertex_distance_for_el(el) 
                for el in self.discr.mesh.elements]
        vertex_distances.sort()
        return 0.6 * vertex_distances[int(0.25*len(vertex_distances))]

    def min_vertex_distance(self):
        return min(self.min_vertex_distance_for_el(el) 
                for el in self.discr.mesh.elements)
