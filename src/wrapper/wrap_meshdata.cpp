// Pyrticle - Particle in Cell in Python
// Main Module of the Python wrapper
// Copyright (C) 2007 Andreas Kloeckner
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.




#include <pyublas/numpy.hpp>
#include <boost/shared_ptr.hpp>
#include "meshdata.hpp"
#include "wrap_helpers.hpp"




namespace python = boost::python;
using namespace pyrticle;



namespace
{
  py_vector get_normal(const mesh_data::face_info &fi)
  { return fi.m_normal; }
  void set_normal(mesh_data::face_info &fi, py_vector n)
  { fi.m_normal = n; }
  py_vector get_face_centroid(const mesh_data::face_info &fi)
  { return fi.m_face_centroid; }
  void set_face_centroid(mesh_data::face_info &fi, py_vector n)
  { fi.m_face_centroid = n; }

}




void expose_meshdata()
{
  {
    typedef mesh_data cl;
    python::class_<cl, boost::noncopyable>("MeshData", python::init<unsigned>())
      .add_static_property("INVALID_ELEMENT", &cl::get_INVALID_ELEMENT)
      .add_static_property("INVALID_AXIS", &cl::get_INVALID_AXIS)

      .DEF_RO_MEMBER(dimensions)

      .DEF_RO_MEMBER(element_info)
      .DEF_SIMPLE_METHOD(set_vertices)
      .DEF_SIMPLE_METHOD(set_nodes)

      .DEF_RO_MEMBER(vertex_adj_element_starts)
      .DEF_RO_MEMBER(vertex_adj_elements)
      .DEF_RO_MEMBER(vertex_adj_periodicity_axes)

      .DEF_RO_MEMBER(periodicities)

      .def("is_in_element", &cl::is_in_element<py_vector>)
      .def("find_containing_element", &cl::find_containing_element<py_vector>)
      ;
  }

  {
    typedef mesh_data::face_info cl;
    python::class_<cl>("FaceInfo")
      // PyUblas member-in-base-class issue: wrap by hand
      .add_property("normal", get_normal, set_normal)
      .DEF_RW_MEMBER(neighbor)
      .DEF_RW_MEMBER(neighbor_periodicity_axis)

      .DEF_RW_MEMBER(face_plane_eqn_rhs)

      .add_property("face_centroid", get_face_centroid, set_face_centroid)
      .DEF_RW_MEMBER(face_radius_from_centroid)
      ;
  }
  {
    typedef mesh_data::element_info cl;
    python::class_<cl>("ElementInfo")
      .DEF_RW_MEMBER(id)
      .DEF_RW_MEMBER(inverse_map)
      .DEF_RW_MEMBER(jacobian)
      .DEF_RW_MEMBER(start)
      .DEF_RW_MEMBER(end)
      .DEF_RW_MEMBER(vertices)

      .DEF_RW_MEMBER(faces)
      ;
  }

  {
    typedef mesh_data::periodicity_axis cl;
    python::class_<cl>("PeriodicityAxis")
      .DEF_RW_MEMBER(min)
      .DEF_RW_MEMBER(max)
      ;
  }

  expose_std_vector<mesh_data::face_info>("FaceInfo");
  expose_std_vector<mesh_data::element_info>("ElementInfo");
  expose_std_vector<mesh_data::periodicity_axis>("PeriodicityAxis");
}
