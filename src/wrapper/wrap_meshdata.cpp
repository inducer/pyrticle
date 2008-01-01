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




#include <boost/shared_ptr.hpp>
#include "meshdata.hpp"
#include "wrap_helpers.hpp"




namespace python = boost::python;
using namespace pyrticle;




void expose_meshdata()
{
  {
    typedef mesh_data cl;
    python::class_<cl, boost::noncopyable>("MeshData", python::init<unsigned>())
      .add_static_property("INVALID_ELEMENT", &cl::get_INVALID_ELEMENT)

      .DEF_RW_MEMBER(nodes)
      .DEF_RO_MEMBER(dimensions)

      .DEF_RO_MEMBER(element_info)
      .DEF_RO_MEMBER(vertices)
      .DEF_RO_MEMBER(nodes)

      .DEF_RO_MEMBER(vertex_adj_element_starts)
      .DEF_RO_MEMBER(vertex_adj_elements)

      .DEF_RO_MEMBER(periodicities)

      .DEF_SIMPLE_METHOD(is_in_element)
      .DEF_SIMPLE_METHOD(find_containing_element)
      ;
  }

  {
    typedef mesh_data::element_info cl;
    python::class_<cl>("ElementInfo")
      .DEF_RW_MEMBER(id)
      .DEF_RW_MEMBER(inverse_map)
      .DEF_RW_MEMBER(start)
      .DEF_RW_MEMBER(end)
      .DEF_RW_MEMBER(vertices)
      .DEF_RW_MEMBER(normals)
      .DEF_RW_MEMBER(neighbors)
      ;
  }

  {
    typedef mesh_data::periodicity_axis cl;
    python::class_<cl>("PeriodicityAxis")
      .DEF_RW_MEMBER(axis)
      .DEF_RW_MEMBER(min)
      .DEF_RW_MEMBER(max)
      ;
  }

  expose_std_vector<mesh_data::element_info>("ElementInfo");
  expose_std_vector<mesh_data::periodicity_axis>("PeriodicityAxis");
  expose_std_vector<hedge::vector>("Vector");
}
