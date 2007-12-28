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




#include "meshdata.hpp"
#include "wrap_helpers.hpp"




namespace python = boost::python;
using namespace pyrticle;




void expose_meshdata()
{
  {
    typedef mesh_data cl;
    python::class_<cl, boost::noncopyable>("MeshData", 
        python::init<unsigned>())
      .add_static_property("INVALID_ELEMENT", &cl::get_INVALID_ELEMENT)

      .DEF_RW_MEMBER(nodes)
      .DEF_RO_MEMBER(dimensions)

      //.DEF_SIMPLE_METHOD(add_local_discretization)
      //.DEF_SIMPLE_METHOD(add_element)
      //.DEF_SIMPLE_METHOD(add_vertex)
      //.DEF_SIMPLE_METHOD(add_nodes)
      .DEF_SIMPLE_METHOD(add_periodicity)

      .DEF_SIMPLE_METHOD(is_in_element)
      .DEF_SIMPLE_METHOD(find_containing_element)
      ;
  }

}
