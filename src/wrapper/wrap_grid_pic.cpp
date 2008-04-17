// Pyrticle - Particle in Cell in Python
// Python wrapper for PIC algorithm
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




#include "wrap_pic.hpp"




using namespace pyrticle;
namespace python = boost::python;




namespace
{
  python::object brick_it_iter(python::object self)
  { return self; }

  bounded_int_vector brick_it_next(grid_reconstructor::brick::iterator &it)
  {
    if (it.at_end())
    {
      PyErr_SetNone(PyExc_StopIteration);
      throw python::error_already_set();
    }

    bounded_int_vector result = *it;
    ++it;
    return result;
  }
}




void expose_grid_pic()
{
  expose_pic_nontarget_pushers_all_dim<grid_reconstructor>();

  {
    typedef grid_reconstructor::brick::iterator cl;
    python::class_<cl>("BrickIterator", 
        python::init<
        grid_reconstructor::brick const &, bounded_int_box const &>(
          python::args("brick", "bounds"))[
        python::with_custodian_and_ward<1,2>()])
      .def("__iter__", brick_it_iter)
      .def("next", brick_it_next)
      ;
  }

  {
    typedef grid_reconstructor::brick cl;
    python::class_<cl>("Brick", 
        python::init<
        grid_reconstructor::brick_number, 
        grid_reconstructor::grid_node_number,
        bounded_vector, bounded_vector, bounded_int_vector>(
          python::args("number", "start_index", "stepwidths", "origin",
            "dimensions")))
      .add_property("number", &cl::number)
      .add_property("start_index", &cl::start_index)
      .add_property("stepwidths", 
          make_function(&cl::stepwidths, 
            python::return_value_policy<python::return_by_value>()))
      .add_property("origin", 
          make_function(&cl::origin, 
            python::return_value_policy<python::return_by_value>()))
      .add_property("dimensions", 
          make_function(&cl::dimensions, 
            python::return_value_policy<python::return_by_value>()))

      .def("__len__", &cl::node_count)
      .DEF_SIMPLE_METHOD(point)
      .DEF_SIMPLE_METHOD(index)
      .DEF_SIMPLE_METHOD(bounding_box)
      .DEF_SIMPLE_METHOD(index_range)
      ;
  }

  {
    typedef grid_reconstructor::element_on_grid cl;
    python::class_<cl>("ElementOnGrid")
      .DEF_RW_MEMBER(grid_nodes)
      .DEF_BYVAL_RW_MEMBER(interpolation_matrix)
      ;
  }

  expose_std_vector<grid_reconstructor::brick>("Brick");
  expose_std_vector<grid_reconstructor::element_on_grid>("ElementOnGrid");
}
