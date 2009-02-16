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




#include "grid.hpp"
#include "wrap_helpers.hpp"




using namespace pyrticle;
namespace python = boost::python;




namespace
{
  python::object brick_it_iter(python::object self)
  { return self; }

  template <class Brick>
  bounded_int_vector brick_it_next(
      typename Brick::iterator &it)
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




void expose_grid()
{
  {
    typedef brick::iterator cl;
    python::class_<cl>("BrickIterator", 
        python::init<
        brick const &, bounded_int_box const &>(
          python::args("brick", "bounds"))[
        python::with_custodian_and_ward<1,2>()])
      .def("__iter__", brick_it_iter)
      .def("next", brick_it_next<brick>)
      ;
  }

  {
    typedef brick cl;
    python::class_<cl>("Brick", 
        python::init<brick_number,
        grid_node_number, bounded_vector, bounded_vector, 
        bounded_int_vector>(
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
      .DEF_SIMPLE_METHOD(split_index)
      .DEF_SIMPLE_METHOD(which_cell)
      .DEF_SIMPLE_METHOD(bounding_box)
      .DEF_SIMPLE_METHOD(cell)
      .DEF_SIMPLE_METHOD(index_range)

      .def("get_iterator",
          (cl::iterator (cl::*)(bounded_box const &) const) &cl::get_iterator)
      .def("get_iterator",
          (cl::iterator (cl::*)(bounded_int_box const &) const) &cl::get_iterator)
      ;
  }
  expose_std_vector<brick>("Brick");

  {
    typedef jiggly_brick::iterator cl;
    python::class_<cl>("JigglyBrickIterator", 
        python::init<
        jiggly_brick const &, 
        bounded_int_box const &>(
          python::args("brick", "bounds"))[
        python::with_custodian_and_ward<1,2>()])
      .def("__iter__", brick_it_iter)
      .def("next", brick_it_next<jiggly_brick>)
      ;
  }

  {
    typedef jiggly_brick cl;
    python::class_<cl, python::bases<brick> >("JigglyBrick", 
        python::init<
        brick_number, grid_node_number,
        bounded_vector, bounded_vector, bounded_int_vector, 
        double>(
          python::args("number", "start_index", "stepwidths", "origin",
            "dimensions", "jiggle_radius")
          )
        )
      ;
  }
  expose_std_vector<jiggly_brick>("JigglyBrick");
}

