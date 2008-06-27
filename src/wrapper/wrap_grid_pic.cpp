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




void expose_grid_pic()
{
  expose_pic_nontarget_pushers_all_dim<grid_reconstructor<jiggly_brick> >();

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

  {
    typedef element_on_grid cl;
    python::class_<cl>("ElementOnGrid")
      .DEF_RW_MEMBER(element_number)
      .DEF_RW_MEMBER(grid_nodes)
      .DEF_BYVAL_RW_MEMBER(weight_factors)
      .DEF_BYVAL_RW_MEMBER(interpolation_matrix)
      .DEF_BYVAL_RW_MEMBER(inverse_interpolation_matrix)
      ;
  }
  expose_std_vector<element_on_grid>("ElementOnGrid");
}
