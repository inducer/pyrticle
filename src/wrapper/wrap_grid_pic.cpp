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




void expose_grid_pic()
{
  expose_pic_nontarget_pushers_all_dim<grid_reconstructor>();

  {
    typedef grid_reconstructor::rec_brick cl;
    python::class_<cl, python::bases<brick> >("RecBrick", 
        python::init<
        grid_reconstructor::brick_number, 
        grid_node_number,
        bounded_vector, bounded_vector, bounded_int_vector>(
          python::args("number", "start_index", "stepwidths", "origin",
            "dimensions")))
      .add_property("number", &cl::number)
      ;
  }

  expose_std_vector<grid_reconstructor::rec_brick>("RecBrick");
  expose_std_vector<grid_reconstructor::element_on_grid>("ElementOnGrid");

  {
    typedef grid_reconstructor::element_on_grid cl;
    python::class_<cl>("ElementOnGrid")
      .DEF_RW_MEMBER(element_number)
      .DEF_RW_MEMBER(grid_nodes)
      .DEF_BYVAL_RW_MEMBER(weight_factors)
      .DEF_BYVAL_RW_MEMBER(interpolation_matrix)
      .DEF_BYVAL_RW_MEMBER(inverse_interpolation_matrix)
      ;
  }
}
