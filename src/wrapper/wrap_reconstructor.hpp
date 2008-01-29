// Pyrticle - Particle in Cell in Python
// Python wrapper for reconstruction bits
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





#ifndef _AYYTYAH_PYRTICLE_WRAP_RECONSTRUCTOR_HPP_INCLUDED
#define _AYYTYAH_PYRTICLE_WRAP_RECONSTRUCTOR_HPP_INCLUDED




#include "wrap_helpers.hpp"
#include "rec_shape.hpp"
#include "rec_advective.hpp"




namespace pyrticle
{
  template <class Wrapper, class PIC>
  void expose_typed_reconstructor(Wrapper &wrp, 
      shape_function_reconstructor::type<PIC> *)
  { 
    typedef shape_function_reconstructor::type<PIC> cl;
    wrp
      .DEF_SIMPLE_METHOD(set_radius)
      ;
  }




  template <class Wrapper, class PIC>
  void expose_typed_reconstructor(Wrapper &wrp, 
      advective_reconstructor::type<PIC> *)
  { 
    using boost::python::arg;

    typedef advective_reconstructor::type<PIC> cl;
    wrp
      .def("setup_advective_reconstructor", 
          &cl::setup_advective_reconstructor,
          (arg("dofs_per_element")))
      .DEF_RW_MEMBER(rho_dof_shift_listener)
      .DEF_SIMPLE_METHOD(add_local_diff_matrix)
      .DEF_SIMPLE_METHOD(add_advective_particle)
      .DEF_SIMPLE_METHOD(get_debug_quantity_on_mesh)
      .DEF_SIMPLE_METHOD(get_advective_particle_rhs)
      .DEF_SIMPLE_METHOD(apply_advective_particle_rhs)
      ;
  }
}




#endif
