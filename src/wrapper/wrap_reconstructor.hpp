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
#include "rec_normshape.hpp"
#include "rec_advective.hpp"
#include "rec_grid.hpp"




namespace pyrticle
{
  template <class Wrapper, class PIC>
  void expose_typed_reconstructor(Wrapper &wrp, 
      shape_function_reconstructor::type<PIC> *)
  { 
    typedef shape_function_reconstructor::type<PIC> cl;
    wrp
      .DEF_RW_MEMBER(shape_function)
      ;
  }




  template <class Wrapper, class PIC>
  void expose_typed_reconstructor(Wrapper &wrp, 
      normalized_shape_function_reconstructor::type<PIC> *)
  { 
    typedef normalized_shape_function_reconstructor::type<PIC> cl;
    wrp
      .DEF_RO_MEMBER(normalization_stats)
      .DEF_RO_MEMBER(centroid_distance_stats)
      .DEF_RO_MEMBER(el_per_particle_stats)
      .DEF_RW_MEMBER(shape_function)
      .DEF_SIMPLE_METHOD(setup_normalized_shape_reconstructor)
      ;
  }




  template <class Wrapper, class PIC>
  void expose_typed_reconstructor(Wrapper &wrp, 
      typename advective_reconstructor::type<PIC> *)
  { 
    using boost::python::arg;

    typedef typename advective_reconstructor::template type<PIC> cl;
    wrp
      .DEF_SIMPLE_METHOD(setup_advective_reconstructor)
      .DEF_RW_MEMBER(rho_dof_shift_listener)

      .DEF_RO_MEMBER(active_elements)

      .DEF_RO_MEMBER(element_activation_counter)
      .DEF_RO_MEMBER(element_kill_counter)

      .DEF_SIMPLE_METHOD(add_local_diff_matrix)
      .DEF_SIMPLE_METHOD(count_advective_particles)
      .DEF_SIMPLE_METHOD(add_advective_particle)
      .DEF_SIMPLE_METHOD(clear_advective_particles)
      .DEF_SIMPLE_METHOD(get_debug_quantity_on_mesh)
      .DEF_SIMPLE_METHOD(get_advective_particle_rhs)
      .DEF_SIMPLE_METHOD(apply_advective_particle_rhs)
      ;
  }




  template <class Wrapper, class PIC>
  void expose_typed_reconstructor(Wrapper &wrp, 
      typename grid_reconstructor::type<PIC> *)
  { 
    typedef grid_reconstructor::type<PIC> cl;
    wrp
      .DEF_RW_MEMBER(shape_function)
      .DEF_RW_MEMBER(bricks)

      .DEF_SIMPLE_METHOD(add_brick)
      .DEF_SIMPLE_METHOD(commit_bricks)

      .DEF_SIMPLE_METHOD(get_grid_rho)
      .DEF_SIMPLE_METHOD(get_grid_j)
      ;
  }
}




#endif
