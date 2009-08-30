// Pyrticle - Particle in Cell in Python
// Python wrapper for particle pusher bits
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





#include "push_monomial.hpp"
#include "push_average.hpp"
#include "wrap_helpers.hpp"
#include "wrap_pic.hpp"
#include "dep_shape.hpp"
#include "dep_normshape.hpp"
#include "dep_advective.hpp"




using namespace  boost::python;
using namespace pyrticle;



namespace
{
  template <class Depositor, class Wrapper>
  void expose_averaging_force_calculator(Wrapper &wrp)
  {
    typedef typename Wrapper::wrapped_type cl;
    typedef typename Depositor::particle_state particle_state;

    if (particle_state::m_vdim == 3)
    {
      wrp
        .def("forces", &cl::template forces< // full-field case
            py_vector, py_vector, py_vector,
            py_vector, py_vector, py_vector, Depositor>,
            args("ex","ey", "ez", "bx", "by", "bz", 
              "particle_state", "pusher_state", 
              "depositor", "depositor_state",
              "velocities", "vis_listener"))
        ;
    }
    else if (particle_state::m_vdim == 2)
    {
      wrp
        /*
        .def("forces", &cl::template forces< // TM case
            zero_vector, zero_vector, py_vector,
            py_vector, py_vector, zero_vector, Depositor>,
            args("ex","ey", "ez", "bx", "by", "bz",
              "particle_state", "pusher_state", 
              "depositor", "depositor_state",
              "velocities", "vis_listener"))
             */
        .def("forces", &cl::template forces< // TE case
            py_vector, py_vector, zero_vector,
            zero_vector, zero_vector, py_vector, Depositor>,
            args("ex","ey", "ez", "bx", "by", "bz", 
              "particle_state", "pusher_state", 
              "depositor", "depositor_state",
              "velocities", "vis_listener"))
        ;
    }
  }




  template <class ParticleState>
  void expose_pushers_for_pstate()
  {
    {
      typedef monomial_particle_pusher<ParticleState> cl;
      class_<cl> wrp(
          ("MonomialPusher"+get_state_class_suffix<ParticleState>()).c_str(), 
          init<const mesh_data &>());

      wrp
        .DEF_RW_MEMBER(local_discretizations)
        .DEF_RW_MEMBER(ldis_indices)
        ;

      if (ParticleState::m_vdim == 3)
      {
        wrp
          .def("forces", &cl::template forces< // full-field case
              py_vector, py_vector, py_vector,
              py_vector, py_vector, py_vector>,
              args("ex","ey", "ez", "bx", "by", "bz", "ps", "velocities", "vis_listener"))
          ;
      }
      else if (ParticleState::m_vdim == 2)
      {
        wrp
          /*
          .def("forces", &cl::template forces< // TM case
              zero_vector, zero_vector, py_vector,
              py_vector, py_vector, zero_vector>,
              args("ex","ey", "ez", "bx", "by", "bz", "ps", "velocities", "vis_listener"))
               */
          .def("forces", &cl::template forces< // TE case
              py_vector, py_vector, zero_vector,
              zero_vector, zero_vector, py_vector>,
              args("ex","ey", "ez", "bx", "by", "bz", "ps", "velocities", "vis_listener"))
          ;
      }
    }
    
    typedef averaging_particle_pusher<ParticleState> cl;
    class_<cl> wrp(
        ("AveragingPusher"+get_state_class_suffix<ParticleState>()).c_str(), 
        init<const mesh_data &, const py_matrix &>());

    EXPOSE_FOR_ALL_TARGET_RECONSTRUCTORS(
        expose_averaging_force_calculator,
        polynomial_shape_function,
        (wrp));

    scope cls_scope = wrp;
    {
      typedef typename cl::pusher_state cl;
      class_<cl>("PusherState")
        .DEF_RW_MEMBER(e_normalization_stats)
        .DEF_RW_MEMBER(b_normalization_stats)
        ;
    }
  }
}




void expose_pusher()
{
  EXPOSE_FOR_ALL_STATE_TYPES(expose_pushers_for_pstate, ());
}

