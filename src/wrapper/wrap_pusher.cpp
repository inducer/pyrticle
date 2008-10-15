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
#include "wrap_helpers.hpp"
#include "wrap_pic.hpp"




using namespace  boost::python;
using namespace pyrticle;



namespace
{
  template <class Pusher, class Wrapper>
  void expose_force_calculator(Wrapper &wrp)
  {
    typedef Pusher cl;

    if (Pusher::particle_state::vdim() == 3)
    {
      wrp
        .def("forces", &cl::template forces< // full-field case
            py_vector, py_vector, py_vector,
            py_vector, py_vector, py_vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("bx"), arg("by"), arg("bz"),
             arg("velocities"), arg("verbose_vis")))
        ;
    }
    else if (Pusher::particle_state::vdim() == 2)
    {
      wrp
        /*
        .def("forces", &cl::template forces< // TM case
            zero_vector, zero_vector, py_vector,
            py_vector, py_vector, zero_vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("bx"), arg("by"), arg("bz"),
             arg("velocities"), arg("verbose_vis")))
             */
        .def("forces", &cl::template forces< // TE case
            py_vector, py_vector, zero_vector,
            zero_vector, zero_vector, py_vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("bx"), arg("by"), arg("bz"),
             arg("velocities"), arg("verbose_vis")))
        ;
    }
  }

  template <class ParticleState>
  void expose_pushers_for_pstate()
  {
    {
      typedef monomial_particle_pusher<ParticleState> cl;
      class_<cl> wrp("MonomialPusher", init<const mesh_data &>())
        ;
      expose_force_calculator<cl>(wrp);
    }
    
  }
}




void expose_pusher()
{
  {
    typedef monomial_basis_function cl;
    class_<cl>("MonomialBasisFunction", init<unsigned, unsigned>())
      .def(init<unsigned, unsigned, unsigned>())
      .def(init<const std::vector<unsigned> & >())
      .def("__call__", 
          (const double (cl::*)(const py_vector &) const)
          &cl::operator())
      ;
  }

  expose_std_vector<monomial_basis_function>(
      "MonomialBasisFunction");

  {
    typedef local_monomial_discretization cl;
    class_<cl>("LocalMonomialDiscretization")
      .DEF_RW_MEMBER(basis)
      .DEF_BYVAL_RW_MEMBER(lu_vandermonde_t)
      .DEF_BYVAL_RW_MEMBER(lu_piv_vandermonde_t)
      ;
  }

  expose_std_vector<local_monomial_discretization>(
      "LocalMonomialDiscretization");

  EXPOSE_FOR_ALL_STATE_TYPES(expose_pushers_for_pstate, ());
}

