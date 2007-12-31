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




namespace python = boost::python;
using namespace pyrticle;




void expose_pusher()
{
  {
    typedef monomial_basis_function cl;
    python::class_<cl>("MonomialBasisFunction", python::init<unsigned, unsigned>())
      .def(python::init<unsigned, unsigned, unsigned>())
      .def(python::init<const std::vector<unsigned> & >())
      .def("__call__", &cl::operator())
      ;
  }

  expose_std_vector<monomial_basis_function>(
      "MonomialBasisFunction");

  {
    typedef local_monomial_discretization cl;
    python::class_<cl>("LocalMonomialDiscretization")
      .DEF_RW_MEMBER(basis)
      .DEF_RW_MEMBER(l_vandermonde_t)
      .DEF_RW_MEMBER(u_vandermonde_t)
      .DEF_RW_MEMBER(p_vandermonde_t)
      ;
  }

  expose_std_vector<local_monomial_discretization>(
      "LocalMonomialDiscretization");
}

