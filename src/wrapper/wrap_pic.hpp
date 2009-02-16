// Pyrticle - Particle in Cell in Python
// Python wrapper for PIC algorithm
// Copyright (C) 2007 Andreas Kloeckner
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or // (at your option) any later version.  // 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.




#ifndef _AFAVCZ_PYRTICLE_WRAP_PIC_HPP_INCLUDED
#define _AFAVCZ_PYRTICLE_WRAP_PIC_HPP_INCLUDED




#include <boost/lexical_cast.hpp>
#include "particle_state.hpp"
#include "diagnostics.hpp"
#include "wrap_helpers.hpp"




namespace python = boost::python;
using namespace pyrticle;




namespace
{
  template <class ParticleState>
  std::string get_state_class_suffix()
  {
    std::string result;
    result += boost::lexical_cast<std::string>(ParticleState::xdim());
    result += boost::lexical_cast<std::string>(ParticleState::vdim());
    return result;
  }

  


#define EXPOSE_FOR_ALL_STATE_TYPES(NAME, ARGLIST) \
  NAME<particle_base_state<2, 2> >ARGLIST; \
  NAME<particle_base_state<3, 3> >ARGLIST;

#define EXPOSE_FOR_ALL_TARGET_RECONSTRUCTORS(NAME, SHAPEFUNC, ARGLIST) \
  NAME<shape_function_depositor<ParticleState, SHAPEFUNC > > ARGLIST; \
  NAME<normalized_shape_function_depositor<ParticleState, SHAPEFUNC > > ARGLIST; \
  NAME<advective_depositor<ParticleState, SHAPEFUNC > > ARGLIST;
}




#endif
