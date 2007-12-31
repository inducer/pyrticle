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




#include "rec_shape.hpp"




namespace pyrticle
{
  template <class Wrapper, class PIC>
  void expose_typed_reconstructor(
      Wrapper &wrp, 
      shape_function_reconstructor::type<PIC> *)
  { }
}




#endif
