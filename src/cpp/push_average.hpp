// Pyrticle - Particle in Cell in Python
// Particle pusher based on weighted averaging
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




#ifndef _AFYAFDFYA_PYRTICLE_PUSH_AVERAGE_HPP_INCLUDED
#define _AFYAFDFYA_PYRTICLE_PUSH_AVERAGE_HPP_INCLUDED




#include "meshdata.hpp"
#include "bases.hpp"
namespace pyrticle
{
  struct averaging_particle_pusher
  {
    template <class PICAlgorithm>
    class type : public pusher_base
    {
      public:
        static const char *get_name()
        { return "Average"; }

        // why all these template arguments? In 2D and 1D,
        // instead of passing a hedge::vector, you may simply
        // pass a zero_vector, and interpolation will know to
        // not even compute anything, but just return zero.
        template <class EX, class EY, class EZ, 
                 class BX, class BY, class BZ>
        hedge::vector forces(
            const EX &ex, const EY &ey, const EZ &ez,
            const BX &bx, const BY &by, const BZ &bz,
            const hedge::vector &velocities,
            bool verbose_vis
            )
        {
        }
}




#endif
