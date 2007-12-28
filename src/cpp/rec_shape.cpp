// Pyrticle - Particle in Cell in Python
// Reconstruction based on shape functions
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




#include "rec_shape.hpp"
#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>





pyrticle::shape_function::shape_function( 
    double radius, unsigned dimensions, double alpha)
: m_alpha(alpha), m_l(radius), 
  m_l_squared(square(radius))
{
  using boost::math::tgamma;
  using boost::math::beta;

  double n = dimensions;

  // see doc/notes.tm
  double sphere_area = 2*pow(M_PI, n/2) / tgamma(n/2);
  m_normalizer = 2 / (sphere_area *
    pow(m_l, n+alpha)*beta(n/2, alpha+1)
    );
}
