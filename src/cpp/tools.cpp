// Pyrticle - Particle in Cell in Python
// Little bits of helpfulness
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





#include "tools.hpp"
#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>




pyrticle::warning_listener *pyrticle::warning_listener::m_singleton = 0;




static inline double sphere_area(double dims)
{
  return 2*pow(M_PI, dims/2) / tgamma(dims/2);
}




pyrticle::polynomial_shape_function::polynomial_shape_function( 
    double radius, unsigned dimensions, double alpha)
: m_alpha(alpha), m_radius(radius), 
  m_radius_squared(square(radius))
{
  using boost::math::tgamma;
  using boost::math::beta;

  const double n = dimensions;

  // see doc/notes.tm
  m_normalizer = 2 / (sphere_area(n) *
    pow(m_radius, n+alpha)*beta(n/2, alpha+1)
    );
}




pyrticle::c_infinity_shape_function::c_infinity_shape_function( 
    double radius, unsigned dimensions, double integral_for_rad1)
: m_radius(radius), m_radius_squared(square(radius))
{
  using boost::math::tgamma;

  const double n = dimensions;

  m_normalizer = 1/(sphere_area(n)*pow(radius, n)*integral_for_rad1);
}
