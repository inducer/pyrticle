// Pyrticle - Particle in Cell in Python
// Stuff we need to store about the field solver's mesh
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




#include <numeric>
#include <boost/foreach.hpp>
#include "meshdata.hpp"




using namespace pyrticle;




const bool is_in_unit_simplex(const hedge::vector &unit_coords)
{
  const double eps = 1e-10;

  BOOST_FOREACH(hedge::vector::value_type ri, unit_coords)
    if (ri < -1-eps)
      return false;

  return std::accumulate(unit_coords.begin(), unit_coords.end(), 
      (double) 0) <= -(signed(unit_coords.size())-2)+eps;
}




const mesh_data::element_number mesh_data::find_containing_element(
    const hedge::vector &pt) const
{
  BOOST_FOREACH(const element_info &el, m_element_info)
    if (is_in_unit_simplex(el.m_inverse_map(pt)))
      return el.m_id;
  return INVALID_ELEMENT;
}

