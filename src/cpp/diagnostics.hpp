// Pyrticle - Particle in Cell in Python
// Diagnostics
// Copyright (C) 2007, 2008 Andreas Kloeckner
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




#ifndef _AFHYYHFA_PYRTICLE_DIAGNOSTICS_HPP_INCLUDED
#define _AFHYYHFA_PYRTICLE_DIAGNOSTICS_HPP_INCLUDED




#include <cmath>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include "pic_algorithm.hpp"




namespace pyrticle
{
  template <class PIC>
  const py_vector kinetic_energies(PIC const &pic)
  {
    const unsigned vdim = pic.get_dimensions_velocity();

    py_vector result(pic.m_particle_count);

    const double c_squared = pic.m_vacuum_c*pic.m_vacuum_c;

    for (particle_number pn = 0; pn < pic.m_particle_count; pn++)
    {
      const unsigned vpstart = vdim*pn;
      const unsigned vpend = vdim*(pn+1);

      const double m = pic.m_masses[pn];
      double p = norm_2(subrange(pic.m_momenta, vpstart, vpend));
      double v = pic.m_vacuum_c*p/sqrt(m*m*c_squared + p*p);
      result[pn] = (p/v-m)*c_squared;
    }
    return result;
  }



  template <class PIC>
  const double total_charge(PIC const &pic)
  {
    double result = 0;
    for (particle_number pn = 0; pn < pic.m_particle_count; pn++)
      result += pic.m_charges[pn];
    return result;
  }




  template <class PIC>
  const py_vector particle_momentum(PIC const &pic)
  {
    const unsigned vdim = pic.get_dimensions_velocity();
    py_vector result(vdim);
    result.clear();

    for (particle_number pn = 0; pn < pic.m_particle_count; pn++)
      result += subrange(pic.m_momenta, vdim*pn, vdim*(pn+1));

    return result;
  }





  template <class PIC>
  const double rms_beam_size(PIC const &pic, unsigned axis)
  {
    if (pic.m_particle_count == 0)
      return 0;

    double result = 0;
    for (particle_number pn = 0; pn < pic.m_particle_count; pn++)
      result += square(pic.m_positions[pn*pic.get_dimensions_pos() + axis]);

    return sqrt(result/pic.m_particle_count);
  }




  template <class PIC>
  const double rms_beam_emittance(PIC const &pic, unsigned axis, unsigned beam_axis)
  {
    if (pic.m_particle_count == 0)
      return 0;

    double mean_x = 0;
    double mean_xp = 0;
    double mean_x_squared = 0;
    double mean_xp_squared = 0;
    double mean_xxp = 0;

    // see doc/notes-2.tm
    for (particle_number pn = 0; pn < pic.m_particle_count; pn++)
    {
      const double x = pic.m_positions[pn*pic.get_dimensions_pos() + axis];

      mean_x += x;
      mean_x_squared += square(x);

      const double px = pic.m_momenta[pn*pic.get_dimensions_velocity() + axis];
      const double pz = pic.m_momenta[pn*pic.get_dimensions_velocity() + beam_axis];

      const double xprime = pz ? px/pz : 0;

      mean_xp += xprime;
      mean_xp_squared += square(xprime);

      mean_xxp += x*xprime;
    }

    mean_x /= pic.m_particle_count;
    mean_xp /= pic.m_particle_count;

    mean_x_squared /= pic.m_particle_count;
    mean_xp_squared /= pic.m_particle_count;

    mean_xxp /= pic.m_particle_count;

    return sqrt(
        mean_x_squared*mean_xp_squared
        - square(mean_x)*mean_xp_squared
        - mean_x_squared*square(mean_xp)
        - (
          square(mean_xxp)
          -2*mean_xxp*mean_x*mean_xp
          )
        );
  }




  template <class PIC>
  const double rms_energy_spread(PIC const &pic)
  {
    if (pic.m_particle_count == 0)
      return 0;
    py_vector energies = kinetic_energies(pic);
    return std_dev(energies.begin(), energies.end());
  }




  template <class PIC>
  const py_vector particle_current(PIC const &pic, py_vector const &velocities,
      double length)
  {
    const unsigned vdim = pic.get_dimensions_velocity();
    py_vector result(vdim);
    result.clear();

    for (particle_number pn = 0; pn < pic.m_particle_count; pn++)
      result += pic.m_charges[pn] * subrange(velocities, vdim*pn, vdim*(pn+1));

    return result / length;
  }
}




#endif
