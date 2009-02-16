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
#include "particle_state.hpp"




namespace pyrticle
{
  template <class ParticleState>
  const py_vector kinetic_energies(ParticleState const &ps, double vacuum_c)
  {
    const unsigned vdim = ps.vdim();

    py_vector result(ps.particle_count);

    const double c_squared = vacuum_c*vacuum_c;

    for (particle_number pn = 0; pn < ps.particle_count; pn++)
    {
      const unsigned vpstart = vdim*pn;
      const unsigned vpend = vdim*(pn+1);

      const double m = ps.masses[pn];
      double p = norm_2(subrange(ps.momenta, vpstart, vpend));
      double v = vacuum_c*p/sqrt(m*m*c_squared + p*p);
      result[pn] = (p/v-m)*c_squared;
    }
    return result;
  }



  template <class ParticleState>
  const double total_charge(ParticleState const &ps)
  {
    double result = 0;
    for (particle_number pn = 0; pn < ps.particle_count; pn++)
      result += ps.charges[pn];
    return result;
  }




  template <class ParticleState>
  const py_vector particle_momentum(ParticleState const &ps)
  {
    const unsigned vdim = ps.vdim();
    py_vector result(vdim);
    result.clear();

    for (particle_number pn = 0; pn < ps.particle_count; pn++)
      result += subrange(ps.momenta, vdim*pn, vdim*(pn+1));

    return result;
  }





  template <class ParticleState>
  const double rms_beam_size(ParticleState const &ps, unsigned axis)
  {
    if (ps.particle_count == 0)
      return 0;

    double result = 0;
    for (particle_number pn = 0; pn < ps.particle_count; pn++)
      result += square(ps.positions[pn*ps.xdim() + axis]);

    return sqrt(result/ps.particle_count);
  }




  template <class ParticleState>
  const double rms_beam_emittance(ParticleState const &ps, unsigned axis, unsigned beam_axis)
  {
    if (ps.particle_count == 0)
      return 0;

    double mean_x = 0;
    double mean_xp = 0;
    double mean_x_squared = 0;
    double mean_xp_squared = 0;
    double mean_xxp = 0;

    // see doc/notes.tm
    for (particle_number pn = 0; pn < ps.particle_count; pn++)
    {
      const double x = ps.positions[pn*ps.xdim() + axis];

      mean_x += x;
      mean_x_squared += square(x);

      const double px = ps.momenta[pn*ps.vdim() + axis];
      const double pz = ps.momenta[pn*ps.vdim() + beam_axis];

      const double xprime = pz ? px/pz : 0;

      mean_xp += xprime;
      mean_xp_squared += square(xprime);

      mean_xxp += x*xprime;
    }

    mean_x /= ps.particle_count;
    mean_xp /= ps.particle_count;

    mean_x_squared /= ps.particle_count;
    mean_xp_squared /= ps.particle_count;

    mean_xxp /= ps.particle_count;

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




  template <class ParticleState>
  const double rms_energy_spread(ParticleState const &ps, double vacuum_c)
  {
    if (ps.particle_count == 0)
      return 0;
    py_vector energies = kinetic_energies(ps, vacuum_c);
    return std_dev(energies.begin(), energies.end());
  }




  template <class ParticleState>
  const py_vector particle_current(ParticleState const &ps, py_vector const &velocities,
      double length)
  {
    const unsigned vdim = ps.vdim();
    py_vector result(vdim);
    result.clear();

    for (particle_number pn = 0; pn < ps.particle_count; pn++)
      result += ps.charges[pn] * subrange(velocities, vdim*pn, vdim*(pn+1));

    return result / length;
  }
}




#endif
