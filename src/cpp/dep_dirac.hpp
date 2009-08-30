// Pyrticle - Particle in Cell in Python
// Deposition based on shape functions
// Copyright (C) 2009 Andreas Kloeckner
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




#ifndef _KAYTHADVOD_PYRTICLE_DEP_DIRAC_HPP_INCLUDED
#define _KAYTHADVOD_PYRTICLE_DEP_DIRAC_HPP_INCLUDED




namespace pyrticle 
{
  template <class ParticleState>
  struct dirac_depositor
  {
    public:
      typedef ParticleState particle_state;

      struct depositor_state
      { };

      const mesh_data &m_mesh_data;




      dirac_depositor(const mesh_data &md)
        : m_mesh_data(md)
      { }
    



      template<class Target>
      void deposit_densities_on_target(
          const depositor_state &ds,
          const particle_state &ps,
          Target &tgt,
          boost::python::slice const &pslice) const
      {
        // element_finder el_finder(m_mesh_data);

        FOR_ALL_SLICE_INDICES_PREP(pslice, ps.particle_count)

        FOR_ALL_SLICE_INDICES_LOOP
        {
          FOR_ALL_SLICE_INDICES_INNER(particle_number, pn);

          element_target<Target> el_target(m_mesh_data, m_shape_function, ps.charges[pn], tgt);

          tgt.begin_particle(pn);
          el_finder(ps, el_target, pn, m_shape_function.radius());
          tgt.end_particle(pn);
        }
      }
  };
}




#endif
