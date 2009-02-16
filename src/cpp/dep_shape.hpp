// Pyrticle - Particle in Cell in Python
// Deposition based on shape functions
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





#ifndef _BADFJAH_PYRTICLE_DEP_SHAPE_HPP_INCLUDED
#define _BADFJAH_PYRTICLE_DEP_SHAPE_HPP_INCLUDED




#include <boost/ref.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/typeof/std/utility.hpp>
#include "tools.hpp"
#include "bases.hpp"
#include "meshdata.hpp"
#include "particle_state.hpp"
#include "element_finder.hpp"




namespace pyrticle 
{
  template <class VecType>
  inline
  bool is_not_near_vertex(VecType const &unit_pt)
  {
    bool not_near_vertex = true;

    double uc_sum = 0;
    BOOST_FOREACH(double uc, unit_pt)
    {
      if (uc > 0)
      {
        not_near_vertex = false;
        break;
      }
      uc_sum += uc;
    }

    return not_near_vertex && (1-0.5*(uc_sum+unit_pt.size()) < 0.5);
  }




  template <class ParticleState, class ShapeFunction>
  struct shape_function_depositor
  {
    private:
      template <class Target>
      class element_target
      {
        private:
          const mesh_data   &m_mesh_data;
          ShapeFunction     m_shape_function;
          const double      m_charge;
          Target            &m_target;
          
        public:
          element_target(
              const mesh_data &md, 
              const ShapeFunction &sf,
              double charge, Target &tgt)
            : m_mesh_data(md), m_shape_function(sf), m_charge(charge), m_target(tgt)
          { }

          void add_shape_on_element(
              const bounded_vector &center,
              const mesh_data::element_number en
              ) const
          {
            const mesh_data::element_info &einfo(
                m_mesh_data.m_element_info[en]);

            const unsigned el_length = einfo.m_end-einfo.m_start;
            dyn_vector el_rho(el_length);

            for (unsigned i = 0; i < el_length; i++)
              el_rho[i] = 
                m_charge * m_shape_function(
                    m_mesh_data.mesh_node(einfo.m_start+i) 
                    - center);

            m_target.add_shape_on_element(en, einfo.m_start, el_rho);
          }
      };

    public:
      typedef ParticleState particle_state;

      struct depositor_state
      { };

      ShapeFunction m_shape_function;
      const mesh_data &m_mesh_data;




      shape_function_depositor(const mesh_data &md)
        : m_mesh_data(md)
      { }
    



      template<class Target>
      void deposit_densities_on_target(
          const depositor_state &ds,
          const particle_state &ps,
          Target &tgt,
          boost::python::slice const &pslice) const
      {
        element_finder el_finder(m_mesh_data);

        FOR_ALL_SLICE_INDICES(particle_number, pn, pslice, ps.particle_count)
        {
          element_target<Target> el_target(m_mesh_data, m_shape_function, ps.charges[pn], tgt);

          tgt.begin_particle(pn);
          el_finder(ps, el_target, pn, m_shape_function.radius());
          tgt.end_particle(pn);
        }
      }
  };
}




#endif
