// Pyrticle - Particle in Cell in Python
// Deposition based on normalized advected shapes
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




#ifndef _AJFYAA_PYRTICLE_DEP_NORMSHAPE_HPP_INCLUDED
#define _AJFYAA_PYRTICLE_DEP_NORMSHAPE_HPP_INCLUDED




#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/typeof/std/utility.hpp>
#include "tools.hpp"
#include "bases.hpp"
#include "meshdata.hpp"
#include "dep_target.hpp"
#include "dep_shape.hpp"




namespace pyrticle 
{
  struct normalized_deposition_stats
  {
    stats_gatherer<double>  m_normalization_stats;
    stats_gatherer<double>  m_centroid_distance_stats;
    stats_gatherer<double>  m_el_per_particle_stats;
  };




  template <class ParticleState, class ShapeFunction>
  struct normalized_shape_function_depositor
  {
    public:
      typedef ParticleState particle_state;

      struct depositor_state
      { 
        normalized_deposition_stats m_stats;
      };

    private:
      template <class Target>
      class normalizing_element_target
      {
        public:
          const normalized_shape_function_depositor &m_dep;
          const particle_state &m_pstate;
          depositor_state &m_dstate;
          Target m_target;
          dyn_vector m_shape_interpolant;
          unsigned m_used_shape_dofs;

          struct shape_element {
            mesh_data::element_number m_el_number;
            unsigned m_el_length;
            unsigned m_my_start_index;
            mesh_data::node_number m_global_start_index;

            shape_element(
                mesh_data::element_number el_number,
                unsigned el_length,
                unsigned my_start_index, 
                mesh_data::node_number global_start_index)
              : 
                m_el_number(el_number),
                m_el_length(el_length),
                m_my_start_index(my_start_index), 
                m_global_start_index(global_start_index)
            { }
          };

          std::vector<shape_element>        m_particle_shape_elements;
          double                            m_integral;

          normalizing_element_target(
              const normalized_shape_function_depositor &dep,
              depositor_state &dstate,
              const particle_state &pstate,
              const dyn_vector &integral_weights,
              Target &target
              )
            : 
              m_dep(dep),
              m_pstate(pstate),
              m_dstate(dstate),
              m_target(target), 
              m_shape_interpolant(100)
          { }




          void begin_particle()
          {
            m_used_shape_dofs = 0;
            m_particle_shape_elements.clear();
            m_integral = 0;
          }




          void add_shape_on_element(
              const bounded_vector &center,
              const mesh_data::element_number en
              )
          {
            const mesh_data::element_info &einfo = m_dep.m_mesh_data.m_element_info[en];
            unsigned element_length = einfo.m_end-einfo.m_start;

            {
              bounded_vector centroid = m_dep.m_mesh_data.element_centroid(en);
              m_dstate.m_stats.m_centroid_distance_stats.add(norm_2(centroid-center));
            }

            // make sure we have enough interpolant dofs available
            while (m_used_shape_dofs + element_length > m_shape_interpolant.size())
            {
              dyn_vector new_shape_interpolant(2*m_shape_interpolant.size());
              subrange(new_shape_interpolant, 0, m_shape_interpolant.size()) =
                m_shape_interpolant;
              new_shape_interpolant.swap(m_shape_interpolant);
            }

            shape_element new_shape_element(
                einfo.m_id,
                element_length,
                m_used_shape_dofs, 
                einfo.m_start);
            m_used_shape_dofs += element_length;

            double el_integral = 0;
            for (unsigned i = 0; i < element_length; i++)
            {
              double shapeval = m_dep.m_shape_function(
                  m_dep.m_mesh_data.mesh_node(i+einfo.m_start)-center);

              m_shape_interpolant[new_shape_element.m_my_start_index+i] 
                = shapeval;
              el_integral += shapeval * m_dep.m_integral_weights[i];
            }
            m_integral += el_integral*einfo.m_jacobian;

            m_particle_shape_elements.push_back(new_shape_element);
          }




          void end_particle(
              particle_number pn,
              const double charge)
          {
            m_dstate.m_stats.m_el_per_particle_stats.add(
                m_particle_shape_elements.size());

            if (m_integral == 0)
            {
              WARN(boost::str(boost::format(
                      "deposited particle mass is zero (particle %d, #elements=%d)") 
                    % pn 
                    % m_particle_shape_elements.size()));
              return;
            }

            const double scale = charge/m_integral;
            m_dstate.m_stats.m_normalization_stats.add(scale);

            m_target.begin_particle(pn);
            BOOST_FOREACH(const shape_element &sel, m_particle_shape_elements)
              m_target.add_shape_on_element(sel.m_el_number,
                  sel.m_global_start_index,
                  scale * subrange(m_shape_interpolant, 
                    sel.m_my_start_index,
                    sel.m_my_start_index+sel.m_el_length));
            m_target.end_particle(pn);
          }
      };




    public:
      // member data --------------------------------------------------------
      ShapeFunction m_shape_function;
      dyn_vector m_integral_weights;
      const mesh_data &m_mesh_data;




      normalized_shape_function_depositor(
          const mesh_data &md,
          const py_matrix &mass_matrix)
        : m_mesh_data(md)
      { 
        m_integral_weights = prod(mass_matrix, 
            boost::numeric::ublas::scalar_vector<double>
            (mass_matrix.size1(), 1));
      }




      template<class Target>
      void deposit_densities_on_target(
          depositor_state &ds,
          const particle_state &ps,
          Target &tgt, boost::python::slice const &pslice) const
      {
        normalizing_element_target<Target> norm_tgt(
            *this, ds, ps, m_integral_weights, tgt);

        element_finder el_finder(m_mesh_data);

        FOR_ALL_SLICE_INDICES(particle_number, pn, 
            pslice, ps.particle_count)
        {
          norm_tgt.begin_particle();
          el_finder(ps, norm_tgt, pn, m_shape_function.radius());
          norm_tgt.end_particle(pn, ps.charges[pn]);
        }
      }
  };
}




#endif
