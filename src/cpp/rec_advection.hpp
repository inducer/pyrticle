// Pyrticle - Particle in Cell in Python
// Reconstruction based on advected shapes
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





#ifndef _AFAYYTAA_PYRTICLE_REC_ADVECTION_HPP_INCLUDED
#define _AFAYYTAA_PYRTICLE_REC_ADVECTION_HPP_INCLUDED




#include <vector>
#include <boost/foreach.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/typeof/std/utility.hpp>
#include "tools.hpp"
#include "meshdata.hpp"
#include "rec_shape.hpp"




namespace pyrticle
{
  struct advection_reconstructor 
  {
    template <class PICAlgorithm>
    class type : public shape_element_finder<PICAlgorithm>
    {
      public:
        static const unsigned max_faces = 4;

        // member types -------------------------------------------------------
        struct active_element
        {
          mesh_data::element_number     m_element;
          hedge::vector                 m_element_values;
          hedge::vector                 m_axis_rhs[PICAlgorithm::dimensions_pos];

          mesh_data::element_number     m_outbound_connections[type::max_faces];
        };

        struct advected_particle
        {
          std::vector<active_element>   m_elements;
          double                        m_radius;
        };

        // member data --------------------------------------------------------

        std::vector<advected_particle>  m_advected_particles;
        unsigned                        m_advection_elements;

        // publicized interface -----------------------------------------------
        type()
          : m_advection_elements(0)
        { }




        static const char *get_name()
        { return "Advection"; }




        void reconstruct_densities(
            hedge::vector &rho, 
            hedge::vector &j,
            const hedge::vector &velocities)
        {

        }




        void reconstruct_j(hedge::vector &j, const hedge::vector &velocities)
        {

        }




        void reconstruct_rho(hedge::vector &rho)
        {

        }




        // particle construction ----------------------------------------------
        template <class Target>
        void add_shape_on_element(
            advected_particle &new_particle,
            const hedge::vector &center,
            const mesh_data::element_number en
            ) const
        {
          const mesh_data::element_info &einfo(
              CONST_PIC_THIS->m_mesh_data.m_element_info[en]);
          /*
          for (unsigned i = einfo.m_start; i < einfo.m_end; i++)
            tgt.add_shape_at_point(i, 
                (*m_shape_function)(CONST_PIC_THIS->m_mesh_data.m_nodes[i]-center));
                */
        }




        void add_advected_particle(double radius, particle_number pn)
        {
          if (pn != m_advected_particles.size())
            throw std::runtime_error("advected particle added out of sequence");

          advected_particle new_particle;
          new_particle.m_radius = radius;

          add_shape(new_particle, pn);

          m_advected_particles.push_back(new_particle);
        }








        hedge::vector get_adv_particle_rhs()
        {
          return hedge::vector();
        }
        



        void apply_adv_particle_rhs(hedge::vector const &rhs)
        {

        }
    };
  };
}




#endif
