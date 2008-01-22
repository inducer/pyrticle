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
#include <boost/array.hpp>
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
          const mesh_data::element_info *m_element_info;
          boost::array<mesh_data::element_number, type::max_faces> m_connections;
          unsigned m_start_index;

          active_element()
            : m_element_info(0)
          {
            for (unsigned i = 0; i < type::max_faces; i++)
              m_connections[i] = mesh_data::INVALID_ELEMENT;
          }
        };

        struct advected_particle
        {
          std::vector<active_element>   m_elements;
          double                        m_radius;
        };




        // member data --------------------------------------------------------
        unsigned                        m_dofs_per_element;
        unsigned                        m_active_elements;
        std::vector<unsigned>           m_freelist;

        std::vector<advected_particle>  m_advected_particles;

        hedge::vector                   m_rho;
        hedge::vector                   m_j; // concatenated components of j




        // public interface ---------------------------------------------------
        type()
          : m_dofs_per_element(0), m_active_elements(0)
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




        void perform_reconstructor_upkeep()
        {
          // retire empty particle subelements 
        }




        // initialization -----------------------------------------------------
        void setup_advection_reconstructor(unsigned dofs_per_element)
        {
          m_dofs_per_element = dofs_per_element;
          resize_state(m_dofs_per_element * 1024);
        }




        // vectors space administration ---------------------------------------
        /* Each element occupies a certain index range in the global
         * state vectors m_rho and m_j (as well as elsewhere). These
         * functions perform allocation and deallocation of space in
         * these vectors.
         */

        void resize_state(unsigned new_size)
        {
          unsigned old_size = m_rho.size();
          unsigned copy_size = std::min(new_size, old_size);

          hedge::vector new_rho(new_size);
          subrange(new_rho, 0, copy_size) = subrange(m_rho, 0, copy_size);
          new_rho.swap(m_rho);

          hedge::vector new_j(new_size * PICAlgorithm::dimensions_pos);
          for (unsigned i = 0; i < PICAlgorithm::dimensions_pos; i++)
            subrange(new_j, i*new_size, i*new_size+copy_size) = 
              subrange(m_rho, i*old_size, i*old_size+copy_size);
          new_j.swap(m_j);
        }

        /** Allocate a space for a new element in the state vector, return
         * the start index.
         */
        unsigned allocate_element()
        {
          if (m_dofs_per_element == 0)
            throw std::runtime_error("tried to allocate element on uninitialized advection reconstructor");

          if (m_freelist.size())
          {
            unsigned result = m_freelist.back();
            m_freelist.pop_back();
            return result*m_dofs_per_element;
          }

          // we're all full, no gaps available.
          // return the past-end spot in the array, reallocate if necessary.
          unsigned avl_space = m_rho.size() / m_dofs_per_element;

          if (m_active_elements == avl_space)
            resize_state(2*m_rho.size());

          return (m_active_elements++)*m_dofs_per_element;
        }

        void deallocate_element(unsigned start_index)
        {
          m_freelist.push_back(start_index/m_dofs_per_element);
          --m_active_elements;
        }

        // particle construction ----------------------------------------------
        void add_shape_on_element(
            advected_particle &new_particle,
            const hedge::vector &center,
            const mesh_data::element_number en
            )
        {
          const mesh_data::element_info &einfo(
              CONST_PIC_THIS->m_mesh_data.m_element_info[en]);

          active_element new_element;
          new_element.m_element_info = &einfo;
          unsigned start = new_element.m_start_index = allocate_element();

          shape_function sf(new_particle.m_radius, PICAlgorithm::dimensions_pos);

          for (unsigned i = einfo.m_start; i < einfo.m_end; i++)
            m_rho[start+i] =
                sf(CONST_PIC_THIS->m_mesh_data.m_nodes[i]-center);
          
          new_particle.m_elements.push_back(new_element);
        }




        void add_advection_particle(double radius, particle_number pn)
        {
          if (pn != m_advected_particles.size())
            throw std::runtime_error("advected particle added out of sequence");

          advected_particle new_particle;
          new_particle.m_radius = radius;

          add_shape(new_particle, pn, radius);

          m_advected_particles.push_back(new_particle);
        }








        hedge::vector get_advection_particle_rhs()
        {
          return hedge::vector();
        }
        



        void apply_advection_particle_rhs(hedge::vector const &rhs)
        {

        }
    };
  };
}




#endif
