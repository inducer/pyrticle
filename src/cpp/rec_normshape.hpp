// Pyrticle - Particle in Cell in Python
// Reconstruction based on normalized advected shapes
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




#ifndef _AJFYAA_PYRTICLE_REC_NORMSHAPE_HPP_INCLUDED
#define _AJFYAA_PYRTICLE_REC_NORMSHAPE_HPP_INCLUDED




#include <boost/foreach.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/typeof/std/utility.hpp>
#include "tools.hpp"
#include "bases.hpp"
#include "meshdata.hpp"
#include "rec_shape.hpp"




namespace pyrticle 
{
  // Note: this does not adhere to the ReconstructionTarget protocol.
  class normalizing_target
  {
    public:
      const hedge::mesh_data &m_mesh_data;
      hedge::vector &m_rho_target_vector;
      hedge::vector &m_j_target_vector;

      hedge::vector m_shape_interpolant;
      unsigned m_used_shape_dofs;

      struct shape_element {
        unsigned m_my_start_index;
        node_number m_global_start_index;
        element_number m_element_number;

        shape_element(unsigned my_start_index, element_number en)
          : m_my_start_index(my_start_index), m_element_number(en)
        { }
      };

      std::vector<shape_element>        m_particle_shape_elements;
      double                            m_integral;

      normalizing_target(
          const mesh_data &md,
          hedge::vector &rho_target_vector, 
          hedge::vector &j_target_vector)
        : m_rho_target_vector(rho_target_vector), 
        m_mesh_data(md), m_j_target_vector(j_target_vector),
        m_shape_interpolant(100)
      { 
        m_target_vector.clear();
      }

      void begin_particle()
      {
        m_used_shape_dofs = 0;
        m_particle_shape_elements.clear();
        m_integral = 0;
      }

      void add_shape_on_element(
          const mesh_data::element_info &einfo,
          const hedge::vector &center,
          const shape_function sf
          )
      {
        unsigned element_length = einfo.m_end-einfo.m_start;

        // make sure we have enough interpolant dofs available
        while (m_used_shape_dofs + element_length > m_shape_interpolant.size())
        {
          hedge::vector new_shape_interpolant(2*m_shape_interpolant.size());
          subrange(new_shape_interpolant, 0, m_shape_interpolant.size()) =
            m_shape_interpolant;
          new_shape_interpolant.swap(m_shape_interpolant);
        }

        shape_element new_shape_element(m_used_shape_dofs, einfo.m_id);
        m_used_shape_dofs += element_length;

        for (unsigned i = 0; i < element_length; i++)
          m_shape_interpolant[new_shape_element.m_my_start_index+i] 
            = sf(m_mesh_data.m_nodes[i]-center);
        
        m_particle_shape_elements.push_back(new_shape_element);
      }

      template <class VelocityVectorExpression>
      void end_particle(
          const double charge,
          const VelocityVectorExpression &p_velocity,
          )
      {

      }

  };




  struct normalized_shape_function_reconstructor
  {
    template <class PICAlgorithm>
    class type : 
      public shape_element_finder<PICAlgorithm>, 
      public reconstructor_base
    {
      public:
        // member data --------------------------------------------------------
        shape_function   m_shape_function;

        // public interface ---------------------------------------------------
        static const char *get_name()
        { return "NormShape"; }

        void set_radius(double radius)
        {
          m_shape_function = shape_function(
              radius, CONST_PIC_THIS->m_mesh_data.m_dimensions);
        }




        void reconstruct_densities(
            hedge::vector &rho, 
            hedge::vector &j,
            const hedge::vector &velocities)
        {
          if (rho.size() != CONST_PIC_THIS->m_mesh_data.m_nodes.size())
            throw std::runtime_error("rho field does not have the correct size");
          if (j.size() != CONST_PIC_THIS->m_mesh_data.m_nodes.size() *
              CONST_PIC_THIS->get_dimensions_velocity())
            throw std::runtime_error("j field does not have the correct size");

          rho_reconstruction_target rho_tgt(rho, CONST_PIC_THIS->m_charges);
          j_reconstruction_target<PICAlgorithm::dimensions_velocity> j_tgt(
              j, CONST_PIC_THIS->m_charges, velocities);

          chained_reconstruction_target
            <rho_reconstruction_target, 
            j_reconstruction_target<PICAlgorithm::dimensions_velocity> >
            tgt(rho_tgt, j_tgt);
          reconstruct_densities_on_target(tgt);

          rho = rho_tgt.result();
        }




        void reconstruct_j(hedge::vector &j, const hedge::vector &velocities)
        {
          if (j.size() != CONST_PIC_THIS->m_mesh_data.m_nodes.size() *
              CONST_PIC_THIS->get_dimensions_velocity())
            throw std::runtime_error("j field does not have the correct size");

          j_reconstruction_target<PICAlgorithm::dimensions_velocity> j_tgt(
              j, CONST_PIC_THIS->m_charges, velocities);

          reconstruct_densities_on_target(j_tgt);
        }




        void reconstruct_rho(hedge::vector &rho)
        {
          if
          if (rho.size() != CONST_PIC_THIS->m_mesh_data.m_nodes.size())
            throw std::runtime_error("rho field does not have the correct size");
          rho_reconstruction_target rho_tgt(rho, CONST_PIC_THIS->m_charges);

          reconstruct_densities_on_target(rho_tgt);
        }




        void perform_reconstructor_upkeep()
        { }




        // inner workings -----------------------------------------------------
        template <class Target>
        void add_shape_on_element(
            Target &tgt, 
            const hedge::vector &center,
            const mesh_data::element_number en
            ) const
        {
          const mesh_data::element_info &einfo(
              CONST_PIC_THIS->m_mesh_data.m_element_info[en]);
        }




        template<class Target>
        void reconstruct_densities_on_target(Target &tgt)
        {
          for (particle_number pn = 0; pn < CONST_PIC_THIS->m_particle_count; ++pn)
          {
            tgt.begin_particle(pn);
            add_shape(tgt, pn, m_shape_function->radius());
          }
        }
    };

  };

}




#endif

