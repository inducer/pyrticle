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
#include <boost/format.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/typeof/std/utility.hpp>
#include "tools.hpp"
#include "bases.hpp"
#include "meshdata.hpp"
#include "rec_shape.hpp"




namespace pyrticle 
{
  /* This does not adhere to the ReconstructionTarget protocol.
   * Instead, it is a layer that sits between the 
   * normalized_shape_function_reconstructor, below, and the 
   * real ReconstructionTarget.
   */
  template <class Target>
  class normalizing_target
  {
    public:
      const mesh_data &m_mesh_data;
      const hedge::vector &m_integral_weights;
      Target m_target;
      hedge::vector m_shape_interpolant;
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

      normalizing_target(
         const mesh_data &mesh_data,
         const hedge::vector &integral_weights,
         Target &target)
        : 
          m_mesh_data(mesh_data),
          m_integral_weights(integral_weights),
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

        shape_element new_shape_element(
            einfo.m_id,
            element_length,
            m_used_shape_dofs, 
            einfo.m_start);
        m_used_shape_dofs += element_length;

        double el_integral = 0;
        for (unsigned i = 0; i < element_length; i++)
        {
          double shapeval = sf(m_mesh_data.m_nodes[i]-center);
          m_shape_interpolant[new_shape_element.m_my_start_index+i] 
            = shapeval;
          el_integral += shapeval * m_integral_weights[i];
        }
        m_integral += el_integral;
        
        m_particle_shape_elements.push_back(new_shape_element);
      }




      void end_particle(
          particle_number pn,
          const double charge)
      {
        if (m_integral == 0)
            throw std::runtime_error(
                str(boost::format("reconstructed particle mass is zero"
                  "(particle %d, #elements=%d)") 
                  % pn 
                  % m_particle_shape_elements.size()).c_str());

        const double scale = charge/m_integral;

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




  struct normalized_shape_function_reconstructor
  {
    template <class PICAlgorithm>
    class type : 
      public shape_element_finder<PICAlgorithm>, 
      public reconstructor_base
    {
      public:
        // member data --------------------------------------------------------
        shape_function                  m_shape_function;
        hedge::vector                   m_integral_weights;



      
        // initialization -----------------------------------------------------
        void setup_normalized_shape_reconstructor(const hedge::matrix &mass_matrix)
        {
          m_integral_weights = prod(mass_matrix, 
              boost::numeric::ublas::scalar_vector<double>
              (mass_matrix.size1(), 1));
        }




        // public interface ---------------------------------------------------
        static const char *get_name()
        { return "NormShape"; }

        void set_radius(double radius)
        {
          m_shape_function = shape_function(
              radius, CONST_PIC_THIS->m_mesh_data.m_dimensions);
        }




        void perform_reconstructor_upkeep()
        { }




        template<class Target>
        void reconstruct_densities_on_target(Target &tgt)
        {
          normalizing_target<Target> norm_tgt(
              CONST_PIC_THIS->m_mesh_data,
              m_integral_weights,
              tgt);

          for (particle_number pn = 0; pn < CONST_PIC_THIS->m_particle_count; ++pn)
          {
            norm_tgt.begin_particle();
            add_shape(norm_tgt, pn, m_shape_function.radius());
            norm_tgt.end_particle(pn, CONST_PIC_THIS->m_charges[pn]);
          }
        }




        // inner workings -----------------------------------------------------
        template <class NormTarget>
        void add_shape_on_element(
            NormTarget &norm_tgt, 
            const hedge::vector &center,
            const mesh_data::element_number en
            ) const
        {
          norm_tgt.add_shape_on_element(
              CONST_PIC_THIS->m_mesh_data.m_element_info[en],
              center,
              m_shape_function);
        }
    };
  };
}




#endif

