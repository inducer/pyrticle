// Pyrticle - Particle in Cell in Python
// Reconstruction based on shape functions
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





#ifndef _BADFJAH_PYRTICLE_REC_SHAPE_HPP_INCLUDED
#define _BADFJAH_PYRTICLE_REC_SHAPE_HPP_INCLUDED




#include <boost/ref.hpp>
#include <boost/foreach.hpp>
#include <boost/unordered_set.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/typeof/std/utility.hpp>
#include "tools.hpp"
#include "bases.hpp"
#include "meshdata.hpp"
#include "pic_algorithm.hpp"




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




  struct shape_function_reconstructor
  {
    template <class PICAlgorithm>
    class type : public target_reconstructor_base<PICAlgorithm>
    {
      private:
        template <class Target>
        class element_target
        {
          private:
            const PICAlgorithm &m_pic_algorithm;
            const double      m_charge;
            Target            &m_target;
            
          public:
            element_target(PICAlgorithm &pic_algorithm, double charge, Target &tgt)
              : m_pic_algorithm(pic_algorithm), m_charge(charge), m_target(tgt)
            { }

            void add_shape_on_element(
                const bounded_vector &center,
                const mesh_data::element_number en
                ) const
            {
              const mesh_data::element_info &einfo(
                  m_pic_algorithm.m_mesh_data.m_element_info[en]);

              const unsigned el_length = einfo.m_end-einfo.m_start;
              dyn_vector el_rho(el_length);

              for (unsigned i = 0; i < el_length; i++)
                el_rho[i] = 
                  m_charge * m_pic_algorithm.m_shape_function(
                      m_pic_algorithm.m_mesh_data.mesh_node(einfo.m_start+i) 
                      - center);

              m_target.add_shape_on_element(en, einfo.m_start, el_rho);
            }
        };

      public:
        // member data --------------------------------------------------------
        shape_function   m_shape_function;

        // public interface ---------------------------------------------------
        static const std::string get_name()
        { return "Shape"; }




        template<class Target>
        void reconstruct_densities_on_target(Target &tgt,
            boost::python::slice const &pslice)
        {
          typename PICAlgorithm::element_finder el_finder(*CONST_PIC_THIS);

          FOR_ALL_SLICE_INDICES(particle_number, pn, 
              pslice, CONST_PIC_THIS->m_particle_count)
          {
            element_target<Target> el_target(
                *CONST_PIC_THIS, 
                CONST_PIC_THIS->m_charges[pn], tgt);

            tgt.begin_particle(pn);
            el_finder(el_target, pn, m_shape_function.radius());
            tgt.end_particle(pn);
          }
        }
    };
  };
}




#endif
