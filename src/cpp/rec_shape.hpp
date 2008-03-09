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
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/typeof/std/utility.hpp>
#include "tools.hpp"
#include "bases.hpp"
#include "meshdata.hpp"




namespace pyrticle 
{
  class shape_function
  {
    public:
      shape_function(double radius=1, unsigned dimensions=1, double alpha=2);

      const double operator()(const hedge::vector &r) const
      {
        double r_squared = inner_prod(r, r);
        if (r_squared > m_l_squared)
          return 0;
        else
          return m_normalizer * pow(m_l-r_squared/m_l, m_alpha);
      }

      const double radius() const
      { return m_l; }

    private:
      double m_normalizer;
      double m_alpha;
      double m_l, m_l_squared;
  };




  inline
  bool is_not_near_vertex(hedge::vector const &unit_pt)
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




  template <class PICAlgorithm>
  class shape_element_finder
  {
    private:
      const PICAlgorithm &m_pic_algorithm;

    public:
      shape_element_finder(const PICAlgorithm &pa)
        : m_pic_algorithm(pa)
      { }

      template <class ElementTarget>
      void add_shape_by_neighbors(
          ElementTarget &target,
          const hedge::vector &pos,
          const mesh_data::element_info &einfo,
          double radius)
      {
        target.add_shape_on_element(pos, einfo.m_id);

        for (unsigned i = 0; i < einfo.m_neighbors.size(); ++i)
        {
          const mesh_data::element_number en = einfo.m_neighbors[i];

          if (en != mesh_data::INVALID_ELEMENT)
          {
            const mesh_data::axis_number per_axis = 
              einfo.m_neighbor_periodicity_axes[i];

            if (per_axis == mesh_data::INVALID_AXIS)
              target.add_shape_on_element(pos, en);
            else
            {
              hedge::vector pos2(pos);
              const mesh_data::periodicity_axis &pa = 
                m_pic_algorithm.m_mesh_data.m_periodicities[per_axis];

              if (pos[per_axis] - radius < pa.m_min)
              {
                pos2[per_axis] += (pa.m_max-pa.m_min);
                target.add_shape_on_element(pos2, en);
              }
              if (pos[per_axis] + radius > pa.m_max)
              {
                pos2[per_axis] -= (pa.m_max-pa.m_min);
                target.add_shape_on_element(pos2, en);
              }
            }
          }
        }
      }




      template <class ElementTarget>
      void add_shape_by_vertex(
          ElementTarget &target,
          const hedge::vector &pos,
          const mesh_data::element_info &einfo,
          double radius)
      {
        // find closest vertex
        mesh_data::vertex_number closest_vertex = 
          mesh_data::INVALID_VERTEX;
        double min_dist = std::numeric_limits<double>::infinity();

        BOOST_FOREACH(mesh_data::vertex_number vi, einfo.m_vertices)
        {
          double dist = norm_2(m_pic_algorithm.m_mesh_data.m_vertices[vi] - pos);
          if (dist < min_dist)
          {
            closest_vertex = vi;
            min_dist = dist;
          }
        }

        const mesh_data &md = m_pic_algorithm.m_mesh_data;

        // go through vertex-adjacent elements
        unsigned start = md.m_vertex_adj_element_starts[closest_vertex];
        unsigned stop = md.m_vertex_adj_element_starts[closest_vertex+1];

        for (unsigned i = start; i < stop; ++i)
        {
          const mesh_data::axis_number per_axis = 
            md.m_vertex_adj_periodicity_axes[i];
          const mesh_data::element_number en = 
            md.m_vertex_adj_elements[i];

          if (per_axis == mesh_data::INVALID_AXIS)
            target.add_shape_on_element(pos, en);
          else
          {
            const mesh_data::periodicity_axis &pa = md.m_periodicities[per_axis];

            hedge::vector pos2(pos);
            if (pos[per_axis] - radius < pa.m_min)
            {
              pos2[per_axis] += (pa.m_max - pa.m_min);
              target.add_shape_on_element(pos2, en);
            }
            if (pos[per_axis] + radius > pa.m_max)
            {
              pos2[per_axis] -= (pa.m_max - pa.m_min);
              target.add_shape_on_element(pos2, en);
            }
          }
        }
      }




      template <class ElementTarget>
      void operator()(ElementTarget &target, particle_number pn, double radius)
      {
        const unsigned dim = m_pic_algorithm.get_dimensions_pos();
        const hedge::vector pos = subrange(
            m_pic_algorithm.m_positions, pn*dim, (pn+1)*dim);
        const mesh_data::element_number containing_el = 
          m_pic_algorithm.m_containing_elements[pn];
        const mesh_data::element_info &einfo(
            m_pic_algorithm.m_mesh_data.m_element_info[containing_el]);

        // We're deciding between RULE A and RULE B below.
        // The decision is made by looking at the barycentric coordinates of
        // the center of the shape function. If that center has all barycentric
        // coordinates 1/2 or smaller, then use the faster RULE A.
        //
        // For this purpose, observe that the unit coordinates are the first
        // barycentric coordinates, and the remaining one is found by calculating
        // 1-sum(lambda_i)
        //
        // Note also that this is not purely a speed tradeoff. RULE B works fairly
        // well alone, but falls flat on its face for near the center of the hypotenuse
        // of a right triangle.

        // FIXME this assumes dimension_pos == dimension_mesh

        if (is_not_near_vertex(einfo.m_inverse_map(pos)))
        {
          // RULE A: we're far enough away from vertices, 
          //m_neighbor_shape_adds.tick();
          add_shape_by_neighbors(target, pos, einfo, radius);
        }
        else
        {
          // RULE B: we're close to a vertex, weight onto all elements
          // adjoining that vertex
          //m_vertex_shape_adds.tick();
          add_shape_by_vertex(target, pos, einfo, radius);
        }
      }
  };





  struct shape_function_reconstructor
  {
    template <class PICAlgorithm>
    class type : public reconstructor_base
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
                const hedge::vector &center,
                const mesh_data::element_number en
                ) const
            {
              const mesh_data::element_info &einfo(
                  m_pic_algorithm.m_mesh_data.m_element_info[en]);

              const unsigned el_length = einfo.m_end-einfo.m_start;
              hedge::vector el_rho(el_length);

              for (unsigned i = 0; i < el_length; i++)
                el_rho[i] = 
                  m_charge * m_pic_algorithm.m_shape_function(
                      m_pic_algorithm.m_mesh_data.m_nodes[einfo.m_start+i] 
                      - center);

              m_target.add_shape_on_element(en, einfo.m_start, el_rho);
            }
        };

      public:
        // member data --------------------------------------------------------
        shape_function   m_shape_function;

        // public interface ---------------------------------------------------
        static const char *get_name()
        { return "Shape"; }




        void perform_reconstructor_upkeep()
        { }




        template<class Target>
        void reconstruct_densities_on_target(Target &tgt)
        {
          shape_element_finder<PICAlgorithm> el_finder(*CONST_PIC_THIS);
          
          for (particle_number pn = 0; pn < CONST_PIC_THIS->m_particle_count; ++pn)
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
