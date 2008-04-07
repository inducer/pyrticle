// Pyrticle - Particle in Cell in Python
// Grid-based reconstruction
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




#ifndef _AFAYYTAA_PYRTICLE_REC_GRID_HPP_INCLUDED
#define _AFAYYTAA_PYRTICLE_REC_GRID_HPP_INCLUDED




#include <vector>
#include <pyublas/numpy.hpp>
#include <boost/foreach.hpp>
#include "rec_shape.hpp"
#include <pyublas/unary_op.hpp>




namespace pyrticle
{
  struct grid_reconstructor 
  {
    template <class PICAlgorithm>
    class type : public reconstructor_base
    {
      public:
        typedef unsigned structured_node_number;

        struct brick
        {
          public:
            structured_node_number m_start_index;
            bounded_vector m_stepwidths;
            bounded_vector m_origin;

            /** This is the number of points in each dimension. If it is 2, then the
             * indices 0 and 1 exist in that dimension.
             */
            bounded_int_vector m_dimensions;

            unsigned node_count()
            {
              unsigned result = 1;
              BOOST_FOREACH(unsigned d, m_dimensions)
                result *= d;
              return result;
            }

            bounded_box bounding_box()
            {
              return std::make_pair(
                  m_origin,
                  m_origin + element_prod(m_dimensions, m_stepwidths) - m_stepwidths
                  );
            }

          private:
            struct int_floor
            {
              typedef double value_type;
              typedef const double &argument_type;
              typedef int result_type;

              static result_type apply(argument_type x)
              {
                return int(floor(x));
              }
            };

            struct int_ceil
            {
              typedef double value_type;
              typedef const double &argument_type;
              typedef int result_type;

              static result_type apply(argument_type x)
              {
                return int(ceil(x));
              }
            };
          public:
            bounded_int_box index_range(const bounded_box &bbox)
            {
              return bounded_int_box(
                  pyublas::unary_op<int_floor>::apply(
                    element_div(bbox.first-m_origin, m_stepwidths)),
                  pyublas::unary_op<int_ceil>::apply(
                    element_div(bbox.second-m_origin, m_stepwidths))
                  );
            }
        };

        std::vector<brick> m_bricks;

        typedef unsigned brick_number;
        typedef unsigned brick_node_number;

        struct element_on_grid
        {
          std::vector<structured_node_number> m_structured_nodes;

          /** The interpolant matrix maps the values (in-order) on the element
           * to the structured grid values at indices m_structured_nodes.
           *
           * The general assumption is that #(element values) < #(structured values),
           * but this may be violated without much harm.
           */
          hedge::matrix m_interpolant_pseudo_inv;
        };

        // member data --------------------------------------------------------
        shape_function   m_shape_function;

        // setup interface ------------------------------------------------------
        void add_brick(
            hedge::vector stepwidths,
            hedge::vector origin,
            pyublas::numpy_vector<unsigned> dims)
        {
          brick new_brick;

          if (m_bricks.size())
            new_brick.m_start_index = m_bricks.back().m_start_index +
              m_bricks.back().node_count();
          else
            new_brick.m_start_index = 0;

          new_brick.m_stepwidths = stepwidths;
          new_brick.m_origin = origin;
          new_brick.m_dimensions = dims;

          m_bricks.push_back(new_brick);
        }




        void commit_bricks()
        {
          BOOST_FOREACH(mesh_data::element_info &el, 
              CONST_PIC_THIS->m_mesh_data.m_element_info)
          {
            element_on_grid eog;
            bounded_box el_bbox = CONST_PIC_THIS->m_mesh_data.element_bounding_box(el.m_id);

            std::vector<bounded_vector> point_coordinates;
            BOOST_FOREACH(brick &brk, m_bricks)
            {
              bounded_box brick_bbox = brk.bounding_box();
              bool does_intersect;
              bounded_box el_brick_box = intersect(brick_bbox, el_bbox, &does_intersect);

              if (!does_intersect)
                continue;















            }
          }
        }




        // public interface -----------------------------------------------------
        static const std::string get_name()
        { return "Grid"; }

        void perform_reconstructor_upkeep()
        { }

        void reconstruct_densities(
            hedge::vector rho, 
            hedge::vector j,
            const hedge::vector &velocities)
        {

        }

        void reconstruct_j(hedge::vector j, const hedge::vector &velocities)
        {

        }

        void reconstruct_rho(hedge::vector rho)
        {

        }
    };
  };
}




#endif

