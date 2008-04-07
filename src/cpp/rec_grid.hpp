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
#include "rec_shape.hpp"




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
          structured_node_number m_start_index;
          hedge::vector m_stepwidths;
          hedge::vector m_origin;
          pyublas::numpy_vector<unsigned> m_size;
        };

        std::vector<brick> m_bricks;

        typedef unsigned brick_number;
        typedef unsigned brick_node_number;

        struct element_info
        {
          std::vector<structured_node_numbers> m_structured_nodes;

          /** The interpolant matrix maps the values (in-order) on the element
           * to the structured grid values at indices m_structured_nodes.
           *
           * The general assumption is that #(element values) < #(structured values),
           * but this may be violated without much harm.
           */
          hedge::matrix m_interpolant_pseudo_inv;
        };


      // setup interface ------------------------------------------------------

      // public interface -----------------------------------------------------
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

