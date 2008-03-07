// Pyrticle - Particle in Cell in Python
// Stuff we need to store about the field solver's mesh
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





#ifndef _BADFJAH_PYRTICLE_MESHDATA_HPP_INCLUDED
#define _BADFJAH_PYRTICLE_MESHDATA_HPP_INCLUDED




#include <hedge/base.hpp>
#include "tools.hpp"




namespace pyrticle 
{
  const bool is_in_unit_simplex(const hedge::vector &unit_coords);




  class mesh_data : boost::noncopyable
  {
    public:
      typedef hedge::element_number element_number;
      typedef hedge::face_number face_number;
      typedef hedge::vertex_number vertex_number;
      typedef hedge::node_index node_number;
      typedef hedge::el_id_vector el_id_vector;
      typedef hedge::vtx_id_vector vtx_id_vector;
      typedef hedge::el_face el_face;

      typedef unsigned axis_number;

      static const element_number INVALID_ELEMENT = hedge::INVALID_ELEMENT;
      static const vertex_number INVALID_VERTEX = hedge::INVALID_VERTEX;
      static const axis_number INVALID_AXIS = UINT_MAX;




      // data structures ------------------------------------------------------
      struct element_info
      {
        element_number                 m_id;
        hedge::affine_map              m_inverse_map;
        double                         m_jacobian;

        unsigned                       m_start, m_end;

        vtx_id_vector                  m_vertices;

        /** The indices for the following lists match up: say at matching
         * indices you find normal n and element index i, then n is the normal
         * of the face leading to element i. By corollary, if there is no 
         * neighbor past that face, INVALID_ELEMENT is stored in that entry.
         */

        std::vector<hedge::vector>      m_normals;
        std::vector<element_number>     m_neighbors;
        std::vector<axis_number>        m_neighbor_periodicity_axes;

        hedge::vector centroid(const std::vector<hedge::vector> &vertex_points) const
        {
          hedge::vector result(vertex_points[m_vertices[0]]);

          for (unsigned i = 1; i < m_vertices.size(); ++i)
            result += vertex_points[m_vertices[i]];

          result /= m_vertices.size();
          return result;
        }
      };


      struct periodicity_axis
      {
        double                  m_min, m_max;
      };

      // data members ---------------------------------------------------------
      const unsigned m_dimensions;

      std::vector<element_info> m_element_info;
      std::vector<hedge::vector> m_vertices, m_nodes;

      /** The following three encode the vertex-adjacent elements in a sort
       * of Compressed-Row-Storage format. */
      std::vector<unsigned> m_vertex_adj_element_starts;
      el_id_vector m_vertex_adj_elements;
      std::vector<axis_number> m_vertex_adj_periodicity_axes;

      std::vector<periodicity_axis> m_periodicities;




      // setup ----------------------------------------------------------------
      mesh_data(unsigned dimensions)
        : m_dimensions(dimensions)
      { }




      static element_number get_INVALID_ELEMENT() { return INVALID_ELEMENT; }
      static axis_number get_INVALID_AXIS() { return INVALID_AXIS; }




      // operations -----------------------------------------------------------
      const bool is_in_element(element_number en, const hedge::vector &pt) const
      {
        const element_info &el(m_element_info[en]);
        hedge::vector uc = el.m_inverse_map(pt);
        return is_in_unit_simplex(uc);
      }

      const element_number find_containing_element(
          const hedge::vector &pt) const;
  };
}




#endif
