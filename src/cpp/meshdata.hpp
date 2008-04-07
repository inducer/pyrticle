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
      struct face_info
      {
        hedge::vector         m_normal;
        element_number        m_neighbor;
        axis_number           m_neighbor_periodicity_axis;

        /** The equation for the hyperplane containing this face is
         * m_normal * x = m_face_plane_eqn_rhs.
         */
        hedge::vector::value_type m_face_plane_eqn_rhs;

        /** These two specify a bounding circle on the face. */
        hedge::vector         m_face_centroid;
        hedge::vector::value_type m_face_radius_from_centroid;
      };

      struct element_info
      {
        element_number                 m_id;
        hedge::affine_map              m_inverse_map;
        double                         m_jacobian;

        unsigned                       m_start, m_end;

        vtx_id_vector                  m_vertices;

        std::vector<face_info>      m_faces;

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
      bounded_vector element_centroid(element_number en) const
      {
        const element_info &ei = m_element_info[en];

        bounded_vector result(m_vertices[ei.m_vertices[0]]);

        for (unsigned i = 1; i < ei.m_vertices.size(); ++i)
          result += m_vertices[ei.m_vertices[i]];

        result /= ei.m_vertices.size();
        return result;
      }

      bounded_box element_bounding_box(element_number en) const
      {
        const element_info &ei = m_element_info[en];
        bounded_vector
          min(m_vertices[ei.m_vertices[0]]), 
          max(m_vertices[ei.m_vertices[0]]);

        for (unsigned vi = 0; vi < ei.m_vertices.size(); ++vi)
        {
          const hedge::vector &vtx = m_vertices[ei.m_vertices[vi]];
          for (unsigned i = 0; i < m_dimensions; ++i)
          {
            if (vtx[i] < min[i]) min[i] = vtx[i];
            if (vtx[i] > max[i]) max[i] = vtx[i];
          }
        }

        return std::make_pair(min, max);
      }

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
