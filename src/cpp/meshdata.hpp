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




#include <numeric>
#include <hedge/base.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include "tools.hpp"




namespace pyrticle 
{
  template <class VecType>
  inline
  const bool is_in_unit_simplex(const VecType &unit_coords, double tolerance=1e-10)
  {
    BOOST_FOREACH(typename VecType::value_type ri, unit_coords)
      if (ri < -1-tolerance)
        return false;

    return std::accumulate(unit_coords.begin(), unit_coords.end(), 
        (double) 0) <= -(signed(unit_coords.size())-2)+tolerance;
  }




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
        bounded_vector        m_normal;
        element_number        m_neighbor;
        axis_number           m_neighbor_periodicity_axis;

        /** The equation for the hyperplane containing this face is
         * m_normal * x = m_face_plane_eqn_rhs.
         */
        double                m_face_plane_eqn_rhs;

        /** These two specify a bounding circle on the face. */
        bounded_vector        m_face_centroid;
        double                m_face_radius_from_centroid;
      };

      struct element_info
      {
        element_number                 m_id;
        hedge::affine_map              m_inverse_map;
        double                         m_norm_forward_map;
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

    private:
      py_vector m_mesh_vertices, m_mesh_nodes;

    public:
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

      void set_vertices(py_vector v)
      { m_mesh_vertices = v; }
      void set_nodes(py_vector n)
      { m_mesh_nodes = n; }

      static element_number get_INVALID_ELEMENT() { return INVALID_ELEMENT; }
      static axis_number get_INVALID_AXIS() { return INVALID_AXIS; }




      // operations -----------------------------------------------------------
      unsigned node_count() const
      { return m_mesh_nodes.size()/m_dimensions; }

      unsigned vertex_count() const
      { return m_mesh_vertices.size()/m_dimensions; }

      typedef boost::numeric::ublas::vector_range<const py_vector> const_mesh_node_type;
      typedef boost::numeric::ublas::vector_range<const py_vector> const_mesh_vertex_type;
      typedef boost::numeric::ublas::vector_range<py_vector> mesh_node_type;
      typedef boost::numeric::ublas::vector_range<py_vector> mesh_vertex_type;

      const_mesh_vertex_type mesh_vertex(vertex_number vn) const
      { return subrange(m_mesh_vertices, vn*m_dimensions, (vn+1)*m_dimensions); }

      mesh_vertex_type mesh_vertex(vertex_number vn)
      { return subrange(m_mesh_vertices, vn*m_dimensions, (vn+1)*m_dimensions); }

      const_mesh_node_type mesh_node(node_number nn) const
      { return subrange(m_mesh_nodes, nn*m_dimensions, (nn+1)*m_dimensions); }

      mesh_node_type mesh_node(node_number nn)
      { return subrange(m_mesh_nodes, nn*m_dimensions, (nn+1)*m_dimensions); }

      bounded_vector element_centroid(element_number en) const
      {
        const element_info &ei = m_element_info[en];

        bounded_vector result(mesh_vertex(ei.m_vertices[0]));

        for (unsigned i = 1; i < ei.m_vertices.size(); ++i)
          result += mesh_vertex(ei.m_vertices[i]);

        result /= ei.m_vertices.size();
        return result;
      }

      bounded_box element_bounding_box(element_number en) const
      {
        const element_info &ei = m_element_info[en];
        bounded_vector
          min(mesh_vertex(ei.m_vertices[0])), 
          max(mesh_vertex(ei.m_vertices[0]));

        for (unsigned vi = 0; vi < ei.m_vertices.size(); ++vi)
        {
          const_mesh_vertex_type vtx = mesh_vertex(ei.m_vertices[vi]);
          for (unsigned i = 0; i < m_dimensions; ++i)
          {
            if (vtx[i] < min[i]) min[i] = vtx[i];
            if (vtx[i] > max[i]) max[i] = vtx[i];
          }
        }

        return std::make_pair(min, max);
      }

      template <class VecType>
      const bool is_in_element(element_number en, const VecType &pt, double tolerance=1e-10) const
      {
        return is_in_unit_simplex(m_element_info[en].m_inverse_map(pt), tolerance);
      }

      template <class VecType>
      const element_number find_containing_element(const VecType &pt) const
      {
        BOOST_FOREACH(const element_info &el, m_element_info)
          if (is_in_unit_simplex(el.m_inverse_map(pt)))
            return el.m_id;
        return INVALID_ELEMENT;
      }
  };
}




#endif
