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




#include "tools.hpp"




namespace pyrticle 
{
  const bool is_in_unit_simplex(const hedge::vector &unit_coords);




  class mesh_data : boost::noncopyable
  {
    public:
      // data structures ------------------------------------------------------
      typedef unsigned element_number;
      typedef unsigned vertex_number;
      typedef std::vector<element_number> el_id_vector;
      typedef std::vector<vertex_number> vtx_id_vector;

      static const element_number INVALID_ELEMENT = UINT_MAX;
      static const vertex_number INVALID_VERTEX = UINT_MAX;

      struct element_info
      {
        element_number                 m_id;
        hedge::affine_map              m_inverse_map;

        unsigned                       m_start, m_end;

        vtx_id_vector                  m_vertices;

        // the indices for the following two lists match up:
        // say at matching indices you find normal n
        // and element index i, then n is the normal
        // of the face leading to element i.

        std::vector<hedge::vector>      m_normals;
        std::vector<element_number>     m_neighbors;
      };


      struct periodicity_axis
      {
        unsigned                m_axis;
        double                  m_min, m_max, m_width;
      };

      // data members ---------------------------------------------------------
      const unsigned m_dimensions;

      std::vector<element_info> m_element_info;
      std::vector<hedge::vector> m_vertices, m_nodes;
      std::vector<el_id_vector> m_vertex_adj_elements;
      std::vector<periodicity_axis> m_periodicities;




      // setup ----------------------------------------------------------------
      mesh_data(unsigned dimensions)
        : m_dimensions(dimensions)
      { }




      static element_number get_INVALID_ELEMENT() { return INVALID_ELEMENT; }




      /*
      void add_local_discretization(python::list basis,
          const hedge::matrix &l_vdmt, 
          const hedge::matrix &u_vdmt, 
          const csr_matrix &p_vdmt)
      {
        local_discretization ldis;

        COPY_PY_LIST(ldis.m_basis, basis, const monomial_basis_function &);

        ldis.m_l_mon_vandermonde_t = l_vdmt;
        ldis.m_u_mon_vandermonde_t = u_vdmt;
        ldis.m_p_mon_vandermonde_t = p_vdmt;
        m_local_discretizations.push_back(ldis);
      }
      */




      /*
      void add_element(const hedge::affine_map &inverse_map, 
          unsigned start, unsigned end,
          python::object vertices,
          python::object normals, 
          python::object neighbors
          )
      {
        element_info ei;
        ei.m_id = m_element_info.size();
        ei.m_inverse_map = inverse_map;

        ei.m_start = start;
        ei.m_end = end;

        std::copy(
            python::stl_input_iterator<vertex_number>(vertices),
            python::stl_input_iterator<vertex_number>(),
            back_inserter(ei.m_vertices));
        std::copy(
            python::stl_input_iterator<const hedge::vector &>(normals),
            python::stl_input_iterator<const hedge::vector &>(),
            back_inserter(ei.m_normals));
        std::copy(
            python::stl_input_iterator<element_number>(neighbors),
            python::stl_input_iterator<element_number>(),
            back_inserter(ei.m_neighbors));

        m_element_info.push_back(ei);
      }
      */




      /*
      void add_vertex(
          vertex_number vn, 
          const hedge::vector &pos,
          python::object adjacent_elements)
      {
        if (vn != m_vertex_adj_elements.size())
          throw std::runtime_error("vertices must be added in increasing order");

        m_vertices.push_back(pos);

        std::auto_ptr<el_id_vector> adj_els_cpp(new el_id_vector);
        std::copy(
            python::stl_input_iterator<element_number>(adjacent_elements),
            python::stl_input_iterator<element_number>(),
            back_inserter(*adj_els_cpp));

        m_vertex_adj_elements.push_back(adj_els_cpp);
      }
      */




      /*
      void add_nodes(unsigned sizehint, python::object iterable)
      {
        m_nodes.reserve(sizehint);
        python::stl_input_iterator<const hedge::vector &> 
          first(iterable), last;
        std::copy(first, last, std::back_inserter(m_nodes));
      }
      */




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
