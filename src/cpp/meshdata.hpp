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
      typedef hedge::el_id_vector el_id_vector;
      typedef hedge::vtx_id_vector vtx_id_vector;
      typedef hedge::el_face el_face;

      static const element_number INVALID_ELEMENT = hedge::INVALID_ELEMENT;
      static const vertex_number INVALID_VERTEX = hedge::INVALID_VERTEX;




      // data structures ------------------------------------------------------
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
        double                  m_min, m_max;
      };

      // data members ---------------------------------------------------------
      const unsigned m_dimensions;

      std::vector<element_info> m_element_info;
      std::vector<hedge::vector> m_vertices, m_nodes;

      std::vector<unsigned> m_vertex_adj_element_starts;
      el_id_vector m_vertex_adj_elements;

      std::vector<periodicity_axis> m_periodicities;




      // setup ----------------------------------------------------------------
      mesh_data(unsigned dimensions)
        : m_dimensions(dimensions)
      { }




      static element_number get_INVALID_ELEMENT() { return INVALID_ELEMENT; }




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
