// Pyrticle - Particle in Cell in Python
// Element finding
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





#ifndef _BADFJAH_PYRTICLE_ELEMENT_FINDER_HPP_INCLUDED
#define _BADFJAH_PYRTICLE_ELEMENT_FINDER_HPP_INCLUDED




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
  class heuristic_element_finder
  {
    private:
      const mesh_data &m_mesh_data;

    public:
      heuristic_element_finder(const mesh_data &md)
        : m_mesh_data(md)
      { }

    private:
      template <class ElementTarget>
      void add_shape_by_neighbors(
          ElementTarget &target,
          const bounded_vector &pos,
          const mesh_data::element_info &einfo,
          double radius)
      {
        target.add_shape_on_element(pos, einfo.m_id);

        BOOST_FOREACH(const mesh_data::face_info &face, einfo.m_faces)
        {
          if (face.m_neighbor != mesh_data::INVALID_ELEMENT)
          {
            const mesh_data::axis_number per_axis = face.m_neighbor_periodicity_axis;

            if (per_axis == mesh_data::INVALID_AXIS)
              target.add_shape_on_element(pos, face.m_neighbor);
            else
            {
              bounded_vector pos2(pos);
              const mesh_data::periodicity_axis &pa = 
                m_mesh_data.m_periodicities[per_axis];

              if (pos[per_axis] - radius < pa.m_min)
              {
                pos2[per_axis] += (pa.m_max-pa.m_min);
                target.add_shape_on_element(pos2, face.m_neighbor);
              }
              if (pos[per_axis] + radius > pa.m_max)
              {
                pos2[per_axis] -= (pa.m_max-pa.m_min);
                target.add_shape_on_element(pos2, face.m_neighbor);
              }
            }
          }
        }
      }




      template <class ElementTarget>
      void add_shape_by_vertex(
          ElementTarget &target,
          const bounded_vector &pos,
          const mesh_data::element_info &einfo,
          double radius)
      {
        // find closest vertex
        mesh_data::vertex_number closest_vertex = 
          mesh_data::INVALID_VERTEX;
        double min_dist = std::numeric_limits<double>::infinity();

        BOOST_FOREACH(mesh_data::vertex_number vi, einfo.m_vertices)
        {
          double dist = norm_2(m_mesh_data.mesh_vertex(vi) - pos);
          if (dist < min_dist)
          {
            closest_vertex = vi;
            min_dist = dist;
          }
        }

        const mesh_data &md = m_mesh_data;

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

            bounded_vector pos2(pos);
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




    public:
      template <class ParticleState, class ElementTarget>
      void operator()(
          const ParticleState &ps, 
          ElementTarget &target, 
          particle_number pn, double radius)
      {
        const unsigned dim = ps.xdim();
        const bounded_vector pos = subrange(
            ps.positions, pn*dim, (pn+1)*dim);
        const mesh_data::element_number containing_el = 
          ps.containing_elements[pn];
        const mesh_data::element_info &einfo(
            m_mesh_data.m_element_info[containing_el]);

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




  class face_based_element_finder
  {
    private:
      const mesh_data &m_mesh_data;

    public:
      face_based_element_finder(const mesh_data &md)
        : m_mesh_data(md)
      { }

    private:
      typedef boost::unordered_set<mesh_data::element_number> el_set_t;

      template <class ElementTarget>
      void recurse(ElementTarget &target, const bounded_vector &pos, double radius,
          el_set_t &el_set, mesh_data::element_number en)
      {
        target.add_shape_on_element(pos, en);
        el_set.insert(en);

        const mesh_data::element_info &einfo(
            m_mesh_data.m_element_info[en]);
        BOOST_FOREACH(const mesh_data::face_info &face, einfo.m_faces)
        {
          if (face.m_neighbor == mesh_data::INVALID_ELEMENT)
            continue;

          if (el_set.find(face.m_neighbor) != el_set.end())
            continue;

          // test 1: necessary for inclusion
          // d(blob, faceplane) < particle_radius?
          double plane_dist = fabs(
              inner_prod(pos, face.m_normal) - face.m_face_plane_eqn_rhs);

          if (plane_dist > radius)
            continue;

          // test 2: sufficient for exclusion
          // d(face_bound_circle, pos) < radius?
          if (norm_2(face.m_face_centroid-pos) 
              > radius+face.m_face_radius_from_centroid)
            continue;

          // treat periodicity
          const mesh_data::axis_number per_axis = face.m_neighbor_periodicity_axis;

          if (per_axis == mesh_data::INVALID_AXIS)
            recurse(target, pos, radius, el_set, face.m_neighbor);
          else
          {
            bounded_vector pos2(pos);
            const mesh_data::periodicity_axis &pa = 
              m_mesh_data.m_periodicities[per_axis];

            if (pos[per_axis] - radius < pa.m_min)
            {
              pos2[per_axis] += (pa.m_max-pa.m_min);
              recurse(target, pos2, radius, el_set, face.m_neighbor);
            }
            if (pos[per_axis] + radius > pa.m_max)
            {
              pos2[per_axis] -= (pa.m_max-pa.m_min);
              recurse(target, pos2, radius, el_set, face.m_neighbor);
            }
          }
        }
      }




    public:
      template <class ParticleState, class ElementTarget>
      void operator()(
          const ParticleState &ps, 
          ElementTarget &target, 
          particle_number pn, double radius)
      {
        const unsigned dim = ps.xdim();
        const bounded_vector pos = subrange(
            ps.m_positions, pn*dim, (pn+1)*dim);

        el_set_t el_set;
        recurse(target, pos, radius, el_set, 
            ps.m_containing_elements[pn]);
      }
  };




  typedef face_based_element_finder element_finder;
}





#endif
