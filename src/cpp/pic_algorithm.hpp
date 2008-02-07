// Pyrticle - Particle in Cell in Python
// Composition of PIC algorithms from reconstructors and
//   particle pushers
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





#ifndef _BADFJAH_PYRTICLE_PIC_ALGORITHM_HPP_INCLUDED
#define _BADFJAH_PYRTICLE_PIC_ALGORITHM_HPP_INCLUDED




#include <boost/foreach.hpp> 
#include <boost/shared_ptr.hpp> 
#include <boost/numeric/ublas/vector_proxy.hpp>
#include "meshdata.hpp"




namespace pyrticle 
{
  template<unsigned DimensionsPos, unsigned DimensionsVelocity>
  struct pic_data
  {
    template <class PICAlgorithm>
    class type : boost::noncopyable
    {
      public:
        mesh_data                         m_mesh_data;

        unsigned                          m_particle_count;

        mesh_data::el_id_vector           m_containing_elements;
        hedge::vector                     m_positions;
        hedge::vector                     m_momenta;
        hedge::vector                     m_charges;
        hedge::vector                     m_masses;

        const double                      m_vacuum_c;

        mutable event_counter             m_find_same;
        mutable event_counter             m_find_by_neighbor;
        mutable event_counter             m_find_by_vertex;
        mutable event_counter             m_find_global;

        boost::shared_ptr<number_shift_listener> m_particle_number_shift_listener;

      public:
        // administrative stuff -----------------------------------------------
        type(unsigned mesh_dimensions, double vacuum_c)
          : m_mesh_data(mesh_dimensions), 
          m_particle_count(0),
          m_vacuum_c(vacuum_c)
        { }

        // dimensions ---------------------------------------------------------
        /* CAVEAT: Use the get_dimensions_XXX functions whenever possible.
         * Otherwise, you will get spurious linker errors about undefined
         * references to these constants.
         *
         * In particular, these expressions should only be used where they're
         * evaluated at compile time.
         */
        static const unsigned dimensions_pos = DimensionsPos;
        static const unsigned dimensions_velocity = DimensionsVelocity;

        static unsigned get_dimensions_pos()
        { return DimensionsPos; }

        static unsigned get_dimensions_velocity()
        { return DimensionsVelocity; }

        // heavy-lifting routines ---------------------------------------------
        const hedge::vector velocities() const
        {
          const unsigned vdim = DimensionsVelocity;

          hedge::vector result(m_particle_count * DimensionsVelocity);

          for (particle_number pn = 0; pn < m_particle_count; pn++)
          {
            unsigned vpstart = vdim*pn;
            unsigned vpend = vdim*(pn+1);

            const double m = m_masses[pn];
            double p = norm_2(subrange(m_momenta, vpstart, vpend));
            double v = m_vacuum_c*p/sqrt(m*m*m_vacuum_c*m_vacuum_c + p*p);
            subrange(result, vpstart, vpend) = v/p*subrange(m_momenta, vpstart, vpend);
          }
          return result;
        }




        void add_rhs(
            hedge::vector const &dx, 
            hedge::vector const &dp)
        {
          noalias(subrange(m_positions, 0, DimensionsPos*m_particle_count))
            += dx;
          noalias(subrange(m_momenta, 0, DimensionsVelocity*m_particle_count))
            += dp;
        }




        mesh_data::element_number find_new_containing_element(particle_number i,
            mesh_data::element_number prev) const
        {
          const unsigned xdim = DimensionsPos;

          const unsigned x_pstart = i*xdim;
          const unsigned x_pend = (i+1)*xdim;
          
          const hedge::vector pt = subrange(m_positions, x_pstart, x_pend);

          if (prev != mesh_data::INVALID_ELEMENT)
          {
            const mesh_data::element_info &prev_el = m_mesh_data.m_element_info[prev];

            // check if we're still in the same element -------------------------
            if (is_in_unit_simplex(prev_el.m_inverse_map(pt)))
            {
              m_find_same.tick();
              return prev;
            }

            // we're not: lookup via normal -------------------------------------
            {
              int closest_normal_idx = -1;
              double max_ip = 0;
              unsigned normal_idx = 0;

              BOOST_FOREACH(const hedge::vector &n, prev_el.m_normals)
              {
                double ip = inner_prod(n, subrange(m_momenta, x_pstart, x_pend));
                if (ip > max_ip)
                {
                  closest_normal_idx = normal_idx;
                  max_ip = ip;
                }
                ++normal_idx;
              }

              if (closest_normal_idx == -1)
                throw std::runtime_error("no best normal found--weird");

              mesh_data::element_number possible_idx =
                prev_el.m_neighbors[closest_normal_idx];

              if (possible_idx != mesh_data::INVALID_ELEMENT)
              {
                const mesh_data::element_info &possible = 
                  m_mesh_data.m_element_info[possible_idx];

                if (is_in_unit_simplex(possible.m_inverse_map(pt)))
                {
                  m_find_by_neighbor.tick();
                  return possible.m_id;
                }
              }
            }

            // look up via closest vertex ---------------------------------------
            {
              mesh_data::vertex_number closest_vertex = 
                mesh_data::INVALID_VERTEX;

              {
                double min_dist = std::numeric_limits<double>::infinity();

                BOOST_FOREACH(mesh_data::vertex_number vi, prev_el.m_vertices)
                {
                  double dist = norm_2(m_mesh_data.m_vertices[vi] - pt);
                  if (dist < min_dist)
                  {
                    closest_vertex = vi;
                    min_dist = dist;
                  }
                }
              }

              // found closest vertex, go through adjacent elements
              BOOST_FOREACH(mesh_data::element_number possible_idx, 
                  std::make_pair(
                    m_mesh_data.m_vertex_adj_elements.begin()
                    + m_mesh_data.m_vertex_adj_element_starts[closest_vertex],
                    m_mesh_data.m_vertex_adj_elements.begin()
                    + m_mesh_data.m_vertex_adj_element_starts[closest_vertex+1]
                    )
                  )
              {
                const mesh_data::element_info &possible = 
                  m_mesh_data.m_element_info[possible_idx];

                if (is_in_unit_simplex(possible.m_inverse_map(pt)))
                {
                  m_find_by_vertex.tick();
                  return possible.m_id;
                }
              }

            }

          }

          // last resort: global search ---------------------------------------
          {
            m_find_global.tick();

            mesh_data::element_number new_el = 
              m_mesh_data.find_containing_element(pt);
            if (new_el != mesh_data::INVALID_ELEMENT)
              return new_el;
          }

          return mesh_data::INVALID_ELEMENT;
        }




        void update_containing_elements()
        {
          for (particle_number i = 0; i < m_particle_count; i++)
          {
            mesh_data::element_number prev = m_containing_elements[i];

            mesh_data::element_number new_el = 
              find_new_containing_element(i, prev);

            // no element found? must be boundary
            if (new_el == mesh_data::INVALID_ELEMENT)
              boundary_hit(i);
            else
              m_containing_elements[i] = new_el;
          }
        }





        void boundary_hit(particle_number pn)
        {
          unsigned x_pstart = pn*DimensionsPos;
          unsigned x_pend = (pn+1)*DimensionsPos;

          bool periodicity_trip = false;

          hedge::vector pt = subrange(m_positions, x_pstart, x_pend);
          BOOST_FOREACH(const mesh_data::periodicity_axis &pa, 
              m_mesh_data.m_periodicities)
          {
            const double xi = pt[pa.m_axis];
            if (xi < pa.m_min)
            {
              pt[pa.m_axis] += (pa.m_max-pa.m_min);
              periodicity_trip = true;
            }
            else if (xi > pa.m_max)
            {
              pt[pa.m_axis] -= (pa.m_max-pa.m_min);
              periodicity_trip = true;
            }
          }

          if (periodicity_trip)
          {
            subrange(m_positions, x_pstart, x_pend) = pt;
            mesh_data::element_number ce = 
              find_new_containing_element(pn, m_containing_elements[pn]);
            if (ce != mesh_data::INVALID_ELEMENT)
            {
              m_containing_elements[pn] = ce;
              return;
            }
          }

          kill_particle(pn);
        }




        void kill_particle(particle_number pn)
        {
          std::cout << "KILL " << pn << std::endl;

          --m_particle_count;
          move_particle(m_particle_count, pn);

          if (m_particle_number_shift_listener.get())
            m_particle_number_shift_listener->note_change_size(m_particle_count);
          PIC_THIS->PICAlgorithm::reconstructor::note_change_size(m_particle_count);
          PIC_THIS->PICAlgorithm::particle_pusher::note_change_size(m_particle_count);
        }




        void move_particle(particle_number from, particle_number to)
        {
          const unsigned xdim = DimensionsPos;
          const unsigned vdim = DimensionsVelocity;

          m_containing_elements[to] = m_containing_elements[from];

          if (m_particle_number_shift_listener.get())
            m_particle_number_shift_listener->note_move(from, to, 1);
          PIC_THIS->PICAlgorithm::reconstructor::note_move(from, to, 1);
          PIC_THIS->PICAlgorithm::particle_pusher::note_move(from, to, 1);

          for (unsigned i = 0; i < xdim; i++)
            m_positions[to*xdim+i] = m_positions[from*xdim+i];

          for (unsigned i = 0; i < vdim; i++)
            m_momenta[to*vdim+i] = m_momenta[from*vdim+i];

          m_charges[to] = m_charges[from];
          m_masses[to] = m_masses[from];
        }
    };

  };
  




  template <class PICData, 
           class Reconstructor, 
           class ParticlePusher>
  class pic : 
    public PICData::template type<pic<PICData, Reconstructor, ParticlePusher> >,
    public Reconstructor::template type<pic<PICData, Reconstructor, ParticlePusher> >,
    public ParticlePusher::template type<pic<PICData, Reconstructor, ParticlePusher> >
  {
    private:
      boost::shared_ptr<visualization_listener> m_vis_listener;

    public:
      typedef typename PICData::template type<pic> pic_data;
      typedef typename Reconstructor::template type<pic> reconstructor;
      typedef typename ParticlePusher::template type<pic> particle_pusher;

      pic(unsigned mesh_dimensions, double vacuum_c)
        : pic_data(mesh_dimensions, vacuum_c)
      { }

      // visualization-related ------------------------------------------------
      void store_vis_vector(
          const char *name,
          const hedge::vector &vec)
      {
        if (m_vis_listener.get())
          m_vis_listener->store_vis_vector(name, vec);
      }

      void set_vis_listener(
          boost::shared_ptr<visualization_listener> listener)
      {
        m_vis_listener = listener;
      }
  };
}




#endif
