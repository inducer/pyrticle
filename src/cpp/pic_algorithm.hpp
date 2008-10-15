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
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include "meshdata.hpp"
#include "rec_target.hpp"




namespace pyrticle 
{
  template <unsigned DimensionsPos, unsigned DimensionsVelocity>
  struct particle_base_state 
  {
    static const unsigned m_xdim = DimensionsPos;
    static const unsigned m_vdim = DimensionsVelocity;

    static unsigned xdim()
    { return DimensionsPos; }

    static unsigned vdim()
    { return DimensionsVelocity; }

    unsigned                          particle_count;

    pyublas::numpy_vector<mesh_data::element_number> containing_elements;
    py_vector                         positions;
    py_vector                         momenta;
    py_vector                         charges;
    py_vector                         masses;

    particle_base_state()
    : particle_count(0)
    { }

    particle_base_state(particle_base_state const &src)
    { assign(src); }

    particle_base_state &operator=(particle_base_state &src)
    { assign(src); }
    
    void assign(particle_base_state const &src)
    {
      particle_count = src.particle_count;
      containing_elements = src.containing_elements.copy();
      positions = src.positions.copy();
      momenta = src.momenta.copy();
      charges = src.charges.copy();
      masses = src.masses.copy();
    }
  };




  class number_shift_listener
  {
    public:
      virtual ~number_shift_listener()
      { }

      virtual void note_change_size(unsigned new_size) const 
      { }
      virtual void note_move(unsigned orig, unsigned dest, unsigned size) const 
      { }
      virtual void note_reset(unsigned start, unsigned size) const 
      { }
  };



  class boundary_hit_listener
  {
    public:
      virtual ~boundary_hit_listener()
      { }

      virtual void note_boundary_hit(particle_number pn) const 
      { }
  };




  struct find_event_counters
  {
    event_counter             find_same;
    event_counter             find_by_neighbor;
    event_counter             find_by_vertex;
    event_counter             find_global;
  };




  // actual functionality -----------------------------------------------------
  template <class ParticleState>
  const py_vector get_velocities(const ParticleState &ps, const double vacuum_c)
  {
    const unsigned vdim = ps.vdim();

    py_vector result(ps.particle_count * vdim);

    for (particle_number pn = 0; pn < ps.particle_count; pn++)
    {
      unsigned vpstart = vdim*pn;
      unsigned vpend = vdim*(pn+1);

      const double m = ps.masses[pn];
      double p = norm_2(subrange(ps.momenta, vpstart, vpend));
      double v = vacuum_c*p/sqrt(m*m*vacuum_c*vacuum_c + p*p);
      subrange(result, vpstart, vpend) = v/p*subrange(ps.momenta, vpstart, vpend);
    }
    return result;
  }




  template <class ParticleState>
  mesh_data::element_number find_new_containing_element(
      const mesh_data &mesh,
      const ParticleState &ps,
      particle_number i,
      mesh_data::element_number prev,
      find_event_counters &counters)
  {
    const unsigned xdim = ps.xdim();

    const unsigned x_pstart = i*xdim;
    const unsigned x_pend = (i+1)*xdim;
    
    const bounded_vector pt = subrange(ps.positions, x_pstart, x_pend);

    if (prev != mesh_data::INVALID_ELEMENT)
    {
      const mesh_data::element_info &prev_el = mesh.element_info[prev];

      // check if we're still in the same element -------------------------
      if (is_in_unit_simplex(prev_el.m_inverse_map(pt)))
      {
        counters.find_same.tick();
        return prev;
      }

      // we're not: lookup via normal -------------------------------------
      {
        int closest_normal_idx = -1;
        double max_ip = 0;
        unsigned normal_idx = 0;

        BOOST_FOREACH(const mesh_data::face_info &f, prev_el.m_faces)
        {
          double ip = inner_prod(
              f.m_normal, 
              subrange(ps.momenta, x_pstart, x_pend));

          if (ip > max_ip)
          {
            closest_normal_idx = normal_idx;
            max_ip = ip;
          }
          ++normal_idx;
        }

        if (closest_normal_idx == -1)
        {
          std::cerr << "face normals:" << std::endl;
          bounded_vector mom = subrange(ps.momenta, x_pstart, x_pend);
          BOOST_FOREACH(const mesh_data::face_info &f, prev_el.m_faces)
          {
            std::cerr 
              << f.m_normal 
              << ", ip with momentum:" << inner_prod(f.m_normal, mom)
              << std::endl;
          }
          throw std::runtime_error(
              str(boost::format("no best normal found--weird (momentum=%1%)") % mom));
        }

        mesh_data::element_number possible_idx =
          prev_el.m_faces[closest_normal_idx].m_neighbor;

        if (possible_idx != mesh_data::INVALID_ELEMENT)
        {
          const mesh_data::element_info &possible = 
            mesh.m_element_info[possible_idx];

          if (is_in_unit_simplex(possible.m_inverse_map(pt)))
          {
            counters.find_by_neighbor.tick();
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
            double dist = norm_2(mesh.mesh_vertex(vi) - pt);
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
              mesh.m_vertex_adj_elements.begin()
              + mesh.m_vertex_adj_element_starts[closest_vertex],
              mesh.m_vertex_adj_elements.begin()
              + mesh.m_vertex_adj_element_starts[closest_vertex+1]
              )
            )
        {
          const mesh_data::element_info &possible = 
            mesh.m_element_info[possible_idx];

          if (is_in_unit_simplex(possible.m_inverse_map(pt)))
          {
            counters.find_by_vertex.tick();
            return possible.m_id;
          }
        }

      }

    }

    // last resort: global search ---------------------------------------
    {
      counters.find_global.tick();

      mesh_data::element_number new_el = 
        mesh.find_containing_element(pt);
      if (new_el != mesh_data::INVALID_ELEMENT)
        return new_el;
    }

    return mesh_data::INVALID_ELEMENT;
  }




  template <class ParticleState>
  void update_containing_elements(
      ParticleState &ps,
      const boundary_hit_listener &bhit_listener
      )
  {
    for (particle_number pn = 0; pn < ps.particle_count;)
    {
      mesh_data::element_number prev = ps.containing_elements[pn];

      mesh_data::element_number new_el = 
        find_new_containing_element(ps, pn, prev);

      // no element found? must be boundary
      if (new_el == mesh_data::INVALID_ELEMENT)
      {
        /* INVARIANT: boundary_hit_listener *must* leave particles
         * with particle number less than pn unchanged.
         */
        boundary_hit(ps, pn, bhit_listener);
      }
      else
      {
        ps.containing_elements[pn] = new_el;
        ++pn;
      }
    }
  }





  template <class ParticleState>
  void boundary_hit(
      const mesh_data &mesh,
      ParticleState &ps, 
      particle_number pn,
      const boundary_hit_listener &bhit_listener)
  {
    unsigned x_pstart = pn*ps.xdim();
    unsigned x_pend = (pn+1)*ps.xdim();

    bool periodicity_trip = false;

    bounded_vector pt = subrange(ps.positions, x_pstart, x_pend);

    mesh_data::axis_number per_axis = 0;
    BOOST_FOREACH(const mesh_data::periodicity_axis &pa, 
        mesh.m_periodicities)
    {
      if (pa.m_min != pa.m_max)
      {
        const double xi = pt[per_axis];
        if (xi < pa.m_min)
        {
          pt[per_axis] += (pa.m_max-pa.m_min);
          periodicity_trip = true;
        }
        else if (xi > pa.m_max)
        {
          pt[per_axis] -= (pa.m_max-pa.m_min);
          periodicity_trip = true;
        }
      }

      ++per_axis;
    }

    if (periodicity_trip)
    {
      subrange(ps.positions, x_pstart, x_pend) = pt;
      mesh_data::element_number ce = 
        find_new_containing_element(pn, ps.containing_elements[pn]);
      if (ce != mesh_data::INVALID_ELEMENT)
      {
        ps.containing_elements[pn] = ce;
        return;
      }
    }

    /* INVARIANT: boundary_hit_listener *must* leave particles
     * with particle number less than pn unchanged.
     */
    bhit_listener.note_boundary_hit(ps, pn);
  }



#if 0
  void kill_particle(
      const particle_dimensionality &pd,
      const particle_state &ps, 
      particle_number pn,
      const number_shift_listener &nshift_listener
      )
  {
    --ps.particle_count;
    if (ps.particle_count)
    {
      // Observe that ps.particle_count was decremented
      // above, so this does indeed move the last valid
      // particle, if it exists.
      move_particle(ps, ps.particle_count, pn, nshift_listener);
    }

    note_change_particle_count(ps.particle_count);
  }




      void note_change_particle_count(particle_state &ps, unsigned particle_count) const
      {
        if (ps.particle_number_shift_listener.get())
          ps.particle_number_shift_listener->note_change_size(ps, ps.particle_count);
        reconstructor::note_change_size(ps, ps.particle_count);
        particle_pusher::note_change_size(ps, ps.particle_count);
      }




      /** Move a particle to a different number.
       *
       * (this has nothing to do with actual particle motion.)
       */
      void move_particle(
          const particle_dimensionality &pd,
          particle_state &ps, 
          particle_number from, particle_number to,
          const number_shift_listener &nshift_listener
          )
      {
        const unsigned xdim = pd.xdim();
        const unsigned vdim = pd.vdim();

        ps.containing_elements[to] = ps.containing_elements[from];

        if (ps.particle_number_shift_listener.get())
          ps.particle_number_shift_listener->note_move(from, to, 1);
        reconstructor::note_move(ps, from, to, 1);
        particle_pusher::note_move(ps, from, to, 1);

        for (unsigned i = 0; i < xdim; i++)
          ps.positions[to*xdim+i] = ps.positions[from*xdim+i];

        for (unsigned i = 0; i < vdim; i++)
          ps.momenta[to*vdim+i] = ps.momenta[from*vdim+i];

        ps.charges[to] = ps.charges[from];
        ps.masses[to] = ps.masses[from];
      }
  };
#endif
}




#endif
