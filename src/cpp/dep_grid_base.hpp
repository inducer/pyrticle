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




#ifndef _AKFAHFHASHDFH_PYRTICLE_DEP_GRID_BASE_HPP_INCLUDED
#define _AKFAHFHASHDFH_PYRTICLE_DEP_GRID_BASE_HPP_INCLUDED




#include "bases.hpp"
#include "meshdata.hpp"
#include "tools.hpp"
#include "grid.hpp"




namespace pyrticle
{
  struct grid_targets
  {
    template <class VecType>
    class rho_target
    {
      private:
        typename VecType::iterator m_target;

      public:
        rho_target(VecType &target_vector)
          : m_target(target_vector.begin())
        { target_vector.clear(); }

        void begin_particle(particle_number pn)
        { }

        void add_shape_value(unsigned vec_idx, double q_shapeval)
        { m_target[vec_idx] += q_shapeval; }

        void end_particle(particle_number pn)
        { }
    };




    /** Reconstruction Target for the current density.
     */
    template<unsigned DimensionsVelocity, class JVecType, class VVecType>
    class j_target
    {
      private:
        typename JVecType::iterator m_target;
        typename VVecType::const_iterator m_velocities;
        double m_scale_factors[DimensionsVelocity];

      public:
        j_target(JVecType &target_vector, const VVecType &velocities)
          : m_target(target_vector.begin()), m_velocities(velocities.begin())
        { 
          target_vector.clear();
          for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
            m_scale_factors[axis] = 0;
        }

        void begin_particle(particle_number pn)
        {
          for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
            m_scale_factors[axis] = m_velocities[pn*DimensionsVelocity+axis];
        }

        void add_shape_value(const unsigned vec_idx, const double q_shapeval)
        { 
          unsigned const base = vec_idx*DimensionsVelocity;
          for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
            m_target[base+axis] += m_scale_factors[axis]*q_shapeval;
        }

        void end_particle(particle_number pn)
        { }
    };





    template <class T1, class T2>
    class chained_target
    {
      private:
        T1 m_target1;
        T2 m_target2;

      public:
        chained_target(T1 &target1, T2 &target2)
          : m_target1(target1), m_target2(target2)
        { }

        void begin_particle(const particle_number pn)
        {
          m_target1.begin_particle(pn);
          m_target2.begin_particle(pn);
        }

        void add_shape_value(unsigned vec_idx, double q_shapeval)
        { 
          m_target1.add_shape_value(vec_idx, q_shapeval);
          m_target2.add_shape_value(vec_idx, q_shapeval);
        }

        void end_particle(const particle_number pn)
        {
          m_target1.end_particle(pn);
          m_target2.end_particle(pn);
        }
    };




    template <class T1, class T2>
    static 
    chained_target<T1, T2> 
    make_chained_target(T1 &target1, T2 &target2)
    {
      return chained_target<T1, T2>(target1, target2);
    }
  };




  struct grid_depositor_base_state 
  {
    boost::numeric::ublas::vector<brick_number> m_particle_brick_numbers;

    void note_move(particle_number from, particle_number to, unsigned size)
    {
      for (unsigned i = 0; i < size; ++i)
        m_particle_brick_numbers[to+i] = m_particle_brick_numbers[from+i];
    }

    void note_change_size(unsigned particle_count)
    {
      unsigned prev_count = m_particle_brick_numbers.size();
      m_particle_brick_numbers.resize(particle_count);

      if (particle_count > prev_count)
        std::fill(
            m_particle_brick_numbers.begin() + prev_count,
            m_particle_brick_numbers.end(),
            0);
    }
  };




  template <class Derived, class ParticleState, class ShapeFunction, class Brick,
           class DepositorState>
  class grid_depositor_base
  {
    public:
      typedef ParticleState particle_state;
      typedef ShapeFunction shape_function;
      typedef Brick brick_type;
      typedef DepositorState depositor_state;


      const mesh_data &m_mesh_data;
      std::vector<brick_type> m_bricks;
      shape_function m_shape_function;




      grid_depositor_base(const mesh_data &md)
        : m_mesh_data(md)
      { }



      unsigned grid_node_count() const
      {
        if (this->m_bricks.size() == 0)
          return 0;
        else
          return 
            this->m_bricks.back().start_index() 
            + this->m_bricks.back().node_count();
      }




      template <class Target>
      void deposit_single_particle_with_cache(
          depositor_state &ds,
          const particle_state &ps,
          Target tgt,
          particle_number pn,
          bounded_vector const &center,
          bounded_box const &particle_box)
      {
        brick_number &bn_cache(ds.m_particle_brick_numbers[pn]);
        const brick_type &last_brick = m_bricks[bn_cache];

        bool is_complete;
        bool does_intersect_cached = 
          static_cast<const Derived *>(this)->deposit_particle_on_one_brick(
              tgt, last_brick, center, particle_box, ps.charges[pn],
              &is_complete);

        if (!is_complete)
        {
          BOOST_FOREACH(brick_type const &brk, m_bricks)
          {
            // don't re-target the cached brick
            if (brk.number() == last_brick.number())
              continue;

            if (static_cast<const Derived *>(this)->deposit_particle_on_one_brick(
                  tgt, brk, center, particle_box, ps.charges[pn])
                && !does_intersect_cached)
            {
              // We did not intersect the cached brick, but we
              // found a brick that we *did* intersect with.
              // Update the cache.
              bn_cache = brk.number();
            }
          }
        }
      }




      template <class Target>
      void deposit_single_particle_without_cache(Target tgt,
          bounded_vector const &center,
          bounded_box const &particle_box,
          const double charge)
      {
        BOOST_FOREACH(brick_type const &brk, m_bricks)
          static_cast<const Derived *>(this)->deposit_particle_on_one_brick(
              tgt, brk, center, particle_box, charge);
      }





      /** This is a set of bit fields. Each member records the fact
       * that a particular combination of periodic transitions has
       * been taken care of.
       *
       * The assumption is that a particle will only be big enough
       * to reach across the boundary normal to each axis exactly
       * once (regardless of whether that's the +X or -X boundary,
       * it is assumed to only touch one of them).
       *
       * Therefore, a combination of periodic transitions may be 
       * represented by a bit field where each bit corresponds to 
       * one axis, and that is exactly what this data type is.
       */
      typedef unsigned axis_bitfield;
      typedef std::vector<axis_bitfield> periodicity_set;

      template <class Target>
      void deposit_periodic_copies(
          depositor_state &ds,
          const particle_state &ps,
          Target tgt, 
          particle_number pn,
          bounded_vector const &center,
          bounded_box const &particle_box,
          axis_bitfield abf,
          periodicity_set &pset)
      {
        using boost::numeric::ublas::zero_vector;

        const unsigned dim_m = m_mesh_data.m_dimensions;

        if (abf == 0)
          deposit_single_particle_with_cache(
              ds, ps, tgt, pn, center, particle_box);
        else
          deposit_single_particle_without_cache(
              tgt, center, particle_box, ps.charges[pn]);
        pset.push_back(abf);

        for (unsigned axis = 0; axis < dim_m; ++axis)
        {
          const mesh_data::periodicity_axis &p_axis(
             m_mesh_data.m_periodicities[axis]);

          if (p_axis.m_min == p_axis.m_max)
            continue;

          const axis_bitfield abf_plus_this = abf | (1<<axis);
          if (std::find(pset.begin(), pset.end(), abf_plus_this) != pset.end())
            continue;

          if (particle_box.m_lower[axis] < p_axis.m_min)
          {
            bounded_vector per_offset = zero_vector<double>(dim_m);
            per_offset[axis] = p_axis.m_max-p_axis.m_min;
            deposit_periodic_copies(ds, ps, tgt,
                pn, center+per_offset,
                bounded_box(
                  particle_box.m_lower+per_offset,
                  particle_box.m_upper+per_offset),
                abf_plus_this, pset);
          }
          else if (particle_box.m_upper[axis] > p_axis.m_max)
          {
            bounded_vector per_offset = zero_vector<double>(dim_m);
            per_offset[axis] = -(p_axis.m_max-p_axis.m_min);
            deposit_periodic_copies(ds, ps, tgt,
                pn, center+per_offset,
                bounded_box(
                  particle_box.m_lower+per_offset,
                  particle_box.m_upper+per_offset),
                abf_plus_this, pset);
          }
        }
      }





      template <class Target>
      void deposit_densities_on_grid_target(
          depositor_state &ds, const particle_state &ps,
          Target tgt, boost::python::slice const &pslice)
      {
        const unsigned dim_x = ps.xdim();
        const unsigned dim_m = m_mesh_data.m_dimensions;

        using boost::numeric::ublas::scalar_vector;
        const scalar_vector<double> shape_extent(
            dim_m, m_shape_function.radius());

        FOR_ALL_SLICE_INDICES(particle_number, pn, 
            pslice, ps.particle_count)
        {
          tgt.begin_particle(pn);
          const bounded_vector center = subrange(
              ps.positions, pn*dim_x, (pn+1)*dim_x);

          bounded_box particle_box(
              center - shape_extent, 
              center + shape_extent);

          periodicity_set pset;
          deposit_periodic_copies(ds, ps, tgt, pn, center, particle_box, 0, pset);
          tgt.end_particle(pn);
        }
      }
  };
}




#endif
