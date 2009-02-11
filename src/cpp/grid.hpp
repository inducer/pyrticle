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




#ifndef _AHTRNAGHDA_PYRTICLE_GRID_HPP_INCLUDED
#define _AHTRNAGHDA_PYRTICLE_GRID_HPP_INCLUDED




#include "bases.hpp"
#include "meshdata.hpp"
#include "tools.hpp"




namespace pyrticle
{
  typedef npy_uint grid_node_number;
  typedef npy_uint16 brick_number;
  typedef npy_uint brick_node_number;




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




  template<class Brick>
  class brick_iterator
  {
    protected:
      const Brick &m_brick;
      const bounded_int_box m_bounds;

      /* these fields must be updated in lock-step */
      bounded_int_vector m_state;
      bounded_vector m_point;
      grid_node_number m_index;

    public:
      brick_iterator(Brick const &brk, const bounded_int_box &bounds)
        : m_brick(brk), m_bounds(bounds), m_state(bounds.m_lower), 
        m_point(brk.point(m_state)),
        m_index(brk.index(m_state))
      { 
        if (m_bounds.is_empty())
          m_state = m_bounds.m_upper;
      }

      virtual ~brick_iterator()
      { }

      const bounded_int_vector &operator*() const
      { return m_state; }

    protected:
      void inc_bottom_half(unsigned i)
      {
        // getting here means that an overflow occurred, and we'll have to
        // update both m_index and m_point from scratch (unless somebody
        // figures out the reset math for those, especially m_index).

        ++i;
        while (i < m_state.size())
        {
          ++m_state[i];

          if (m_state[i] < m_bounds.m_upper[i])
          {
            m_point = m_brick.point(m_state);
            m_index = m_brick.index(m_state);
            return;
          }

          m_state[i] = m_bounds.m_lower[i];
          ++i;
        }

        // we overflowed everything back to the origin, meaning we're done.
        // no need to update m_point and m_index.

        m_state = m_bounds.m_upper;
        return;
      }


    public:
      virtual brick_iterator &operator++()
      {
        const unsigned i = 0; // silently assuming non-zero size here

        ++m_state[i];

        if (m_state[i] >= m_bounds.m_upper[i])
        {
          m_state[i] = m_bounds.m_lower[i];

          // split off the non-default case bottom half of this routine in the
          // hope that at least the fast default case will get inlined.
          inc_bottom_half(i); 
        }
        else
        {
          m_point[i] += m_brick.stepwidths()[i];
          m_index += m_brick.strides()[i];
        }

        return *this;
      }

      bool at_end() const
      { return m_state[0] >= m_bounds.m_upper[0]; }

      grid_node_number index() const
      { return m_index; }

      const bounded_vector &point() const
      { return m_point; }
  };




  /** A brick defines a cell-centered grid on a certain, box-like region 
   * of space.
   * 
   * .--------------------------------.
   * | * | * | ...                | * |
   * |---+---+--------------------+---|
   * | * | * | ...                | * |
   * `--------------------------------'
   *
   * The thin line above denotes the area said to be "covered" by the 
   * brick (mostly for the assignment of extra points, see below).
   * Each "*" symbol denotes one grid point.
   *
   * m_stepwidths gives the distance along each axis between gridpoints.
   * m_dimensions gives the number of cells/grid points in each direction.
   * m_origin gives the lower left hand corner of the box.
   * m_origin_plus_half = m_origin + m_stepwidths/2, useful for 
   * easy computation of the points.
   */
  struct brick
  {
    protected:
      brick_number m_number;
      grid_node_number m_start_index;
      bounded_vector m_stepwidths;
      bounded_vector m_origin;
      bounded_vector m_origin_plus_half;

      /** This is the number of points in each dimension. If it is 2, then the
       * indices 0 and 1 exist in that dimension.
       */
      bounded_int_vector m_dimensions;
      bounded_int_vector m_strides;

    public:
      brick(
          brick_number number,
          grid_node_number start_index,
          bounded_vector stepwidths,
          bounded_vector origin,
          bounded_int_vector dimensions)
        : m_number(number), m_start_index(start_index), m_dimensions(dimensions)
      {
        unsigned d = m_dimensions.size();

        for (unsigned i = 0; i < d; ++i)
          if (stepwidths[i] < 0)
          {
            stepwidths[i] *= -1;
            origin[i] -= stepwidths[i]*dimensions[i];
          }

        m_stepwidths = stepwidths;
        m_origin = origin;
        m_origin_plus_half = origin+stepwidths/2;

        // This ordering is what Visit expects by default and calls
        // "row-major" (I suppose Y-major).

        m_strides.resize(d);
        unsigned i = 0;
        unsigned current_stride = 1;
        while (i < m_dimensions.size())
        {
          m_strides[i] = current_stride;
          current_stride *= m_dimensions[i];
          ++i;
        }
      }

      virtual ~brick() 
      { }

      brick_number number() const
      { return m_number; }

      grid_node_number start_index() const
      { return m_start_index; }

      bounded_vector const &stepwidths() const
      { return m_stepwidths; }

      bounded_vector const &origin() const
      { return m_origin; }

      bounded_int_vector const &dimensions() const
      { return m_dimensions; }

      bounded_int_vector const &strides() const
      { return m_strides; }

      unsigned node_count() const
      {
        unsigned result = 1;
        BOOST_FOREACH(unsigned d, m_dimensions)
          result *= d;
        return result;
      }

      double cell_volume() const
      {
        double result = 1;
        BOOST_FOREACH(double dx, m_stepwidths)
          result *= dx;
        return result;
      }

      virtual bounded_vector point(const bounded_int_vector &idx) const
      { return m_origin_plus_half + element_prod(idx, m_stepwidths); }

      grid_node_number index(const bounded_int_vector &idx) const
      { return m_start_index + inner_prod(idx, m_strides); }

      bounded_int_vector split_index(grid_node_number gnn) const
      {
        bounded_int_vector result(m_dimensions.size());

        int i = m_dimensions.size();
        while (i >= 0)
        {
          result[i] = gnn/m_strides[i];
          gnn -= result[i]*m_strides[i];
          --i;
        }

        return result;
      }

      bounded_int_vector which_cell(const bounded_vector &pt) const
      {
        bounded_int_vector result = pyublas::unary_op<int_floor>::apply(
              element_div(pt-m_origin, m_stepwidths));
        
        for (unsigned i = 0; i < result.size(); ++i)
          if (result[i] < 0 || result[i] >= m_dimensions[i])
            throw std::invalid_argument("point is out of this brick's bounds");

        return result;
      }

      bounded_box bounding_box() const
      {
        return bounded_box(
            m_origin,
            m_origin + element_prod(m_dimensions, m_stepwidths)
            );
      }

      bounded_box cell (const bounded_int_vector &idx) const
      {
        return bounded_box(
            m_origin + element_prod(idx, m_stepwidths),
            m_origin + element_prod(idx, m_stepwidths) + m_stepwidths);
      }

    private:
      struct int_ceil
      {
        typedef double value_type;
        typedef const double &argument_type;
        typedef int result_type;

        static result_type apply(argument_type x)
        {
          return int(ceil(x));
        }
      };

      struct int_floor
      {
        typedef double value_type;
        typedef const double &argument_type;
        typedef int result_type;

        static result_type apply(argument_type x)
        {
          return int(floor(x));
        }
      };

    public:
      bounded_int_box index_range(const bounded_box &bbox) const
      {
        return bounded_int_box(
            pyublas::unary_op<int_floor>::apply(
              element_div(bbox.m_lower-m_origin, m_stepwidths)),
            pyublas::unary_op<int_ceil>::apply(
              element_div(bbox.m_upper-m_origin, m_stepwidths))
            );
      }

      typedef brick_iterator<brick> iterator;
      
      iterator get_iterator(bounded_int_box const &bounds) const
      { return iterator(*this, bounds); }

      iterator get_iterator(bounded_box const &bounds) const
      { return iterator(*this, index_range(bounds)); }
  };




  template <class PICAlgorithm, class Brick>
  class grid_reconstructor_base : public reconstructor_base
  {
    public:
      // member data --------------------------------------------------------
      shape_function   m_shape_function;

      std::vector<Brick> m_bricks;
      boost::numeric::ublas::vector<brick_number> m_particle_brick_numbers;





      // internals ------------------------------------------------------------
      unsigned grid_node_count() const
      {
        if (m_bricks.size() == 0)
          return 0;
        else
          return m_bricks.back().start_index() + m_bricks.back().node_count();
      }




      template <class Target>
      void reconstruct_single_particle_with_cache(Target tgt,
          particle_number pn,
          bounded_vector const &center,
          bounded_box const &particle_box)
      {
        const double charge = CONST_PIC_THIS->m_charges[pn];
        brick_number &bn_cache(m_particle_brick_numbers[pn]);
        const Brick &last_brick = m_bricks[bn_cache];

        bool is_complete;
        bool does_intersect_cached = 
          PIC_THIS->reconstruct_particle_on_one_brick(
              tgt, last_brick, center, particle_box, charge,
              &is_complete);

        if (!is_complete)
        {
          BOOST_FOREACH(Brick const &brk, m_bricks)
          {
            // don't re-target the cached brick
            if (brk.number() == last_brick.number())
              continue;

            if (PIC_THIS->reconstruct_particle_on_one_brick(
                  tgt, brk, center, particle_box, charge)
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
      void reconstruct_single_particle_without_cache(Target tgt,
          particle_number pn,
          bounded_vector const &center,
          bounded_box const &particle_box)
      {
        BOOST_FOREACH(Brick const &brk, m_bricks)
          PIC_THIS->reconstruct_particle_on_one_brick(
              tgt, brk, center, particle_box, 
              CONST_PIC_THIS->m_charges[pn]);
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
      void reconstruct_periodic_copies(Target tgt, 
          particle_number pn,
          bounded_vector const &center,
          bounded_box const &particle_box,
          axis_bitfield abf,
          periodicity_set &pset)
      {
        using boost::numeric::ublas::zero_vector;

        const unsigned dim_m = CONST_PIC_THIS->m_mesh_data.m_dimensions;

        if (abf == 0)
          reconstruct_single_particle_with_cache(
              tgt, pn, center, particle_box);
        else
          reconstruct_single_particle_without_cache(
              tgt, pn, center, particle_box);
        pset.push_back(abf);

        for (unsigned axis = 0; axis < dim_m; ++axis)
        {
          mesh_data::periodicity_axis &p_axis(
             CONST_PIC_THIS->m_mesh_data.m_periodicities[axis]);

          if (p_axis.m_min == p_axis.m_max)
            continue;

          const axis_bitfield abf_plus_this = abf | (1<<axis);
          if (std::find(pset.begin(), pset.end(), abf_plus_this) != pset.end())
            continue;

          if (particle_box.m_lower[axis] < p_axis.m_min)
          {
            bounded_vector per_offset = zero_vector<double>(dim_m);
            per_offset[axis] = p_axis.m_max-p_axis.m_min;
            reconstruct_periodic_copies(tgt,
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
            reconstruct_periodic_copies(tgt,
                pn, center+per_offset,
                bounded_box(
                  particle_box.m_lower+per_offset,
                  particle_box.m_upper+per_offset),
                abf_plus_this, pset);
          }
        }
      }




      template <class Target>
      void reconstruct_densities_on_grid_target(Target tgt,
            boost::python::slice const &pslice)
      {
        const unsigned dim_x = CONST_PIC_THIS->get_dimensions_pos();
        const unsigned dim_m = CONST_PIC_THIS->m_mesh_data.m_dimensions;

        using boost::numeric::ublas::scalar_vector;
        const scalar_vector<double> shape_extent(
            dim_m, m_shape_function.radius());

        FOR_ALL_SLICE_INDICES(particle_number, pn, 
            pslice, CONST_PIC_THIS->m_particle_count)
        {
          tgt.begin_particle(pn);
          const bounded_vector center = subrange(
              CONST_PIC_THIS->m_positions, pn*dim_x, (pn+1)*dim_x);

          bounded_box particle_box(
              center - shape_extent, 
              center + shape_extent);

          periodicity_set pset;
          reconstruct_periodic_copies(tgt, pn, center, particle_box, 0, pset);
          tgt.end_particle(pn);
        }
      }




      // particle numbering notifications -----------------------------------
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
}




#endif
