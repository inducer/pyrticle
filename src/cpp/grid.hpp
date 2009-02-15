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
#include <boost/random.hpp>




namespace pyrticle
{
  typedef npy_uint grid_node_number;
  typedef npy_uint16 brick_number;
  typedef npy_uint brick_node_number;




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
          if (m_dimensions[i] == 0)
            throw std::logic_error("zero dimensions for brick not supported");
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

        int i = m_dimensions.size()-1;
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




  struct jiggly_brick : public brick
  {
    private:
      static const unsigned point_origin_count = 3;
      typedef brick super;
      std::vector<bounded_vector> m_point_origins;
      std::vector<bounded_vector> m_axis0_offsets;

    public:
      jiggly_brick(
          brick_number number,
          grid_node_number start_index,
          bounded_vector stepwidths,
          bounded_vector origin,
          bounded_int_vector dimensions,
          double jiggle_radius)
        : super(number, start_index, stepwidths, origin, dimensions)
      { 
        boost::variate_generator<
          boost::mt19937, 
          boost::uniform_real<> > offset_rng(
              boost::mt19937(), boost::uniform_real<>(-jiggle_radius, jiggle_radius));

        for (unsigned i = 0; i < point_origin_count; ++i)
        {
          bounded_vector offset(m_origin.size());
          for (unsigned j = 0; j < m_origin.size(); ++j)
            offset[j] = offset_rng() * m_stepwidths[j];
          m_point_origins.push_back(m_origin_plus_half + offset);
        }

        for (unsigned i = 0; i < point_origin_count; ++i)
        {
          bounded_vector axis0_offset(m_origin.size());
          axis0_offset = m_point_origins[
            (i+1)%point_origin_count] - m_point_origins[i];
          axis0_offset[0] += m_stepwidths[0];
          m_axis0_offsets.push_back(axis0_offset);
        }
      }

      bounded_vector point(const bounded_int_vector &idx) const
      { 
        return m_point_origins[origin_index(idx)]
          + element_prod(idx, m_stepwidths); 
      }

      unsigned origin_index(const bounded_int_vector &idx) const
      { 
        return std::accumulate(idx.begin(), idx.end(), 0) 
          % point_origin_count;
      }

      class iterator : public brick_iterator<jiggly_brick>
      {
        private:
          typedef brick_iterator<jiggly_brick> super;

        protected:
          unsigned m_origin_index;

        public:
          iterator(jiggly_brick const &brk, 
              const bounded_int_box &bounds)
            : super(brk, bounds)
          { 
            update_origin_index();
          }

        protected:
          void update_origin_index()
          { m_origin_index = m_brick.origin_index(m_state); }

        public:
          iterator &operator++()
          {
            ++m_state[0];

            if (m_state[0] < m_bounds.m_upper[0])
            {
              m_origin_index = (m_origin_index + 1) % point_origin_count;
              m_point += m_brick.m_axis0_offsets[m_origin_index];
              m_index += m_brick.strides()[0];
            }
            else
            {
              m_state[0] = m_bounds.m_lower[0];

              // split off the non-default case bottom half of this routine in the
              // hope that at least the fast default case will get inlined.
              inc_bottom_half(0); 
              update_origin_index();
            }
            return *this;
          }
      };

      friend class iterator;

      iterator get_iterator(bounded_int_box const &bounds) const
      { return iterator(*this, bounds); }

      iterator get_iterator(bounded_box const &bounds) const
      { return iterator(*this, index_range(bounds)); }
  };
}




#endif
