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




#include "tools.hpp"




namespace pyrticle
{
  typedef unsigned grid_node_number;
  typedef unsigned brick_node_number;




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
   * easy computation
   */
  struct brick
  {
    protected:
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
          grid_node_number start_index,
          bounded_vector stepwidths,
          bounded_vector origin,
          bounded_int_vector dimensions)
        : m_start_index(start_index), m_dimensions(dimensions)
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

      bounded_box bounding_box() const
      {
        return bounded_box(
            m_origin,
            m_origin + element_prod(m_dimensions, m_stepwidths)
            );
      }

    private:
      struct int_round
      {
        typedef double value_type;
        typedef const double &argument_type;
        typedef int result_type;

        static result_type apply(argument_type x)
        {
          return int(round(x));
        }
      };

    public:
      bounded_int_box index_range(const bounded_box &bbox) const
      {
        return bounded_int_box(
            pyublas::unary_op<int_round>::apply(
              element_div(bbox.m_lower-m_origin, m_stepwidths)),
            pyublas::unary_op<int_round>::apply(
              element_div(bbox.m_upper-m_origin, m_stepwidths))
            );
      }

      typedef brick_iterator<brick> iterator;
      
      iterator get_iterator(bounded_int_box const &bounds) const
      { return iterator(*this, bounds); }

      iterator get_iterator(bounded_box const &bounds) const
      { return iterator(*this, index_range(bounds)); }
  };
}




#endif
