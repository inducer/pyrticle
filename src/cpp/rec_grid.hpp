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




#ifndef _AFAYYTAA_PYRTICLE_REC_GRID_HPP_INCLUDED
#define _AFAYYTAA_PYRTICLE_REC_GRID_HPP_INCLUDED




#include <vector>
#include <limits>
#include <pyublas/numpy.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/lapack/gesvd.hpp>
#include <boost/numeric/bindings/lapack/gesdd.hpp>
#include <boost/numeric/bindings/blas/blas2.hpp>
#include <pyublas/elementwise_op.hpp>
#include "rec_shape.hpp"




namespace pyrticle
{
  struct grid_reconstructor 
  {
    typedef unsigned grid_node_number;
    typedef unsigned brick_number;
    typedef unsigned brick_node_number;




    // specialized targets ----------------------------------------------------
    class rho_target
    {
      private:
        dyn_vector &m_target_vector;

      public:
        rho_target(dyn_vector &target_vector)
          : m_target_vector(target_vector)
        { m_target_vector.clear(); }

        void begin_particle(particle_number pn)
        { }

        void add_shape_value(grid_node_number gnn, double q_shapeval)
        { m_target_vector[gnn] += q_shapeval; }

        void end_particle(particle_number pn)
        { }
    };




    /** Reconstruction Target for the current density.
     */
    template<unsigned DimensionsVelocity>
    class j_target
    {
      private:
        dyn_vector &m_target_vector;
        const dyn_vector &m_velocities;
        double m_scale_factors[DimensionsVelocity];

      public:
        j_target(
            dyn_vector &target_vector, 
            const dyn_vector &velocities)
          : m_target_vector(target_vector), m_velocities(velocities)
        { 
          m_target_vector.clear();
          for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
            m_scale_factors[axis] = 0;
        }

        void begin_particle(particle_number pn)
        {
          for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
            m_scale_factors[axis] = m_velocities[pn*DimensionsVelocity+axis];
        }

        void add_shape_value(const grid_node_number gnn, const double q_shapeval)
        { 
          unsigned const base = gnn*DimensionsVelocity;
          for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
            m_target_vector[base+axis] += m_scale_factors[axis]*q_shapeval;
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

        void add_shape_value(grid_node_number gnn, double q_shapeval)
        { 
          m_target1.add_shape_value(gnn, q_shapeval);
          m_target2.add_shape_value(gnn, q_shapeval);
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




    // brick ------------------------------------------------------------------
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
      private:
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
          : m_number(number), m_start_index(start_index), 
          m_stepwidths(stepwidths), m_origin(origin), 
          m_origin_plus_half(origin+stepwidths/2),
          m_dimensions(dimensions)
        {
          // This ordering is what Visit expects by default and calls
          // "row-major" (I suppose Y-major).

          m_strides.resize(m_dimensions.size());
          unsigned i = 0;
          unsigned current_stride = 1;
          while (i < m_dimensions.size())
          {
            m_strides[i] = current_stride;
            current_stride *= m_dimensions[i];
            ++i;
          }
        }

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

        unsigned node_count() const
        {
          unsigned result = 1;
          BOOST_FOREACH(unsigned d, m_dimensions)
            result *= d;
          return result;
        }

        bounded_vector point(const bounded_int_vector &idx) const
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




        class iterator
        {
          private:
            const brick &m_brick;
            const bounded_int_box m_bounds;

            /* these fields must be updated in lock-step */
            bounded_int_vector m_state;
            bounded_vector m_point;
            grid_node_number m_index;

          public:
            iterator(brick const &brk, const bounded_int_box &bounds)
              : m_brick(brk), m_bounds(bounds), m_state(bounds.m_lower), 
              m_point(brk.point(m_state)),
              m_index(brk.index(m_state))
            { 
            }

            const bounded_int_vector &operator*() const
            { return m_state; }

          private:
            iterator &inc_bottom_half(unsigned i)
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
                  return *this;
                }

                m_state[i] = m_bounds.m_lower[i];
                ++i;
              }

              // we overflowed everything back to the origin, meaning we're done.
              // no need to update m_point and m_index.

              m_state = m_bounds.m_upper;
              return *this;
            }


          public:
            iterator &operator++()
            {
              const unsigned i = 0; // silently assuming non-zero size here

              ++m_state[i];

              if (m_state[i] >= m_bounds.m_upper[i])
              {
                m_state[i] = m_bounds.m_lower[i];

                // split off the non-default case bottom half of this routine in the
                // hope that at least the fast default case will get inlined.
                return inc_bottom_half(i); 
              }

              m_point[i] += m_brick.m_stepwidths[i];
              m_index += m_brick.m_strides[i];
              return *this;
            }

            bool at_end() const
            { return m_state[0] >= m_bounds.m_upper[0]; }

            grid_node_number index() const
            { return m_index; }

            const bounded_vector &point() const
            { return m_point; }
        };

        iterator get_iterator(bounded_int_box const &bounds) const
        {
          return iterator(*this, bounds);
        }

        iterator get_iterator(bounded_box const &bounds) const
        {
          return iterator(*this, index_range(bounds));
        }
    };




    struct element_on_grid
    {
      std::vector<grid_node_number> m_grid_nodes;

      /** The interpolant matrix maps the values (in-order) on the element
       * to the structured grid values at indices m_grid_nodes.
       *
       * The general assumption is that #(element nodes) < #(grid points in element),
       * but this may be violated without much harm.
       */
      dyn_fortran_matrix m_interpolation_matrix;
    };




    // main reconstructor type ------------------------------------------------
    template <class PICAlgorithm>
    class type : public reconstructor_base
    {
      public:
        // member data --------------------------------------------------------
        shape_function   m_shape_function;

        std::vector<brick> m_bricks;
        std::vector<element_on_grid> m_elements_on_grid;
        boost::numeric::ublas::vector<brick_number> m_particle_brick_numbers;

        /** Each brick may have a number of "extra" points to resolve
         * situations where the structured nodes situated on an element
         * are not enough to separate element modes. m_extra_points
         * describes this set of points. Points are subsequent in
         * m_extra_points, each point taking up a number of scalars
         * equal to the mesh dimension.
         *
         * m_extra_point_brick_starts[brick_number] contains the point
         * index in m_extra_points at which the extra points for
         * brick_number are stored, up to (but not including)
         * m_extra_point_brick_starts[brick_number+1].
         *
         * extra point numbering starts at m_first_extra_point.
         */
        grid_node_number m_first_extra_point;
        dyn_vector m_extra_points;
        std::vector<unsigned> m_extra_point_brick_starts;

        // internals ------------------------------------------------------------
        unsigned grid_node_count() const
        {
          const unsigned extra_pts = m_extra_points.size() 
            / CONST_PIC_THIS->m_mesh_data.m_dimensions;
          if (m_bricks.size() == 0)
            return extra_pts+0;
          else
            return extra_pts+m_bricks.back().start_index() + m_bricks.back().node_count();
        }




        template <class Target>
        bool reconstruct_particle_on_one_brick(Target tgt, 
            const brick &brk, 
            const bounded_vector &center,
            const bounded_box &particle_box,
            double charge,
            bool *complete = 0) const
        {
          const bounded_box intersect_box = 
            brk.bounding_box().intersect(particle_box);

          if (complete)
            *complete = intersect_box == particle_box;

          if (intersect_box.is_empty())
            return false;

          const bounded_int_box particle_brick_index_box = 
            brk.index_range(intersect_box);

          brick::iterator it(brk, particle_brick_index_box);

          while (!it.at_end())
          {
            tgt.add_shape_value(it.index(), 
                charge*m_shape_function(center-it.point()));
            ++it;
          }

          // now treat extra points for this brick
          const unsigned extra_stop 
            = m_extra_point_brick_starts[brk.number()+1];
          const grid_node_number fep = m_first_extra_point;
          const unsigned mdim = CONST_PIC_THIS->m_mesh_data.m_dimensions;

          for (unsigned extra_i = m_extra_point_brick_starts[brk.number()];
              extra_i < extra_stop; ++extra_i)
            tgt.add_shape_value(fep+extra_i, 
                charge*m_shape_function(
                  center-subrange(
                    m_extra_points, 
                    extra_i*mdim,
                    (extra_i+1)*mdim)));
          
          return true;
        }




        template <class Target>
        void reconstruct_single_particle_with_cache(Target tgt,
            particle_number pn,
            bounded_vector const &center,
            bounded_box const &particle_box)
        {
          const double charge = CONST_PIC_THIS->m_charges[pn];
          brick_number &bn_cache(m_particle_brick_numbers[pn]);
          const brick &last_brick = m_bricks[bn_cache];

          bool is_complete;
          bool does_intersect_cached = 
            reconstruct_particle_on_one_brick(
                tgt, last_brick, center, particle_box, charge,
                &is_complete);

          if (!is_complete)
          {
            BOOST_FOREACH(brick const &brk, m_bricks)
            {
              // don't re-target the cached brick
              if (brk.number() == last_brick.number())
                continue;

              if (reconstruct_particle_on_one_brick(
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
          BOOST_FOREACH(brick const &brk, m_bricks)
            reconstruct_particle_on_one_brick(
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
        void reconstruct_densities_on_grid_target(Target tgt)
        {
          const unsigned dim_x = CONST_PIC_THIS->get_dimensions_pos();
          const unsigned dim_m = CONST_PIC_THIS->m_mesh_data.m_dimensions;

          using boost::numeric::ublas::scalar_vector;
          const scalar_vector<double> shape_extent(
              dim_m, m_shape_function.radius());

          for (particle_number pn = 0; pn < CONST_PIC_THIS->m_particle_count; ++pn)
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




        template <class FromVec, class ToVec>
        void remap_grid_to_mesh(const FromVec &from, ToVec &to, 
            const unsigned offset=0, const unsigned increment=1) const
        {
          BOOST_FOREACH(const mesh_data::element_info &el, 
              CONST_PIC_THIS->m_mesh_data.m_element_info)
          {
            const element_on_grid &eog(m_elements_on_grid[el.m_id]);

            dyn_vector grid_values(eog.m_grid_nodes.size());

            {
              dyn_vector::iterator gv_it = grid_values.begin();
              BOOST_FOREACH(grid_node_number gnn, eog.m_grid_nodes)
                *gv_it++ = from[offset + gnn*increment];
            }

            /*
            noalias(subrange(to, el.m_start, el.m_end)) = prod(
                eog.m_interpolation_matrix, grid_values);
                */
            {
              const dyn_fortran_matrix &matrix = eog.m_interpolation_matrix;
              using namespace boost::numeric::bindings;
              using blas::detail::gemv;
              gemv(
                  'N',
                  eog.m_interpolation_matrix.size1(),
                  eog.m_interpolation_matrix.size2(),
                  /*alpha*/ 1,
                  traits::matrix_storage(matrix),
                  traits::leading_dimension(matrix),

                  traits::vector_storage(grid_values), /*incx*/ 1,

                  /*beta*/ 0,
                  traits::vector_storage(to) + el.m_start*increment + offset, 
                  /*incy*/ increment);
            }
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




        py_vector get_grid_rho()
        {
          dyn_vector grid_rho(grid_node_count());

          rho_target rho_tgt(grid_rho);
          reconstruct_densities_on_grid_target(rho_tgt);
          return grid_rho;
        }




        py_vector get_grid_j(py_vector const &velocities)
        {
          dyn_vector grid_j(
              CONST_PIC_THIS->get_dimensions_velocity()
              * grid_node_count());

          dyn_vector velocities_copy(velocities);
          j_target<PICAlgorithm::dimensions_velocity> 
            j_tgt(grid_j, velocities_copy);
          reconstruct_densities_on_grid_target(j_tgt);
          return grid_j;
        }




        // public interface -----------------------------------------------------
        static const std::string get_name()
        { return "Grid"; }

        void perform_reconstructor_upkeep()
        { }

        void reconstruct_densities(
            py_vector rho, py_vector j, const py_vector &velocities)
        {
          if (rho.size() != PIC_THIS->m_mesh_data.node_count())
            throw std::runtime_error("rho field does not have the correct size");
          if (j.size() != PIC_THIS->m_mesh_data.node_count() *
              PIC_THIS->get_dimensions_velocity())
            throw std::runtime_error("j field does not have the correct size");

          const unsigned gnc = grid_node_count();
          const unsigned vdim = CONST_PIC_THIS->get_dimensions_velocity();

          dyn_vector v(velocities);
          dyn_vector grid_rho(gnc);
          dyn_vector grid_j(vdim * gnc);

          rho_target rho_tgt(grid_rho);
          typedef j_target<PICAlgorithm::dimensions_velocity> j_tgt_t;
          j_tgt_t j_tgt(grid_j, v);
          chained_target<rho_target, j_tgt_t>
              tgt(rho_tgt, j_tgt);

          reconstruct_densities_on_grid_target(tgt);

          remap_grid_to_mesh(grid_rho, rho);
          for (unsigned i = 0; i < vdim; ++i)
            remap_grid_to_mesh(grid_j, j, i, vdim);
        }




        void reconstruct_j(py_vector j, const py_vector &velocities)
        {
          if (j.size() != PIC_THIS->m_mesh_data.node_count() *
              PIC_THIS->get_dimensions_velocity())
            throw std::runtime_error("j field does not have the correct size");

          const unsigned gnc = grid_node_count();
          const unsigned vdim = CONST_PIC_THIS->get_dimensions_velocity();

          dyn_vector v(velocities);
          dyn_vector grid_j(vdim * gnc);

          typedef j_target<PICAlgorithm::dimensions_velocity> j_tgt_t;
          j_tgt_t j_tgt(grid_j, v);

          reconstruct_densities_on_grid_target(j_tgt);

          for (unsigned i = 0; i < vdim; ++i)
            remap_grid_to_mesh(grid_j, j, i, vdim);
        }




        void reconstruct_rho(py_vector rho)
        {
          if (rho.size() != PIC_THIS->m_mesh_data.node_count())
            throw std::runtime_error("rho field does not have the correct size");

          dyn_vector grid_rho(grid_node_count());

          rho_target rho_tgt(grid_rho);
          reconstruct_densities_on_grid_target(rho_tgt);

          remap_grid_to_mesh(grid_rho, rho);
        }
    };
  };
}




#endif
