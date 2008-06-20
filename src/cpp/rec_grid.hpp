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
#include <boost/random.hpp>
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
#include "grid.hpp"




namespace pyrticle
{
  struct grid_reconstructor 
  {
    // specialized targets ----------------------------------------------------
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

        void add_shape_value(grid_node_number gnn, double q_shapeval)
        { m_target[gnn] += q_shapeval; }

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

        void add_shape_value(const grid_node_number gnn, const double q_shapeval)
        { 
          unsigned const base = gnn*DimensionsVelocity;
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




    typedef unsigned brick_number;




    struct rec_brick : public brick
    {
      private:
        static const unsigned point_origin_count = 3;
        brick_number m_number;
        typedef brick super;
        std::vector<bounded_vector> m_point_origins;
        std::vector<bounded_vector> m_axis0_offsets;

      public:
        rec_brick(
            brick_number number,
            grid_node_number start_index,
            bounded_vector stepwidths,
            bounded_vector origin,
            bounded_int_vector dimensions)
          : super(start_index, stepwidths, origin, dimensions),
          m_number(number)
        { 
          boost::variate_generator<
            boost::mt19937, 
            boost::uniform_real<> > offset_rng(
                boost::mt19937(), boost::uniform_real<>(-0.1,0.1));

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

        brick_number number() const
        { return m_number; }

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

        class iterator : public brick_iterator<rec_brick>
        {
          private:
            typedef brick_iterator<rec_brick> super;

          protected:
            unsigned m_origin_index;

          public:
            iterator(rec_brick const &brk, 
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




    struct element_on_grid
    {
      mesh_data::element_number m_element_number;
      std::vector<grid_node_number> m_grid_nodes;
      py_vector m_weight_factors;

      /** The interpolant matrix maps the values (in-order) on the element
       * to the structured grid values at indices m_grid_nodes.
       *
       * The general assumption is that #(element nodes) < #(grid points in element).
       */
      dyn_fortran_matrix m_interpolation_matrix;
      py_fortran_matrix m_inverse_interpolation_matrix;
    };




    // main reconstructor type ------------------------------------------------
    template <class PICAlgorithm>
    class type : public reconstructor_base
    {
      public:
        // member data --------------------------------------------------------
        shape_function   m_shape_function;

        std::vector<rec_brick> m_bricks;
        std::vector<element_on_grid> m_elements_on_grid;
        mutable unsigned m_max_el_grid_values;
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

        /** Average groups for continuity enforcement. 
         *
         * For each averaging group consisting of node indices,
         * the corresponding nodes of the result of elementwise remapping
         * are averaged and written back.
         *
         * It is implied that the first average group starts at 0.
         */
        std::vector<mesh_data::node_number> m_average_groups;
        std::vector<unsigned> m_average_group_starts;





        // construction -------------------------------------------------------
        type()
          : m_max_el_grid_values(0)
        { }




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




        py_vector find_points_in_element(element_on_grid &eog, double scaled_tolerance) const
        {
          const unsigned mdims = CONST_PIC_THIS->m_mesh_data.m_dimensions;
          const mesh_data::element_info &el = 
            CONST_PIC_THIS->m_mesh_data.m_element_info[
            eog.m_element_number];

          bounded_box el_bbox = CONST_PIC_THIS->m_mesh_data.element_bounding_box(
              eog.m_element_number);

          using boost::numeric::ublas::scalar_vector;
          const scalar_vector<double> tolerance_vec(mdims, scaled_tolerance);
          el_bbox.m_lower -= tolerance_vec;
          el_bbox.m_upper += tolerance_vec;

          unsigned gnc = grid_node_count();

          std::vector<bounded_vector> points;
          std::vector<double> weights;

          // For each element, find all structured points inside the element.
          BOOST_FOREACH(rec_brick const &brk, m_bricks)
          {
            const double dV = brk.cell_volume();

            bounded_box brk_and_el = brk.bounding_box().intersect(el_bbox);
            if (brk_and_el.is_empty())
                continue;

            const bounded_int_box el_brick_index_box = 
              brk.index_range(brk_and_el);

            rec_brick::iterator it(brk, el_brick_index_box);

            while (!it.at_end())
            {
              bool in_el = true;

              bounded_vector point = it.point();

              BOOST_FOREACH(const mesh_data::face_info &f, el.m_faces)
                if (inner_prod(f.m_normal, point) 
                    - f.m_face_plane_eqn_rhs > scaled_tolerance)
                {
                  in_el = false;
                  break;
                }

              if (in_el)
              {
                points.push_back(point);

                grid_node_number gni = it.index();
                if (gni >= gnc)
                  throw std::runtime_error("rec_grid: structured point index out of bounds");
                eog.m_grid_nodes.push_back(gni);
                weights.push_back(sqrt(dV));
              }

              ++it;
            }
          }

          eog.m_weight_factors.resize(weights.size());
          std::copy(weights.begin(), weights.end(), eog.m_weight_factors.begin());

          npy_intp dims[] = { points.size(), mdims };
          py_vector points_copy(2, dims);
          for (unsigned i = 0; i < points.size(); ++i)
            subrange(points_copy, mdims*i, mdims*(i+1)) = points[i];
          return points_copy;
        }




        template <class Target>
        bool reconstruct_particle_on_one_brick(Target tgt, 
            const rec_brick &brk, 
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

          rec_brick::iterator it(brk, particle_brick_index_box);

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
          const rec_brick &last_brick = m_bricks[bn_cache];

          bool is_complete;
          bool does_intersect_cached = 
            reconstruct_particle_on_one_brick(
                tgt, last_brick, center, particle_box, charge,
                &is_complete);

          if (!is_complete)
          {
            BOOST_FOREACH(rec_brick const &brk, m_bricks)
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
          BOOST_FOREACH(rec_brick const &brk, m_bricks)
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




        void remap_grid_to_mesh(const py_vector from, py_vector to, 
            const unsigned offset=0, const unsigned increment=1) const
        {
          const py_vector::const_iterator from_it = from.begin();

          if (m_max_el_grid_values == 0)
          {
            BOOST_FOREACH(const element_on_grid &eog, m_elements_on_grid)
              m_max_el_grid_values = std::max<unsigned>(
                  m_max_el_grid_values,
                  eog.m_grid_nodes.size());
          }

          dyn_vector grid_values(m_max_el_grid_values);

          BOOST_FOREACH(const element_on_grid &eog, m_elements_on_grid)
          {
            const mesh_data::element_info &el = 
              CONST_PIC_THIS->m_mesh_data.m_element_info[
              eog.m_element_number];

            // pick values off the grid
            const py_vector::const_iterator weights = eog.m_weight_factors.begin();

            for (unsigned i = 0; i < eog.m_grid_nodes.size(); ++i)
              grid_values[i] = weights[i]
                * from_it[offset + eog.m_grid_nodes[i]*increment];

            // and apply the interpolation matrix
            {
              const dyn_fortran_matrix &matrix = eog.m_interpolation_matrix;
              using namespace boost::numeric::bindings;
              using blas::detail::gemv;
              gemv(
                  'N',
                  matrix.size1(),
                  matrix.size2(),
                  /*alpha*/ 1,
                  traits::matrix_storage(matrix),
                  traits::leading_dimension(matrix),

                  traits::vector_storage(grid_values), /*incx*/ 1,

                  /*beta*/ 1,
                  traits::vector_storage(to) + el.m_start*increment + offset, 
                  /*incy*/ increment);
            }
          }

          // cross-element continuity enforcement
          const py_vector::iterator to_it = to.begin();

          std::vector<unsigned>::const_iterator 
            ag_starts_first = m_average_group_starts.begin(),
            ag_starts_last = m_average_group_starts.end();

          unsigned ag_start = 0;
          unsigned ag_end = *ag_starts_first++;

          while (true)
          {
            double avg = 0;
            const std::pair<
              std::vector<mesh_data::node_number>::const_iterator,
              std::vector<mesh_data::node_number>::const_iterator
              > ag_range(
                m_average_groups.begin()+ag_start,
                m_average_groups.begin()+ag_end);
            BOOST_FOREACH(mesh_data::node_number nn, ag_range)
              avg += to_it[nn];
            avg /= (ag_end-ag_start);

            BOOST_FOREACH(mesh_data::node_number nn, ag_range)
              to_it[nn] = avg;

            if (ag_starts_first == ag_starts_last)
              break;

            ag_start = ag_end;
            ag_end = *ag_starts_first++;
          }
        }




        void remap_residual(const py_vector &from, py_vector to, 
            const unsigned offset=0, const unsigned increment=1) const
        {
          unsigned max_el_size = 0;
          BOOST_FOREACH(const mesh_data::element_info &el, 
              CONST_PIC_THIS->m_mesh_data.m_element_info)
            max_el_size = std::max(max_el_size, el.m_end-el.m_start);
          dyn_vector mesh_values(max_el_size);

          const py_vector::const_iterator from_it = from.begin();
          const py_vector::iterator to_it = to.begin();

          BOOST_FOREACH(const element_on_grid &eog, m_elements_on_grid)
          {
            if (!eog.m_inverse_interpolation_matrix.is_valid())
              break;

            const py_vector::const_iterator weights = eog.m_weight_factors.begin();
            dyn_vector grid_values(eog.m_grid_nodes.size());

            // pick values off the grid
            for (unsigned i = 0; i < eog.m_grid_nodes.size(); ++i)
              grid_values[i] = weights[i]
                * from_it[offset + eog.m_grid_nodes[i]*increment];

            // apply the interpolation matrix
            {
              const dyn_fortran_matrix &matrix = eog.m_interpolation_matrix;
              using namespace boost::numeric::bindings;
              using blas::detail::gemv;
              gemv(
                  'N',
                  matrix.size1(),
                  matrix.size2(),
                  /*alpha*/ 1,
                  traits::matrix_storage(matrix),
                  traits::leading_dimension(matrix),

                  traits::vector_storage(grid_values), /*incx*/ 1,

                  /*beta*/ 0,
                  traits::vector_storage(mesh_values), 
                  /*incy*/ 1);
            }

            // apply the inverse interpolation matrix
            {
              const dyn_fortran_matrix &matrix = eog.m_inverse_interpolation_matrix;
              using namespace boost::numeric::bindings;
              using blas::detail::gemv;
              gemv(
                  'N',
                  matrix.size1(),
                  matrix.size2(),
                  /*alpha*/ 1,
                  traits::matrix_storage(matrix),
                  traits::leading_dimension(matrix),

                  traits::vector_storage(mesh_values), 
                  /*incx*/ 1,

                  /*beta*/ 0,
                  traits::vector_storage(grid_values), /*incy*/ 1
                  );
            }

            // add squared residuals back onto the grid
            for (unsigned i = 0; i < eog.m_grid_nodes.size(); ++i)
              to_it[offset + eog.m_grid_nodes[i]*increment] +=
                square(grid_values[i] / weights[i]
                    - from_it[offset + eog.m_grid_nodes[i]*increment]);
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




        // gridded output -----------------------------------------------------
        boost::tuple<py_vector, py_vector> 
          reconstruct_grid_densities(const py_vector &velocities)
        {
          const unsigned gnc = grid_node_count();
          const unsigned vdim = CONST_PIC_THIS->get_dimensions_velocity();

          py_vector grid_rho(gnc);
          npy_intp dims[] = { gnc, vdim };
          py_vector grid_j(2, dims);

          rho_target<py_vector> rho_tgt(grid_rho);
          typedef j_target<PICAlgorithm::dimensions_velocity, py_vector, py_vector> 
            j_tgt_t;
          j_tgt_t j_tgt(grid_j, velocities);
          chained_target<rho_target<py_vector>, j_tgt_t>
              tgt(rho_tgt, j_tgt);

          reconstruct_densities_on_grid_target(tgt);

          return boost::make_tuple(grid_rho, grid_j);
        }




        py_vector reconstruct_grid_j(py_vector const &velocities)
        {
          const unsigned gnc = grid_node_count();
          const unsigned vdim = CONST_PIC_THIS->get_dimensions_velocity();

          npy_intp dims[] = { gnc, vdim };
          py_vector grid_j(2, dims);

          j_target<PICAlgorithm::dimensions_velocity, py_vector, py_vector> 
            j_tgt(grid_j, velocities);
          reconstruct_densities_on_grid_target(j_tgt);
          return grid_j;
        }




        py_vector reconstruct_grid_rho()
        {
          py_vector grid_rho(grid_node_count());

          rho_target<py_vector> rho_tgt(grid_rho);
          reconstruct_densities_on_grid_target(rho_tgt);
          return grid_rho;
        }




        // public interface -----------------------------------------------------
        static const std::string get_name()
        { return "Grid"; }
    };
  };
}




#endif
