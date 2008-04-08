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
#include "rec_shape.hpp"
#include <pyublas/unary_op.hpp>
#include <boost/python.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/bindings/lapack/gesvd.hpp>




namespace pyrticle
{
  typedef unsigned grid_node_number;




  namespace grid_targets 
  {
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

        void add_shape_value(grid_node_number gnn, double q_shapeval)
        { 
          for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
            m_target_vector[gnn*DimensionsVelocity+axis] += m_scale_factors[axis]*q_shapeval;
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
    inline
    chained_target<T1, T2> 
    make_chained_target(T1 &target1, T2 &target2)
    {
      return chained_target<T1, T2>(target1, target2);
    }
  }




  struct grid_reconstructor 
  {
    template <class PICAlgorithm>
    class type : public reconstructor_base
    {
      public:
        class brick_iterator
        {
          private:
            const bounded_int_box &m_bounds;
            bounded_int_vector m_state;

          public:
            brick_iterator(const bounded_int_box &bounds)
              : m_bounds(bounds), m_state(bounds.first)
            { }

            const bounded_int_vector &operator*() const
            { return m_state; }

            brick_iterator operator++()
            {
              int i = m_state.size();

              while (i > 0)
              {
                --i;
                ++m_state[i];
                if (m_state[i] < m_bounds.second[i])
                  return *this;

                m_state[i] = m_bounds.first[i];
              }

              return *this;
            }

            bool at_end() const
            { return m_state[0] >= m_bounds.second[0]; }
        };




        struct brick
        {
          public:
            grid_node_number m_start_index;
            bounded_vector m_stepwidths;
            bounded_vector m_origin;

            /** This is the number of points in each dimension. If it is 2, then the
             * indices 0 and 1 exist in that dimension.
             */
            bounded_int_vector m_dimensions;

            unsigned node_count() const
            {
              unsigned result = 1;
              BOOST_FOREACH(unsigned d, m_dimensions)
                result *= d;
              return result;
            }

            bounded_vector point(const bounded_int_vector &idx) const
            { return m_origin + element_prod(idx, m_stepwidths); }

            grid_node_number index(const bounded_int_vector &idx) const
            { return m_start_index + inner_prod(idx, m_dimensions); }

            bounded_box bounding_box() const
            {
              return std::make_pair(
                  m_origin,
                  m_origin + element_prod(m_dimensions, m_stepwidths) - m_stepwidths
                  );
            }

          private:
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

          public:
            bounded_int_box index_range(const bounded_box &bbox) const
            {
              return bounded_int_box(
                  pyublas::unary_op<int_floor>::apply(
                    element_div(bbox.first-m_origin, m_stepwidths)),
                  pyublas::unary_op<int_ceil>::apply(
                    element_div(bbox.second-m_origin, m_stepwidths))
                  );
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
          dyn_matrix m_interpolant_pseudo_inv;
        };




        typedef unsigned brick_number;
        typedef unsigned brick_node_number;




        // member data --------------------------------------------------------
        shape_function   m_shape_function;

        std::vector<brick> m_bricks;
        std::vector<element_on_grid> m_elements_on_grid;
        boost::numeric::ublas::vector<brick_number> m_particle_brick_numbers;

        // setup interface ------------------------------------------------------
        void add_brick(
            py_vector stepwidths,
            py_vector origin,
            pyublas::numpy_vector<unsigned> dims)
        {
          brick new_brick;
          new_brick.m_start_index = grid_node_count();
          new_brick.m_stepwidths = stepwidths;
          new_brick.m_origin = origin;
          new_brick.m_dimensions = dims;

          m_bricks.push_back(new_brick);
        }




        void commit_bricks(py_matrix inv_vander_t, boost::python::object basis)
        {
          BOOST_FOREACH(const mesh_data::element_info &el, 
              CONST_PIC_THIS->m_mesh_data.m_element_info)
          {
            element_on_grid eog;
            bounded_box el_bbox = CONST_PIC_THIS->m_mesh_data.element_bounding_box(el.m_id);

            // for each element, and each brick, figure it out the points in the brick
            // that fall inside the element. Add them to point_coordinates and 
            // eog.m_grid_nodes.
            std::vector<bounded_vector> point_coordinates;

            BOOST_FOREACH(const brick &brk, m_bricks)
            {
              bounded_box brick_bbox = brk.bounding_box();
              bool does_intersect;
              bounded_int_box el_brick_index_box = 
                brk.index_range(
                    intersect(brick_bbox, el_bbox, &does_intersect));

              if (!does_intersect)
                continue;

              brick_iterator it(el_brick_index_box);

              while (!it.at_end())
              {
                bounded_vector point = brk.point(*it);
                if (CONST_PIC_THIS->m_mesh_data.is_in_element(el.m_id, point))
                {
                  point_coordinates.push_back(point);
                  eog.m_grid_nodes.push_back(brk.index(*it));
                }

                ++it;
              }
            }

            // build the interpolant matrix that maps a element nodal values to
            // values on the brick grid.
            unsigned 
              rows = point_coordinates.size(),
                   cols = el.m_end-el.m_start;
            dyn_fortran_matrix interpolant_mat(rows, cols);
            interpolant_mat.clear();

            for (unsigned i = 0; i < point_coordinates.size(); ++i)
            {
              py_vector pt(point_coordinates[i]);

              using boost::python::object;

              dyn_vector basis_values_at_pt(el.m_end-el.m_start);
              BOOST_FOREACH(object basis_func,
                  std::make_pair(
                    boost::python::stl_input_iterator<object>(basis),
                    boost::python::stl_input_iterator<object>()
                    ))
                basis_values_at_pt[i] = boost::python::extract<double>(
                    pt.to_python().get());

              row(interpolant_mat, i) = prod(
                  inv_vander_t, basis_values_at_pt);
            }

            // now find the pseudoinverse of the interpolant matrix
            using boost::numeric::bindings::lapack::gesvd;

            dyn_fortran_matrix u(rows, rows), vt(cols, cols);
            dyn_vector s(std::min(rows, cols));
            u.clear();
            vt.clear();
            s.clear();

            gesvd(interpolant_mat, s, u, vt);

            // from http://en.wikipedia.org/wiki/Moore-Penrose_pseudoinverse
            // Matlab apparently uses this threshold
            const double threshold = 
              std::numeric_limits<double>::epsilon() *
              std::max(rows, cols) *
              norm_inf(s);

            boost::numeric::ublas::diagonal_matrix<double> inv_s(s.size());
            for (unsigned i = 0; i < s.size(); ++i)
              if (fabs(s[i]) > threshold)
                inv_s(i, i) = 1/s[i];
              else
                inv_s(i, i) = 0;

            noalias(eog.m_interpolant_pseudo_inv) = prod(
                trans(vt), dyn_matrix(prod(inv_s, trans(u))));

            m_elements_on_grid.push_back(eog);
          }
        }




        unsigned grid_node_count() const
        {
          if (m_bricks.size() == 0)
            return 0;
          else
            return m_bricks.back().m_start_index + m_bricks.back().node_count();
        }




        template <class Target>
        void reconstruct_particle_on_one_brick(Target tgt, 
            const brick &brk, 
            const bounded_vector &center,
            const bounded_box &intersect_box)
        {
          bounded_int_box particle_brick_index_box = brk.index_range(intersect_box);

          brick_iterator it(particle_brick_index_box);

          while (!it.at_end())
          {
            tgt.add_shape_value(brk.index(*it), 
                m_shape_function(center- brk.point(*it)));
            ++it;
          }
        }




        template <class Target>
        void reconstruct_densities_on_grid_target(Target tgt)
        {
          const unsigned dim = CONST_PIC_THIS->get_dimensions_pos();

          for (particle_number pn = 0; pn < CONST_PIC_THIS->m_particle_count; ++pn)
          {
            tgt.begin_particle(pn);

            const bounded_vector center = subrange(
                CONST_PIC_THIS->m_positions, pn*dim, (pn+1)*dim);

            using boost::numeric::ublas::scalar_vector;

            bounded_box particle_box(
                center - scalar_vector<double>(dim, m_shape_function.radius()),
                center + scalar_vector<double>(dim, m_shape_function.radius()));

            const brick_number last_brick_number = m_particle_brick_numbers[pn];
            const brick &last_brick = m_bricks[last_brick_number];

            bool does_intersect;
            bounded_box brick_and_particle = intersect(
                last_brick.bounding_box(), particle_box, &does_intersect);

            bool consider_other_bricks = !does_intersect;
            bool reset_brick_cache = !does_intersect;
            if (does_intersect)
            {
              reconstruct_particle_on_one_brick(
                  tgt, last_brick, center, brick_and_particle);

              if (brick_and_particle != particle_box)
                consider_other_bricks = true;
            }

            if (consider_other_bricks)
            {
              brick_number bn = 0;
              BOOST_FOREACH(brick const &brk, m_bricks)
              {
                // don't re-target the cached brick
                if (bn == last_brick_number)
                  continue;

                bool does_intersect;
                bounded_box brick_and_particle = intersect(
                    brk.bounding_box(), particle_box, &does_intersect);

                if (does_intersect)
                {
                  if (reset_brick_cache)
                    m_particle_brick_numbers[pn] = bn;

                  reconstruct_particle_on_one_brick(
                      tgt, brk, center, brick_and_particle);
                }

                ++bn;
              }
            }

            tgt.end_particle(pn);
          }
        }




        template <class FromVec, class ToVec>
        void remap_grid_to_mesh(const FromVec &from, ToVec &to) const
        {
          BOOST_FOREACH(const mesh_data::element_info &el, 
              CONST_PIC_THIS->m_mesh_data.m_element_info)
          {
            const element_on_grid &eog(m_elements_on_grid[el.m_id]);

            dyn_vector grid_values(eog.m_grid_nodes.size());

            {
              dyn_vector::iterator gv_it = grid_values.begin();
              BOOST_FOREACH(grid_node_number gnn, eog.m_grid_nodes)
                *gv_it++ = from[gnn];
            }

            subrange(to, el.m_start, el.m_end) = prod(
                eog.m_interpolant_pseudo_inv, grid_values);
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

          using namespace grid_targets;

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
          {
            boost::numeric::ublas::vector_slice<py_vector> j_mesh_slice(
                j, boost::numeric::ublas::slice(
                  0, vdim, PIC_THIS->m_mesh_data.node_count()));

            remap_grid_to_mesh(
                subslice(grid_j, 0, vdim, gnc),
                j_mesh_slice);
          }
        }




        void reconstruct_j(py_vector j, const py_vector &velocities)
        {
          if (j.size() != PIC_THIS->m_mesh_data.node_count() *
              PIC_THIS->get_dimensions_velocity())
            throw std::runtime_error("j field does not have the correct size");

          using namespace grid_targets;

          const unsigned gnc = grid_node_count();
          const unsigned vdim = CONST_PIC_THIS->get_dimensions_velocity();

          dyn_vector v(velocities);
          dyn_vector grid_j(vdim * gnc);

          typedef j_target<PICAlgorithm::dimensions_velocity> j_tgt_t;
          j_tgt_t j_tgt(grid_j, v);

          reconstruct_densities_on_grid_target(j_tgt);

          for (unsigned i = 0; i < vdim; ++i)
          {
            boost::numeric::ublas::vector_slice<py_vector> j_mesh_slice(
                j, boost::numeric::ublas::slice(
                  0, vdim, PIC_THIS->m_mesh_data.node_count()));

            remap_grid_to_mesh(
                subslice(grid_j, 0, vdim, gnc),
                j_mesh_slice);
          }
        }




        void reconstruct_rho(py_vector rho)
        {
          if (rho.size() != PIC_THIS->m_mesh_data.node_count())
            throw std::runtime_error("rho field does not have the correct size");

          using namespace grid_targets;

          const unsigned gnc = grid_node_count();

          dyn_vector grid_rho(gnc);

          rho_target rho_tgt(grid_rho);
          reconstruct_densities_on_grid_target(rho_tgt);

          remap_grid_to_mesh(grid_rho, rho);
        }
    };
  };
}




#endif
