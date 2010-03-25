// Pyrticle - Particle in Cell in Python
// Grid-based deposition
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
#include "dep_shape.hpp"
#include "grid.hpp"
#include "dep_grid_base.hpp"




namespace pyrticle
{
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




  template <class ParticleState, class ShapeFunction, class Brick>
  struct grid_depositor : 
    public grid_targets, 
    public grid_depositor_base<grid_depositor<ParticleState, ShapeFunction, Brick>,
      ParticleState, ShapeFunction, Brick, grid_depositor_base_state>
  {
    private:
      typedef grid_depositor_base<
        grid_depositor<ParticleState, ShapeFunction, Brick>, 
        ParticleState, ShapeFunction, Brick, grid_depositor_base_state>
          base_type;

    public:
      typedef typename base_type::particle_state particle_state;
      typedef typename base_type::brick_type brick_type;
      typedef typename base_type::depositor_state depositor_state;

      // member data --------------------------------------------------------
      std::vector<element_on_grid> m_elements_on_grid;
      mutable unsigned m_max_el_grid_values;

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
      std::vector<npy_uint32> m_extra_point_brick_starts;

      /** Average groups for continuity enforcement. 
       *
       * For each averaging group consisting of node indices,
       * the corresponding nodes of the result of elementwise remapping
       * are averaged and written back.
       *
       * It is implied that the first average group starts at 0.
       */
      std::vector<mesh_data::node_number> m_average_groups;
      std::vector<npy_uint32> m_average_group_starts;





      // construction -------------------------------------------------------
      grid_depositor(const mesh_data &md)
        : base_type(md), m_max_el_grid_values(0)
      { }




      // internals ------------------------------------------------------------
      unsigned grid_node_count_with_extra() const
      {
        const unsigned extra_pts = m_extra_points.size() 
          / this->m_mesh_data.m_dimensions;
        return extra_pts + this->grid_node_count();
      }




      py_vector find_points_in_element(element_on_grid &eog, double scaled_tolerance) const
      {
        const unsigned mdims = this->m_mesh_data.m_dimensions;
        const mesh_data::element_info &el = 
          this->m_mesh_data.m_element_info[
          eog.m_element_number];

        bounded_box el_bbox = this->m_mesh_data.element_bounding_box(
            eog.m_element_number);

        using boost::numeric::ublas::scalar_vector;
        const scalar_vector<double> tolerance_vec(mdims, scaled_tolerance);
        el_bbox.m_lower -= tolerance_vec;
        el_bbox.m_upper += tolerance_vec;

        unsigned gnc = grid_node_count_with_extra();

        std::vector<bounded_vector> points;
        std::vector<double> weights;

        // For each element, find all structured points inside the element.
        BOOST_FOREACH(Brick const &brk, this->m_bricks)
        {
          const double dV = brk.cell_volume();

          bounded_box brk_and_el = brk.bounding_box().intersect(el_bbox);
          if (brk_and_el.is_empty())
              continue;

          const bounded_int_box el_brick_index_box = 
            brk.index_range(brk_and_el);

          typename Brick::iterator it(brk, el_brick_index_box);

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




      // 1 particle + 1 brick -------------------------------------------------
      template <class Target>
      bool deposit_particle_on_one_brick(Target tgt, 
          const Brick &brk, 
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

        typename Brick::iterator it(brk, particle_brick_index_box);

        while (!it.at_end())
        {
          tgt.add_shape_value(it.index(), 
              charge*this->m_shape_function(center-it.point()));
          ++it;
        }

        // now treat extra points for this brick
        const unsigned extra_stop 
          = m_extra_point_brick_starts[brk.number()+1];
        const grid_node_number fep = m_first_extra_point;
        const unsigned mdim = this->m_mesh_data.m_dimensions;

        for (unsigned extra_i = m_extra_point_brick_starts[brk.number()];
            extra_i < extra_stop; ++extra_i)
          tgt.add_shape_value(fep+extra_i, 
              charge*this->m_shape_function(
                center-subrange(
                  m_extra_points, 
                  extra_i*mdim,
                  (extra_i+1)*mdim)));
        
        return true;
      }




      // other stuff ----------------------------------------------------------
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
            this->m_mesh_data.m_element_info[eog.m_element_number];

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
        const py_vector::iterator to_it = to.begin() + offset;

        std::vector<npy_uint32>::const_iterator 
          ag_starts_first = m_average_group_starts.begin(),
          ag_starts_last = m_average_group_starts.end();

        npy_uint32 ag_start = 0;
        npy_uint32 ag_end = *ag_starts_first++;

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
            avg += to_it[nn*increment];
          avg /= (ag_end-ag_start);

          BOOST_FOREACH(mesh_data::node_number nn, ag_range)
            to_it[nn*increment] = avg;

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
            this->m_mesh_data.m_element_info)
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




      // gridded output -----------------------------------------------------
      boost::tuple<py_vector, py_vector> 
        deposit_grid_densities(
            depositor_state &ds,
            const particle_state &ps,
            const py_vector &velocities,
            boost::python::slice const &pslice)
      {
        const unsigned gnc = grid_node_count_with_extra();
        const unsigned vdim = ps.vdim();

        py_vector grid_rho(gnc);
        npy_intp dims[] = { gnc, vdim };
        py_vector grid_j(2, dims);

        rho_target<py_vector> rho_tgt(grid_rho);
        typedef j_target<particle_state::m_vdim, py_vector, py_vector> 
          j_tgt_t;
        j_tgt_t j_tgt(grid_j, velocities);
        chained_target<rho_target<py_vector>, j_tgt_t>
            tgt(rho_tgt, j_tgt);

        deposit_densities_on_grid_target(ds, ps, tgt, pslice);

        return boost::make_tuple(grid_rho, grid_j);
      }




      py_vector deposit_grid_j(
          depositor_state &ds,
          const particle_state &ps,
          py_vector const &velocities,
          boost::python::slice const &pslice)
      {
        const unsigned gnc = grid_node_count_with_extra();
        const unsigned vdim = ps.vdim();

        npy_intp dims[] = { gnc, vdim };
        py_vector grid_j(2, dims);

        j_target<particle_state::m_vdim, py_vector, py_vector> 
          j_tgt(grid_j, velocities);
        deposit_densities_on_grid_target(ds, ps, j_tgt, pslice);
        return grid_j;
      }




      py_vector deposit_grid_rho(
          depositor_state &ds,
          const particle_state &ps,
          boost::python::slice const &pslice)
      {
        py_vector grid_rho(grid_node_count_with_extra());

        rho_target<py_vector> rho_tgt(grid_rho);
        deposit_densities_on_grid_target(ds, ps, rho_tgt, pslice);
        return grid_rho;
      }
  };
}




#endif
