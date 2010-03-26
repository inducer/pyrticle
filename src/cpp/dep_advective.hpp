// Pyrticle - Particle in Cell in Python
// Deposition based on advected shapes
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





#ifndef _AFAYYTAA_PYRTICLE_REC_ADVECTIVE_HPP_INCLUDED
#define _AFAYYTAA_PYRTICLE_REC_ADVECTIVE_HPP_INCLUDED




#include <vector>
#include <numeric>
#include <algorithm>
#include <boost/array.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/blas/blas3.hpp>
#include <boost/typeof/std/utility.hpp>
#include <boost/unordered_map.hpp>
#include <pyublas/elementwise_op.hpp>
#include <hedge/face_operators.hpp>
#include "tools.hpp"
#include "meshdata.hpp"
#include "dep_target.hpp"
#include "dep_shape.hpp"
#include "element_finder.hpp"




namespace pyrticle
{
  template <class ParticleState, class ShapeFunction>
  struct advective_depositor
  {
    public:
      // member types -------------------------------------------------------
      typedef ShapeFunction shape_function;
      typedef ParticleState particle_state;

      static const unsigned max_faces = 4;

      struct active_element
      {
        const mesh_data::element_info *m_element_info;
        boost::array<mesh_data::element_number,
          advective_depositor::max_faces> m_connections;
        unsigned m_start_index;
        unsigned m_min_life;

        active_element()
          : m_element_info(0)
        {
          for (unsigned i = 0; i < advective_depositor::max_faces; i++)
            m_connections[i] = mesh_data::INVALID_ELEMENT;
        }
      };

      struct advected_particle
      {
        shape_function                m_shape_function;
        std::vector<active_element>   m_elements;

        active_element *find_element(mesh_data::element_number en)
        {
          if (en == mesh_data::INVALID_ELEMENT)
            return 0;
          BOOST_FOREACH(active_element &el, m_elements)
          {
            if (el.m_element_info->m_id == en)
              return &el;
          }
          return 0;
        }

        const active_element *find_element(mesh_data::element_number en) const
        {
          if (en == mesh_data::INVALID_ELEMENT)
            return 0;
          BOOST_FOREACH(const active_element &el, m_elements)
          {
            if (el.m_element_info->m_id == en)
              return &el;
          }
          return 0;
        }

      };




      // particle state for advective -----------------------------------------
      struct depositor_state
      {
        unsigned                        m_active_elements;
        std::vector<unsigned>           m_freelist;

        std::vector<advected_particle>  m_advected_particles;
        dyn_vector                      m_rho;

        boost::shared_ptr<number_shift_listener> m_rho_dof_shift_listener;

        event_counter m_element_activation_counter, m_element_kill_counter;




        depositor_state()
          : m_active_elements(0)
        { }

        template <class NewRhoExpr>
        depositor_state(
            const depositor_state &src,
            const NewRhoExpr &new_rho,
            boost::shared_ptr<number_shift_listener> rho_dof_sl)
          : m_active_elements(src.m_active_elements),
          m_freelist(src.m_freelist),
          m_advected_particles(src.m_advected_particles),
          m_rho(new_rho),
          m_rho_dof_shift_listener(rho_dof_sl),
          m_element_activation_counter(src.m_element_activation_counter),
          m_element_kill_counter(src.m_element_kill_counter)
        { }

        virtual ~depositor_state()
        { }

        unsigned count_advective_particles() const
        {
          return m_advected_particles.size();
        }

        void resize_rho(unsigned new_size)
        {
          unsigned old_size = m_rho.size();
          unsigned copy_size = std::min(new_size, old_size);

          dyn_vector new_rho(new_size);
          subrange(new_rho, 0, copy_size) = subrange(m_rho, 0, copy_size);
          new_rho.swap(m_rho);
        }


        void clear()
        {
          m_advected_particles.clear();
          m_freelist.clear();
          m_active_elements = 0;
        }
      };




      // member data --------------------------------------------------------
      const mesh_data                 &m_mesh_data;

      unsigned                        m_faces_per_element;
      unsigned                        m_dofs_per_element;

      py_matrix                       m_mass_matrix;
      dyn_vector                      m_integral_weights;
      dyn_fortran_matrix              m_inverse_mass_matrix;
      py_matrix                       m_face_mass_matrix;
      dyn_vector                      m_face_integral_weights;
      dyn_fortran_matrix              m_filter_matrix;

      std::vector<dyn_fortran_matrix> m_local_diff_matrices;

      typedef hedge::face_pair<hedge::straight_face> face_pair_type;
      typedef hedge::face_group<face_pair_type> face_group_type;

      boost::shared_ptr<face_group_type > m_int_face_group;
      boost::shared_ptr<face_group_type > m_bdry_face_group;

      struct face_pair_locator
      {
        face_group_type const *m_face_group;
        face_pair_type const *m_face_pair;

        face_pair_locator()
          : m_face_group(0), m_face_pair(0)
        { }

        face_pair_locator(
            face_group_type     const &face_group,
            face_pair_type      const &face_pair
            )
          : m_face_group(&face_group), m_face_pair(&face_pair)
        { }
      };

      boost::unordered_map<mesh_data::el_face, face_pair_locator> m_el_face_to_face_pair_locator;

      double m_activation_threshold;
      double m_kill_threshold;
      double m_upwind_alpha;




      // initialization -----------------------------------------------------
      advective_depositor(
          const mesh_data &md,
          unsigned faces_per_element,
          unsigned dofs_per_element,
          const py_matrix &mass_matrix,
          const py_matrix &inverse_mass_matrix,
          const py_matrix &filter_matrix,
          const py_matrix &face_mass_matrix,
          boost::shared_ptr<face_group_type> int_face_group,
          boost::shared_ptr<face_group_type> bdry_face_group,
          double activation_threshold,
          double kill_threshold,
          double upwind_alpha
          )
        : m_mesh_data(md), m_faces_per_element(0), m_dofs_per_element(0),
        m_activation_threshold(0), m_kill_threshold(0),
        m_upwind_alpha(1)
      {
        m_faces_per_element = faces_per_element;
        m_dofs_per_element = dofs_per_element;

        m_mass_matrix = mass_matrix;
        m_integral_weights = prod(m_mass_matrix,
            boost::numeric::ublas::scalar_vector<double>
            (m_mass_matrix.size1(), 1));
        m_inverse_mass_matrix = inverse_mass_matrix;

        m_filter_matrix = filter_matrix;

        m_face_mass_matrix = face_mass_matrix;
        m_face_integral_weights = prod(m_face_mass_matrix,
            boost::numeric::ublas::scalar_vector<double>
            (m_face_mass_matrix.size1(), 1));

        m_int_face_group = int_face_group;
        m_bdry_face_group = bdry_face_group;

        // build m_el_face_to_face_pair_locator
        BOOST_FOREACH(const hedge::face_pair<hedge::straight_face> &fp,
            int_face_group->face_pairs)
        {
          const hedge::straight_face &f = fp.int_side;
          m_el_face_to_face_pair_locator
            [std::make_pair(f.element_id, f.face_id)] =
            face_pair_locator(*int_face_group, fp);

          const hedge::straight_face &ext_f = fp.ext_side;
          m_el_face_to_face_pair_locator
            [std::make_pair(ext_f.element_id, ext_f.face_id)] =
            face_pair_locator(*int_face_group, fp);
        }

        BOOST_FOREACH(const face_pair_type &fp,
            bdry_face_group->face_pairs)
        {
          const hedge::straight_face &f = fp.int_side;
          m_el_face_to_face_pair_locator
            [std::make_pair(f.element_id, f.face_id)] =
            face_pair_locator(*bdry_face_group, fp);
        }

        m_activation_threshold = activation_threshold;
        m_kill_threshold = kill_threshold;
        m_upwind_alpha = upwind_alpha;
      }




      void add_local_diff_matrix(unsigned coordinate, const py_matrix &dmat)
      {
        if (coordinate != m_local_diff_matrices.size())
          throw std::runtime_error("local diff matrices added out of order");

        m_local_diff_matrices.push_back(dmat);
      }




      // convenience ----------------------------------------------------------
      unsigned get_dimensions_mesh() const
      { return m_mesh_data.m_dimensions; }




      // main driver ----------------------------------------------------------
      template<class Target>
      void deposit_densities_on_target(
          const depositor_state &ds,
          const ParticleState &ps,
          Target &tgt, boost::python::slice const &pslice) const
      {
        FOR_ALL_SLICE_INDICES(pslice, ps.particle_count)
        {
          FOR_ALL_SLICE_INDICES_INNER(particle_number, pn);

          tgt.begin_particle(pn);
          BOOST_FOREACH(const active_element &el,
              ds.m_advected_particles[pn].m_elements)
            tgt.add_shape_on_element(
                el.m_element_info->m_id,
                el.m_element_info->m_start,
                subrange(ds.m_rho, el.m_start_index, el.m_start_index+m_dofs_per_element));
          tgt.end_particle(pn);
        }
      }




      py_vector get_debug_quantity_on_mesh(
          depositor_state &ds,
          const particle_state &ps,
          const std::string &qty,
          py_vector const &velocities)
      {
        if (qty == "rhs")
          return map_particle_space_to_mesh_space(ds,
            get_advective_particle_rhs(ds, ps, velocities));
        if (qty == "active_elements")
          return get_active_elements(ds, ps);
        else if (qty == "fluxes")
          return map_particle_space_to_mesh_space(ds,
              calculate_fluxes(ds, ps, velocities));
        else if (qty == "minv_fluxes")
          return map_particle_space_to_mesh_space(ds,
              apply_elementwise_inverse_mass_matrix(ds,
              calculate_fluxes(ds, ps, velocities)));
        else if (qty == "local_div")
          return map_particle_space_to_mesh_space(ds,
              calculate_local_div(ds, ps, velocities));
        else
          throw std::runtime_error("invalid debug quantity");
      }




      void perform_depositor_upkeep(
          depositor_state &ds,
          ParticleState &ps)
      {
        // retire empty particle subelements
        if (m_kill_threshold == 0)
          throw std::runtime_error("zero kill threshold");

        particle_number pn = 0;
        BOOST_FOREACH(advected_particle &p, ds.m_advected_particles)
        {
          double particle_charge = fabs(ps.charges[pn]);
          for (unsigned i_el = 0; i_el < p.m_elements.size(); ++i_el)
          {
            active_element &el = p.m_elements[i_el];

            if (el.m_min_life)
              --el.m_min_life;

            const double element_charge = element_l1(
                el.m_element_info->m_jacobian,
                subrange(
                  ds.m_rho,
                  el.m_start_index,
                  el.m_start_index+m_dofs_per_element));

            if (el.m_min_life == 0  && element_charge / particle_charge < m_kill_threshold)
            {
              // retire this element
              const hedge::element_number_t en = el.m_element_info->m_id;

              // kill connections
              for (hedge::face_number_t fn = 0; fn < m_faces_per_element; ++fn)
              {
                hedge::element_number_t connected_en = el.m_connections[fn];
                if (connected_en != hedge::INVALID_ELEMENT)
                {
                  active_element &connected_el = *p.find_element(connected_en);

                  for (hedge::face_number_t cfn = 0; cfn < m_faces_per_element; ++cfn)
                  {
                    if (connected_el.m_connections[cfn] == en)
                      connected_el.m_connections[cfn] = hedge::INVALID_ELEMENT;
                  }
                }
              }

              deallocate_element(ds, el.m_start_index);

              // kill the element
              p.m_elements.erase(p.m_elements.begin()+i_el);
            }
            else
              ++i_el;
          }

          ++pn;
        }
      }




      void kill_advected_particle(depositor_state &ds, particle_number pn)
      {
        BOOST_FOREACH(active_element &el,
            ds.m_advected_particles[pn].m_elements)
          deallocate_element(ds, el.m_start_index);
      }




      void note_move(depositor_state &ds,
          particle_number from, particle_number to, unsigned size)
      {
        for (unsigned i = 0; i < size; ++i)
        {
          kill_advected_particle(ds, to+i);
          ds.m_advected_particles[to+i] = ds.m_advected_particles[from+i];
        }
      }




      void note_change_size(depositor_state &ds, unsigned particle_count)
      {
        for (particle_number pn = particle_count;
            pn < ds.m_advected_particles.size();
            ++pn)
          kill_advected_particle(ds, pn);

        ds.m_advected_particles.resize(particle_count);
      }




      // initialization -----------------------------------------------------
      void dump_particle(advected_particle const &p) const
      {
        std::cout << "particle, radius " << p.m_shape_function.radius() << std::endl;
        unsigned i_el = 0;
        BOOST_FOREACH(const active_element &el, p.m_elements)
        {
          std::cout << "#" << el.m_element_info->m_id << " cnx:(";
          for (unsigned fn = 0; fn < m_faces_per_element; ++fn)
            if (el.m_connections[fn] == hedge::INVALID_ELEMENT)
              std::cout << "X" << ',';
            else
              std::cout << el.m_connections[fn]  << ',';

          std::cout << ")" << std::endl;

          ++i_el;
        }
      }




      // vectors space administration ---------------------------------------
      /* Each element occupies a certain index range in the state
       * vector m_rho (as well as elsewhere). These functions perform
       * allocation and deallocation of space in these vectors.
       */

      /** Allocate a space for a new element in the state vector, return
       * the start index.
       */
      unsigned allocate_element(depositor_state &ds)
      {
        if (m_dofs_per_element == 0)
          throw std::runtime_error("tried to allocate element on uninitialized advection depositor");

        ds.m_element_activation_counter.tick();

        if (ds.m_freelist.size())
        {
          ++ds.m_active_elements;
          unsigned result = ds.m_freelist.back();
          ds.m_freelist.pop_back();
          return result*m_dofs_per_element;
        }

        // there are no gaps available.
        // return the past-end spot in the array, reallocate if necessary.
        unsigned avl_space = ds.m_rho.size() / m_dofs_per_element;

        if (ds.m_active_elements == avl_space)
        {
          ds.resize_rho(std::max(
                dyn_vector::size_type(m_dofs_per_element*1024),
                2*ds.m_rho.size()));
          if (ds.m_rho_dof_shift_listener.get())
            ds.m_rho_dof_shift_listener->note_change_size(ds.m_rho.size());
        }

        return (ds.m_active_elements++)*m_dofs_per_element;
      }


      void deallocate_element(depositor_state &ds, unsigned start_index)
      {
        if (start_index % m_dofs_per_element != 0)
          throw std::runtime_error("invalid advective element deallocation");

        const unsigned el_index = start_index/m_dofs_per_element;
        --ds.m_active_elements;

        ds.m_element_kill_counter.tick();

        // unless we're deallocating the last element, add it to the freelist.
        if (el_index != ds.m_active_elements+ds.m_freelist.size())
          ds.m_freelist.push_back(el_index);

        if (ds.m_rho_dof_shift_listener.get())
          ds.m_rho_dof_shift_listener->note_reset(start_index, m_dofs_per_element);
      }



      template <class VecType>
      py_vector map_particle_space_to_mesh_space(
          const depositor_state &ds,
          VecType const &pspace) const
      {
        py_vector result(
            m_dofs_per_element*m_mesh_data.m_element_info.size());
        result.clear();

        BOOST_FOREACH(const advected_particle &p, ds.m_advected_particles)
        {
          BOOST_FOREACH(const active_element &el, p.m_elements)
          {
            const mesh_data::element_info &einfo = *el.m_element_info;
            noalias(subrange(result, einfo.m_start, einfo.m_end)) +=
              subrange(pspace, el.m_start_index, el.m_start_index+m_dofs_per_element);
          }
        }

        return result;
      }




      py_vector get_active_elements(
          const depositor_state &ds,
          particle_state const &ps) const
      {
        py_vector result(
            m_dofs_per_element*m_mesh_data.m_element_info.size());
        result.clear();

        BOOST_FOREACH(const advected_particle &p, ds.m_advected_particles)
        {
          BOOST_FOREACH(const active_element &el, p.m_elements)
          {
            const mesh_data::element_info &einfo = *el.m_element_info;
            noalias(subrange(result, einfo.m_start, einfo.m_end)) +=
              boost::numeric::ublas::scalar_vector<double>(m_dofs_per_element, 1);
          }
        }

        return result;
      }




      // particle construction ----------------------------------------------
    private:
      struct advected_particle_element_target
      {
        private:
          advective_depositor &m_depositor;
          depositor_state &m_dep_state;
          advected_particle &m_particle;

        public:
          advected_particle_element_target(
              advective_depositor &dep,
              depositor_state &ds,
              advected_particle &p)
            : m_depositor(dep), m_dep_state(ds), m_particle(p)
          { }

          void add_shape_on_element(
              const bounded_vector &center,
              const mesh_data::element_number en
              )
          {
            const mesh_data::element_info &einfo(
                m_depositor.m_mesh_data.m_element_info[en]);

            active_element new_element;
            new_element.m_element_info = &einfo;
            unsigned start = new_element.m_start_index =
              m_depositor.allocate_element(m_dep_state);
            new_element.m_min_life = 0;

            for (unsigned i = 0; i < m_depositor.m_dofs_per_element; ++i)
              m_dep_state.m_rho[start+i] =
                m_particle.m_shape_function(
                    m_depositor.m_mesh_data.mesh_node(einfo.m_start+i)-center);

            m_particle.m_elements.push_back(new_element);
          }
      };




    public:
      void add_advective_particle(
          depositor_state &ds,
          const ParticleState &ps,
          shape_function sf, particle_number pn)
      {
        if (pn != ds.m_advected_particles.size())
          throw std::runtime_error("advected particle added out of sequence");

        advected_particle new_particle;
        new_particle.m_shape_function = sf;

        element_finder el_finder(m_mesh_data);

        advected_particle_element_target el_tgt(*this, ds, new_particle);
        el_finder(ps, el_tgt, pn, sf.radius());

        // make connections
        BOOST_FOREACH(active_element &el, new_particle.m_elements)
        {
          const mesh_data::element_info &einfo = *el.m_element_info;

          unsigned fn = 0;
          BOOST_FOREACH(const mesh_data::face_info &f, einfo.m_faces)
          {
            if (new_particle.find_element(f.m_neighbor))
              el.m_connections[fn] = f.m_neighbor;
            ++fn;
          }
        }

        // scale so the amount of charge is correct
        std::vector<double> unscaled_masses;
        BOOST_FOREACH(active_element &el, new_particle.m_elements)
          unscaled_masses.push_back(element_integral(
                el.m_element_info->m_jacobian,
                subrange(ds.m_rho,
                  el.m_start_index,
                  el.m_start_index+m_dofs_per_element)));

        const double charge = ps.charges[pn];
        const double total_unscaled_mass = std::accumulate(
            unscaled_masses.begin(), unscaled_masses.end(), double(0));

        double scale;
        if (total_unscaled_mass == 0)
        {
          WARN(boost::str(boost::format("deposited initial particle mass is zero"
                  "(particle %d, #elements=%d)") % pn % new_particle.m_elements.size()));
          scale = charge;
        }
        else
          scale = charge / total_unscaled_mass;

        BOOST_FOREACH(active_element &el, new_particle.m_elements)
          subrange(ds.m_rho,
              el.m_start_index,
              el.m_start_index+m_dofs_per_element)
          *= scale;

        ds.m_advected_particles.push_back(new_particle);
      }




      // rhs calculation ----------------------------------------------------
      py_vector calculate_local_div(
          depositor_state &ds,
          const ParticleState &ps,
          py_vector const &velocities) const
      {
        const unsigned dofs = ds.m_rho.size();
        const unsigned active_contiguous_elements =
          ds.m_active_elements + ds.m_freelist.size();

        py_vector local_div(dofs);
        local_div.clear();

        // calculate local rst derivatives ----------------------------------
        dyn_vector rst_derivs(get_dimensions_mesh()*dofs);
        rst_derivs.clear();
        using namespace boost::numeric::bindings;
        using blas::detail::gemm;

        for (unsigned loc_axis = 0; loc_axis < get_dimensions_mesh(); ++loc_axis)
        {
          const dyn_fortran_matrix &matrix = m_local_diff_matrices.at(loc_axis);

          gemm(
              'N', // "matrix" is row-major
              'N', // a contiguous array of vectors is column-major
              matrix.size1(),
              active_contiguous_elements,
              matrix.size2(),
              /*alpha*/ 1,
              /*a*/ traits::matrix_storage(matrix),
              /*lda*/ matrix.size2(),
              /*b*/ traits::vector_storage(ds.m_rho),
              /*ldb*/ m_dofs_per_element,
              /*beta*/ 1,
              /*c*/ traits::vector_storage(rst_derivs) + loc_axis*dofs,
              /*ldc*/ m_dofs_per_element
              );
        }

        // combine them into local part of dot(v, grad rho) -----------------
        {
          particle_number pn = 0;
          BOOST_FOREACH(const advected_particle &p, ds.m_advected_particles)
          {
            bounded_vector v = subrange(velocities,
                ps.vdim()*pn, ps.vdim()*(pn+1));

            BOOST_FOREACH(const active_element &el, p.m_elements)
            {
              for (unsigned loc_axis = 0; loc_axis < get_dimensions_mesh(); ++loc_axis)
              {
                double coeff = 0;
                for (unsigned glob_axis = 0; glob_axis < get_dimensions_mesh(); ++glob_axis)
                  coeff += -v[glob_axis] *
                    el.m_element_info->m_inverse_map.matrix()(loc_axis, glob_axis);

                subrange(local_div,
                    el.m_start_index,
                    el.m_start_index + m_dofs_per_element) += coeff *
                  subrange(rst_derivs,
                      loc_axis*dofs + el.m_start_index,
                      loc_axis*dofs + el.m_start_index + m_dofs_per_element);
              }
            }
            ++pn;
          }
        }

        return local_div;
      }




      py_vector calculate_fluxes(
          depositor_state &ds,
          const ParticleState &ps,
          py_vector const &velocities)
      {
        if (m_activation_threshold == 0)
          throw std::runtime_error("zero activation threshold");

        py_vector fluxes(ds.m_rho.size());
        fluxes.clear();

        particle_number pn = 0;
        BOOST_FOREACH(advected_particle &p, ds.m_advected_particles)
        {
          const double shape_peak =
            p.m_shape_function(
                boost::numeric::ublas::zero_vector<double>(get_dimensions_mesh()))
            * ps.charges[pn];

          const bounded_vector v = subrange(velocities,
              ps.vdim()*pn, ps.vdim()*(pn+1));

          for (unsigned i_el = 0; i_el < p.m_elements.size(); ++i_el)
          {
            active_element const *el = &p.m_elements[i_el];

            for (hedge::face_number_t fn = 0; fn < m_faces_per_element; ++fn)
            {
              const mesh_data::element_number en = el->m_element_info->m_id;

              /* Find correct face instance
               *
               * A face_pair represents both sides of a face. It points
               * to one or two hedge::face instances in its face_group that
               * carry information about each side of the face.
               *
               * The "ext" side of the face_pair may be unpopulated because
               * of a boundary.
               *
               * First, we need to identify which side of the face (en,fn)
               * identifies, guarding against an unpopulated "ext" side.
               */
              bool is_boundary;
              const face_pair_type::int_side_type *flux_face;
              const face_pair_type::int_side_type *external_flux_face;
              hedge::index_lists_t::const_iterator idx_list;
              hedge::index_lists_t::const_iterator ext_idx_list;

              {
                const face_pair_locator &fp_locator = map_get(
                    m_el_face_to_face_pair_locator,
                    std::make_pair(en, fn));

                const face_group_type &fg(*fp_locator.m_face_group);
                const face_pair_type &fp(*fp_locator.m_face_pair);

                const face_pair_type::int_side_type &flux_face_a = fp.int_side;

                const bool is_face_b = en != flux_face_a.element_id;

                is_boundary =
                  fp.ext_side.element_id == hedge::INVALID_ELEMENT;

                const face_pair_type::ext_side_type *flux_face_b = &fp.ext_side;

                if (is_boundary && is_face_b)
                  throw std::runtime_error("looking for non-existant cross-boundary element");

                if (is_face_b && en != flux_face_b->element_id)
                  throw std::runtime_error("el/face lookup failed");

                flux_face = is_face_b ? flux_face_b : &flux_face_a;
                external_flux_face = is_face_b ? &flux_face_a : flux_face_b;

                idx_list = fg.index_list(flux_face->face_index_list_number);
                ext_idx_list = fg.index_list(external_flux_face->face_index_list_number);
              }

              // Find information about this face
              const double n_dot_v = inner_prod(v, flux_face->normal);
              const bool inflow = n_dot_v <= 0;
              bool active = el->m_connections[fn] != mesh_data::INVALID_ELEMENT;

              if (is_boundary && active)
                throw std::runtime_error("detected boundary non-connection as active");

              const double int_coeff =
                flux_face->face_jacobian*(-n_dot_v)*(
                    m_upwind_alpha*(1 - (inflow ? 0 : 1))
                    +
                    (1-m_upwind_alpha)*0.5);
              const double ext_coeff =
                flux_face->face_jacobian*(-n_dot_v)*(
                    m_upwind_alpha*-(inflow ? 1 : 0)
                    +
                    (1-m_upwind_alpha)*-0.5);

              const mesh_data::node_number this_base_idx = el->m_start_index;

              // activate outflow, if necessary -----------------------------
              if (!is_boundary && !active && !inflow)
              {
                const unsigned face_length = m_face_mass_matrix.size1();

                double max_density = 0;
                for (unsigned i = 0; i < face_length; i++)
                  max_density = std::max(max_density,
                      fabs(ds.m_rho[this_base_idx+idx_list[i]]));

                // std::cout << max_density << ' ' << shape_peak << std::endl;
                if (max_density > m_activation_threshold*fabs(shape_peak))
                {
                  // yes, activate the external element

                  const hedge::element_number_t ext_en = external_flux_face->element_id;

                  const mesh_data::element_info &ext_einfo(
                      m_mesh_data.m_element_info[ext_en]);

                  active_element ext_element;
                  ext_element.m_element_info = &ext_einfo;

                  unsigned start = ext_element.m_start_index = allocate_element(ds);
                  subrange(ds.m_rho, start, start+m_dofs_per_element) =
                    boost::numeric::ublas::zero_vector<double>(m_dofs_per_element);

                  if (ds.m_rho.size() != fluxes.size())
                  {
                    // allocate_element enlarged the size of the state vector
                    // fluxes needs to be changed as well.
                    py_vector new_fluxes(ds.m_rho.size());
                    new_fluxes.clear();
                    noalias(subrange(new_fluxes, 0, fluxes.size())) = fluxes;
                    fluxes.swap(new_fluxes);
                  }

                  ext_element.m_min_life = 10;

                  // update connections
                  hedge::face_number_t ext_fn = 0;
                  BOOST_FOREACH(const mesh_data::face_info &ext_face, ext_einfo.m_faces)
                  {
                    const hedge::element_number_t ext_neigh_en = ext_face.m_neighbor;
                    active_element *ext_neigh_el = p.find_element(ext_neigh_en);
                    if (ext_neigh_el)
                    {
                      /* We found an active neighbor of our "external" element.
                       *
                       * Notation:
                       *        *
                       *       / \
                       *      /ext_neigh
                       *     *-----*
                       *    / \ext/
                       *   / el\ /
                       *  *-----*
                       *
                       * el: The element currently under consideration in the
                       *   "big" loop.
                       * ext: The "external" element that we just decided to
                       *   activate.
                       * ext_neigh: Neighbor of ext, also part of this
                       *   advected_particle
                       */

                       // First, tell ext that ext_neigh exists.
                      ext_element.m_connections[ext_fn] = ext_neigh_en;

                      // Next, tell ext_neigh that ext exists.
                      const mesh_data::element_info &ext_neigh_einfo(
                          m_mesh_data.m_element_info[ext_neigh_en]);

                      mesh_data::face_number ext_index_in_ext_neigh = 0;
                      for (;ext_index_in_ext_neigh < ext_neigh_einfo.m_faces.size()
                          ;++ext_index_in_ext_neigh)
                        if (ext_neigh_einfo.m_faces[ext_index_in_ext_neigh].m_neighbor
                            == ext_en)
                          break;

                      if (ext_index_in_ext_neigh == ext_neigh_einfo.m_faces.size())
                        throw std::runtime_error("ext not found in ext_neigh");

                      ext_neigh_el->m_connections[ext_index_in_ext_neigh] = ext_en;
                    }

                    ++ext_fn;
                  }

                  p.m_elements.push_back(ext_element);

                  // modification of m_elements might have invalidated el,
                  // refresh it

                  el = &p.m_elements[i_el];

                  active = true;
                }
              }

              // treat fluxes between active elements -----------------------
              if (active)
              {
                const active_element *ext_el = p.find_element(el->m_connections[fn]);
                const mesh_data::node_number ext_base_idx = ext_el->m_start_index;

                if (ext_el == 0)
                {
                  dump_particle(p);
                  throw std::runtime_error(
                      boost::str(boost::format("external element %d of (el:%d,face:%d) for active connection not found")
                      % el->m_connections[fn] % en % fn).c_str());
                }

                const unsigned face_length = m_face_mass_matrix.size1();

                for (unsigned i = 0; i < face_length; i++)
                {
                  const int ili = this_base_idx+idx_list[i];

                  hedge::index_lists_t::const_iterator ilj_iterator = idx_list;
                  hedge::index_lists_t::const_iterator oilj_iterator = ext_idx_list;

                  double res_ili_addition = 0;

                  for (unsigned j = 0; j < face_length; j++)
                  {
                    const double fmm_entry = m_face_mass_matrix(i, j);

                    const int ilj = this_base_idx+*ilj_iterator++;
                    const int oilj = ext_base_idx+*oilj_iterator++;

                    res_ili_addition +=
                      ds.m_rho[ilj]*int_coeff*fmm_entry
                      +ds.m_rho[oilj]*ext_coeff*fmm_entry;
                  }

                  fluxes[ili] += res_ili_addition;
                }
              }

              // handle zero inflow from inactive neighbors -----------------
              else if (inflow)
              {
                const unsigned face_length = m_face_mass_matrix.size1();

                for (unsigned i = 0; i < face_length; i++)
                {
                  const int ili = this_base_idx+idx_list[i];

                  hedge::index_lists_t::const_iterator ilj_iterator = idx_list;

                  double res_ili_addition = 0;

                  for (unsigned j = 0; j < face_length; j++)
                    res_ili_addition += ds.m_rho[this_base_idx+*ilj_iterator++]
                      *int_coeff
                      *m_face_mass_matrix(i, j);

                  fluxes[ili] += res_ili_addition;
                }
              }
            }

          }
          ++pn;
        }

        return fluxes;
      }




      py_vector apply_elementwise_inverse_mass_matrix(
          depositor_state &ds,
          py_vector const &operand) const
      {
        py_vector result(ds.m_rho.size());
        result.clear();

        const unsigned active_contiguous_elements =
          ds.m_active_elements + ds.m_freelist.size();

        using namespace boost::numeric::bindings;
        using blas::detail::gemm;

        const dyn_fortran_matrix &matrix = m_inverse_mass_matrix;
        gemm(
            'T', // "matrix" is row-major
            'N', // a contiguous array of vectors is column-major
            matrix.size1(),
            active_contiguous_elements,
            matrix.size2(),
            /*alpha*/ 1,
            /*a*/ traits::matrix_storage(matrix),
            /*lda*/ matrix.size2(),
            /*b*/ traits::vector_storage(operand),
            /*ldb*/ m_dofs_per_element,
            /*beta*/ 1,
            /*c*/ traits::vector_storage(result),
            /*ldc*/ m_dofs_per_element
            );

        // perform jacobian scaling
        BOOST_FOREACH(const advected_particle &p, ds.m_advected_particles)
          BOOST_FOREACH(const active_element &el, p.m_elements)
          {
            subrange(result,
                el.m_start_index,
                el.m_start_index+m_dofs_per_element) *=
            1/el.m_element_info->m_jacobian;
          }

        return result;
      }




      py_vector get_advective_particle_rhs(
          depositor_state &ds,
          const ParticleState &ps,
          py_vector const &velocities)
      {
        // calculate_fluxes may resize the state vector--calculate it first,
        // everything else later.
        py_vector fluxes = calculate_fluxes(ds, ps, velocities);

        return calculate_local_div(ds, ps, velocities)
        - apply_elementwise_inverse_mass_matrix(ds, fluxes);
      }




      depositor_state *apply_advective_particle_rhs(
          depositor_state &ds,
          const ParticleState &ps,
          py_vector const &rhs,
          boost::shared_ptr<number_shift_listener> rho_dof_sl)
      {
        if (m_filter_matrix.size1() && m_filter_matrix.size2())
        {
          using namespace boost::numeric::bindings;
          using blas::detail::gemm;

          const unsigned active_contiguous_elements =
            ds.m_active_elements + ds.m_freelist.size();

          const dyn_fortran_matrix &matrix = m_filter_matrix;

          dyn_vector new_rho(ds.m_rho.size());
          new_rho.size();

          gemm(
              'N',
              'N', // a contiguous array of vectors is column-major
              matrix.size1(),
              active_contiguous_elements,
              matrix.size2(),
              /*alpha*/ 1,
              /*a*/ traits::matrix_storage(matrix),
              /*lda*/ matrix.size2(),
              /*b*/ traits::vector_storage(rhs),
              /*ldb*/ m_dofs_per_element,
              /*beta*/ 1,
              /*c*/ traits::vector_storage(new_rho),
              /*ldc*/ m_dofs_per_element
              );
          return new depositor_state(ds, new_rho, rho_dof_sl);
        }
        else
          return new depositor_state(ds, ds.m_rho + rhs, rho_dof_sl);
      }




      template <class VectorExpression>
      double element_integral(double jacobian, const VectorExpression &ve)
      {
        return jacobian * inner_prod(m_integral_weights, ve);
      }




      template <class VectorExpression>
      double element_l1(double jacobian, const VectorExpression &ve)
      {
        return jacobian * inner_prod(m_integral_weights,
            pyublas::unary_op<pyublas::unary_ops::fabs>::apply(ve));
      }
  };
}




#endif
