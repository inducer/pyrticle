// Pyrticle - Particle in Cell in Python
// Reconstruction based on advected shapes
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
#include <hedge/face_operators.hpp>
#include "tools.hpp"
#include "meshdata.hpp"
#include "rec_shape.hpp"




namespace pyrticle
{
  struct advective_reconstructor 
  {
    template <class PICAlgorithm>
    class type : 
      public shape_element_finder<PICAlgorithm>,
      public reconstructor_base
    {
      public:
        static const unsigned max_faces = 4;

        // member types -------------------------------------------------------
        struct active_element
        {
          const mesh_data::element_info *m_element_info;
          boost::array<mesh_data::element_number, type::max_faces> m_connections;
          unsigned m_start_index;
          unsigned m_min_life;

          active_element()
            : m_element_info(0)
          {
            for (unsigned i = 0; i < type::max_faces; i++)
              m_connections[i] = mesh_data::INVALID_ELEMENT;
          }
        };

        struct advected_particle
        {
          std::vector<active_element>   m_elements;
          double                        m_radius;

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




        // member data --------------------------------------------------------
        unsigned                        m_dimensions_mesh;

        unsigned                        m_faces_per_element;
        unsigned                        m_dofs_per_element;
        unsigned                        m_active_elements;
        std::vector<unsigned>           m_freelist;

        hedge::matrix                   m_mass_matrix;
        hedge::vector                   m_integral_weights;
        hedge::matrix                   m_inverse_mass_matrix;
        hedge::matrix                   m_face_mass_matrix;
        hedge::vector                   m_face_integral_weights;

        std::vector<hedge::matrix>      m_local_diff_matrices;

        boost::shared_ptr<hedge::face_group> m_int_face_group;
        boost::shared_ptr<hedge::face_group> m_bdry_face_group;

        struct face_pair_locator
        {
          hedge::face_group     const *m_face_group;
          hedge::face_pair      const *m_face_pair;

          face_pair_locator()
            : m_face_group(0), m_face_pair(0)
          { }
          face_pair_locator(
              hedge::face_group     const &face_group,
              hedge::face_pair      const &face_pair
              )
            : m_face_group(&face_group), m_face_pair(&face_pair)
          { }

        };

        boost::unordered_map<mesh_data::el_face, face_pair_locator> m_el_face_to_face_pair_locator;

        std::vector<advected_particle>  m_advected_particles;

        hedge::vector                   m_rho;

        boost::shared_ptr<number_shift_listener> m_rho_dof_shift_listener;

        event_counter m_element_activation_counter, m_element_kill_counter;

        double m_activation_threshold;
        double m_kill_threshold;
        double m_upwind_alpha;




        // public interface ---------------------------------------------------
        type()
          : m_faces_per_element(0), 
          m_dofs_per_element(0), m_active_elements(0),
          m_activation_threshold(0), m_kill_threshold(0)
        { }




        unsigned get_dimensions_mesh() const
        { return CONST_PIC_THIS->m_mesh_data.m_dimensions; }




        static const char *get_name()
        { return "Advective"; }




        void reconstruct_densities(
            hedge::vector &rho, 
            hedge::vector &j,
            const hedge::vector &velocities)
        {

        }




        void reconstruct_j(hedge::vector &j, const hedge::vector &velocities)
        {
          const unsigned dim = get_dimensions_mesh();
          particle_number pn = 0;
          BOOST_FOREACH(advected_particle &p, m_advected_particles)
          {
            hedge::vector v = subrange(velocities, 
                PICAlgorithm::get_dimensions_velocity()*pn,
                PICAlgorithm::get_dimensions_velocity()*(pn+1));

            BOOST_FOREACH(active_element &el, p.m_elements)
            {
              const mesh_data::element_info &einfo = *el.m_element_info;
              for (unsigned i = 0; i < dim; ++i)
              {
                noalias(subslice(j, 
                      dim*einfo.m_start, 
                      dim,
                      m_dofs_per_element)) += v[i] *
                  subrange(m_rho, 
                      el.m_start_index, 
                      el.m_start_index+m_dofs_per_element);

              }
            }
            ++pn;
          }
        }




        void reconstruct_rho(hedge::vector &rho)
        {
          rho = map_particle_space_to_mesh_space(m_rho);
        }




        hedge::vector get_debug_quantity_on_mesh(
            const std::string &qty, 
            hedge::vector const &velocities)
        {
          if (qty == "rhs")
            return map_particle_space_to_mesh_space(
              get_advective_particle_rhs(velocities));
          if (qty == "active_elements")
            return get_active_elements();
          else if (qty == "fluxes")
            return map_particle_space_to_mesh_space(
                calculate_fluxes(velocities));
          else if (qty == "minv_fluxes")
            return map_particle_space_to_mesh_space(
                apply_elementwise_inverse_mass_matrix(
                calculate_fluxes(velocities)));
          else if (qty == "local_div")
            return map_particle_space_to_mesh_space(
                calculate_local_div(velocities));
          else
            throw std::runtime_error("invalid debug quantity");
        }




        void perform_reconstructor_upkeep()
        {
          // retire empty particle subelements 
          if (m_kill_threshold == 0)
            throw std::runtime_error("zero kill threshold");

          particle_number pn = 0;
          BOOST_FOREACH(advected_particle &p, m_advected_particles)
          {
            double particle_charge = fabs(CONST_PIC_THIS->m_charges[pn]);
            for (unsigned i_el = 0; i_el < p.m_elements.size(); ++i_el)
            {
              active_element &el = p.m_elements[i_el];

              if (el.m_min_life)
                --el.m_min_life;

              const double element_charge = element_l1(
                  el.m_element_info->m_jacobian,
                  subrange(
                    m_rho,
                    el.m_start_index,
                    el.m_start_index+m_dofs_per_element));

              if (el.m_min_life == 0  && element_charge / particle_charge < m_kill_threshold)
              {
                // retire this element
                const hedge::element_number en = el.m_element_info->m_id;

                // kill connections
                for (hedge::face_number fn = 0; fn < m_faces_per_element; ++fn)
                {
                  hedge::element_number connected_en = el.m_connections[fn];
                  if (connected_en != hedge::INVALID_ELEMENT)
                  {
                    active_element &connected_el = *p.find_element(connected_en);

                    for (hedge::face_number cfn = 0; cfn < m_faces_per_element; ++cfn)
                    {
                      if (connected_el.m_connections[cfn] == en)
                        connected_el.m_connections[cfn] = hedge::INVALID_ELEMENT;
                    }
                  }
                }

                deallocate_element(el.m_start_index);

                // kill the element
                p.m_elements.erase(p.m_elements.begin()+i_el);
              }
              else
                ++i_el;
            }

            ++pn;
          }
        }




        void note_move(particle_number from, particle_number to, unsigned size)
        {
          for (unsigned i = 0; i < size; ++i)
          {
            BOOST_FOREACH(active_element &el, m_advected_particles[to+i].m_elements)
              deallocate_element(el.m_start_index);

            m_advected_particles[to+i] = m_advected_particles[from+i];

          }
        }




        void note_change_size(unsigned particle_count)
        {
          m_advected_particles.resize(particle_count);
        }




        // initialization -----------------------------------------------------
        void setup_advective_reconstructor(
            unsigned faces_per_element, 
            unsigned dofs_per_element,
            const hedge::matrix &mass_matrix,
            const hedge::matrix &inverse_mass_matrix,
            const hedge::matrix &face_mass_matrix,
            boost::shared_ptr<hedge::face_group> int_face_group,
            boost::shared_ptr<hedge::face_group> bdry_face_group,
            double activation_threshold,
            double kill_threshold,
            double upwind_alpha
            )
        {
          m_faces_per_element = faces_per_element;
          m_dofs_per_element = dofs_per_element;
          resize_state(m_dofs_per_element * 1024);

          m_mass_matrix = mass_matrix;
          m_integral_weights = prod(m_mass_matrix, 
              boost::numeric::ublas::scalar_vector<double>
              (m_mass_matrix.size1(), 1));
          m_inverse_mass_matrix = inverse_mass_matrix;

          m_face_mass_matrix = face_mass_matrix;
          m_face_integral_weights = prod(m_face_mass_matrix, 
              boost::numeric::ublas::scalar_vector<double>
              (m_face_mass_matrix.size1(), 1));

          m_int_face_group = int_face_group;
          m_bdry_face_group = bdry_face_group;

          // build m_el_face_to_face_pair_locator
          BOOST_FOREACH(const hedge::face_pair &fp, int_face_group->face_pairs)
          {
            hedge::fluxes::face &f = int_face_group->flux_faces[fp.flux_face_index];
            m_el_face_to_face_pair_locator
              [std::make_pair(f.element_id, f.face_id)] = 
              face_pair_locator(*int_face_group, fp);

            hedge::fluxes::face &opp_f = int_face_group->flux_faces[fp.opp_flux_face_index];
            m_el_face_to_face_pair_locator
              [std::make_pair(opp_f.element_id, opp_f.face_id)] = 
              face_pair_locator(*int_face_group, fp);
          }

          BOOST_FOREACH(const hedge::face_pair &fp, bdry_face_group->face_pairs)
          {
            hedge::fluxes::face &f = bdry_face_group->flux_faces[fp.flux_face_index];
            m_el_face_to_face_pair_locator
              [std::make_pair(f.element_id, f.face_id)] = 
              face_pair_locator(*bdry_face_group, fp);
          }

          m_activation_threshold = activation_threshold;
          m_kill_threshold = kill_threshold;
          m_upwind_alpha = upwind_alpha;
        }




        void add_local_diff_matrix(unsigned coordinate, const hedge::matrix &dmat)
        {
          if (coordinate != m_local_diff_matrices.size())
            throw std::runtime_error("local diff matrices added out of order");

          m_local_diff_matrices.push_back(dmat);
        }




          void dump_particle(advected_particle const &p) const
          {
            std::cout << "particle, radius " << p.m_radius << std::endl;
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
        /* Each element occupies a certain index range in the global state
         * vector m_rho (as well as elsewhere). These functions perform
         * allocation and deallocation of space in these vectors.
         */

        void resize_state(unsigned new_size)
        {
          unsigned old_size = m_rho.size();
          unsigned copy_size = std::min(new_size, old_size);

          hedge::vector new_rho(new_size);
          subrange(new_rho, 0, copy_size) = subrange(m_rho, 0, copy_size);
          new_rho.swap(m_rho);
        }

        /** Allocate a space for a new element in the state vector, return
         * the start index.
         */
        unsigned allocate_element()
        {
          if (m_dofs_per_element == 0)
            throw std::runtime_error("tried to allocate element on uninitialized advection reconstructor");

          m_element_activation_counter.tick();

          if (m_freelist.size())
          {
            ++m_active_elements;
            unsigned result = m_freelist.back();
            m_freelist.pop_back();
            return result*m_dofs_per_element;
          }

          // we're all full, no gaps available.
          // return the past-end spot in the array, reallocate if necessary.
          unsigned avl_space = m_rho.size() / m_dofs_per_element;

          if (m_active_elements == avl_space)
            resize_state(2*m_rho.size());

          return (m_active_elements++)*m_dofs_per_element;
        }

        void deallocate_element(unsigned start_index)
        {
          if (start_index % m_dofs_per_element != 0)
            throw std::runtime_error("invalid advective element deallocation");

          const unsigned el_index = start_index/m_dofs_per_element;
          --m_active_elements;

          m_element_kill_counter.tick();

          // unless we're deallocating the last element, add it to the freelist.
          if (el_index != m_active_elements+m_freelist.size())
            m_freelist.push_back(el_index);

          if (m_rho_dof_shift_listener.get())
            m_rho_dof_shift_listener->note_reset(start_index, m_dofs_per_element);
        }




        hedge::vector map_particle_space_to_mesh_space(hedge::vector const &pspace) const
        {
          hedge::vector result(
              m_dofs_per_element*CONST_PIC_THIS->m_mesh_data.m_element_info.size());
          result.clear();

          BOOST_FOREACH(const advected_particle &p, m_advected_particles)
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




        hedge::vector get_active_elements() const
        {
          hedge::vector result(
              m_dofs_per_element*CONST_PIC_THIS->m_mesh_data.m_element_info.size());
          result.clear();

          BOOST_FOREACH(const advected_particle &p, m_advected_particles)
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
        void add_shape_on_element(
            advected_particle &new_particle,
            const hedge::vector &center,
            const mesh_data::element_number en
            )
        {
          const mesh_data::element_info &einfo(
              CONST_PIC_THIS->m_mesh_data.m_element_info[en]);

          active_element new_element;
          new_element.m_element_info = &einfo;
          unsigned start = new_element.m_start_index = allocate_element();
          new_element.m_min_life = 0;

          shape_function sf(new_particle.m_radius, get_dimensions_mesh());

          for (unsigned i = 0; i < m_dofs_per_element; ++i)
            m_rho[start+i] =
              sf(CONST_PIC_THIS->m_mesh_data.m_nodes[einfo.m_start+i]-center);

          new_particle.m_elements.push_back(new_element);
        }




        void add_advective_particle(double radius, particle_number pn)
        {
          if (pn != m_advected_particles.size())
            throw std::runtime_error("advected particle added out of sequence");

          advected_particle new_particle;
          new_particle.m_radius = radius;

          add_shape(new_particle, pn, radius);

          // make connections
          BOOST_FOREACH(active_element &el, new_particle.m_elements)
          {
            const mesh_data::element_info &einfo = *el.m_element_info;

            unsigned fn = 0;
            BOOST_FOREACH(mesh_data::element_number neighbor, einfo.m_neighbors)
            {
              if (new_particle.find_element(neighbor))
                el.m_connections[fn] = neighbor;
              ++fn;
            }
          }

          // scale so the amount of charge is correct
          std::vector<double> unscaled_masses;
          BOOST_FOREACH(active_element &el, new_particle.m_elements)
            unscaled_masses.push_back(element_integral(
                  el.m_element_info->m_jacobian,
                  subrange(m_rho, 
                    el.m_start_index, 
                    el.m_start_index+m_dofs_per_element)));


          const double charge = CONST_PIC_THIS->m_charges[pn];
          const double total_unscaled_mass = std::accumulate(
              unscaled_masses.begin(), unscaled_masses.end(), double(0));

          if (total_unscaled_mass == 0)
            throw std::runtime_error("total reconstructed mass is zero");

          BOOST_FOREACH(active_element &el, new_particle.m_elements)
            subrange(m_rho, 
                el.m_start_index, 
                el.m_start_index+m_dofs_per_element) 
            *= charge / total_unscaled_mass;

          m_advected_particles.push_back(new_particle);
        }




        // rhs calculation ----------------------------------------------------
        hedge::vector calculate_local_div(hedge::vector const &velocities) const
        {
          const unsigned dofs = m_rho.size();
          const unsigned active_contiguous_elements = 
            m_active_elements + m_freelist.size();

          hedge::vector local_div(dofs);
          local_div.clear();

          // calculate local rst derivatives ----------------------------------
          hedge::vector rst_derivs(get_dimensions_mesh()*dofs);
          rst_derivs.clear();
          using namespace boost::numeric::bindings;
          using blas::detail::gemm;

          for (unsigned loc_axis = 0; loc_axis < get_dimensions_mesh(); ++loc_axis)
          {
            const hedge::matrix &matrix = m_local_diff_matrices.at(loc_axis);

            gemm(
                'T', // "matrix" is row-major
                'N', // a contiguous array of vectors is column-major
                matrix.size1(),
                active_contiguous_elements,
                matrix.size2(),
                /*alpha*/ 1,
                /*a*/ traits::matrix_storage(matrix), 
                /*lda*/ matrix.size2(),
                /*b*/ traits::vector_storage(m_rho), 
                /*ldb*/ m_dofs_per_element,
                /*beta*/ 1,
                /*c*/ traits::vector_storage(rst_derivs) + loc_axis*dofs, 
                /*ldc*/ m_dofs_per_element
                );
          }

          // combine them into local part of dot(v, grad rho) -----------------
          {
            particle_number pn = 0;
            BOOST_FOREACH(const advected_particle &p, m_advected_particles)
            {
              hedge::vector v = subrange(velocities, 
                  PICAlgorithm::dimensions_velocity*pn,
                  PICAlgorithm::dimensions_velocity*(pn+1));

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




        hedge::vector calculate_fluxes(hedge::vector const &velocities)
        {
          if (m_activation_threshold == 0)
            throw std::runtime_error("zero activation threshold");

          hedge::vector fluxes(m_rho.size());
          fluxes.clear();

          typedef hedge::vector::value_type scalar_t;
          particle_number pn = 0;
          BOOST_FOREACH(advected_particle &p, m_advected_particles)
          {
            const shape_function sf(p.m_radius, get_dimensions_mesh());
            const double shape_peak = sf(
                boost::numeric::ublas::zero_vector<double>(get_dimensions_mesh()))
              *CONST_PIC_THIS->m_charges[pn];

            const hedge::vector v = subrange(velocities, 
                PICAlgorithm::dimensions_velocity*pn,
                PICAlgorithm::dimensions_velocity*(pn+1));

            for (unsigned i_el = 0; i_el < p.m_elements.size(); ++i_el)
            {
              active_element const *el = &p.m_elements[i_el];

              for (hedge::face_number fn = 0; fn < m_faces_per_element; ++fn)
              {
                const mesh_data::element_number en = el->m_element_info->m_id;

                /* Find correct fluxes::face instance
                 *
                 * A face_pair represents both sides of a face. It points
                 * to one or two hedge::fluxes::face instances in its face_group that
                 * carry information about each side of the face.
                 *
                 * The "opp" side of the face_pair may be unpopulated because
                 * of a boundary.
                 *
                 * First, we need to identify which side of the face (en,fn)
                 * identifies, guarding against an unpopulated "opp" side.
                 */
                bool is_boundary;
                const hedge::fluxes::face *flux_face;
                const hedge::fluxes::face *opposite_flux_face;
                const hedge::index_list *idx_list;
                const hedge::index_list *opp_idx_list;

                {
                  const face_pair_locator &fp_locator = map_get(
                      m_el_face_to_face_pair_locator,
                      std::make_pair(en, fn));

                  const hedge::face_group &fg(*fp_locator.m_face_group);
                  const hedge::face_pair &fp(*fp_locator.m_face_pair);

                  const hedge::fluxes::face &flux_face_a = fg.flux_faces[
                    fp.flux_face_index];

                  const bool is_face_b = en != flux_face_a.element_id;

                  is_boundary = 
                    fp.opp_flux_face_index == hedge::face_pair::INVALID_INDEX;

                  const hedge::fluxes::face *flux_face_b = 0;
                  if (!is_boundary)
                    flux_face_b = &fg.flux_faces[fp.opp_flux_face_index];

                  if (is_boundary && is_face_b)
                    throw std::runtime_error("looking for non-existant cross-boundary element");

                  if (is_face_b && en != flux_face_b->element_id)
                    throw std::runtime_error("el/face lookup failed");

                  flux_face = is_face_b ? flux_face_b : &flux_face_a;
                  opposite_flux_face = is_face_b ? &flux_face_a : flux_face_b;

                  idx_list = is_face_b ?
                      &fg.index_lists[fp.opp_face_index_list_number]
                      : &fg.index_lists[fp.face_index_list_number];
                  opp_idx_list = is_face_b ?
                      &fg.index_lists[fp.face_index_list_number]
                      : &fg.index_lists[fp.opp_face_index_list_number];
                }

                // Find information about this face
                const scalar_t n_dot_v = inner_prod(v, flux_face->normal);
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

                const mesh_data::node_index this_base_idx = el->m_start_index;

                // activate outflow, if necessary -----------------------------
                if (!is_boundary && !active && !inflow)
                {
                  const unsigned face_length = m_face_mass_matrix.size1();
                  assert(face_length == idx_list->size());

                  scalar_t max_density = 0;
                  for (unsigned i = 0; i < face_length; i++)
                    max_density = std::max(max_density, 
                        fabs(m_rho[this_base_idx+(*idx_list)[i]]));

                  // std::cout << max_density << ' ' << shape_peak << std::endl;
                  if (max_density > m_activation_threshold*fabs(shape_peak))
                  {
                    // yes, activate the opposite element

                    const hedge::element_number opp_en = opposite_flux_face->element_id;

                    const mesh_data::element_info &opp_einfo(
                        CONST_PIC_THIS->m_mesh_data.m_element_info[opp_en]);

                    active_element opp_element;
                    opp_element.m_element_info = &opp_einfo;

                    unsigned start = opp_element.m_start_index = allocate_element();
                    subrange(m_rho, start, start+m_dofs_per_element) = 
                      boost::numeric::ublas::zero_vector<double>(m_dofs_per_element);

                    opp_element.m_min_life = 10;

                    // update connections
                    hedge::face_number opp_fn = 0;
                    BOOST_FOREACH(hedge::element_number opp_neigh_en, opp_einfo.m_neighbors)
                    {
                      active_element *opp_neigh_el = p.find_element(opp_neigh_en);
                      if (opp_neigh_el)
                      {
                        /* We found an active neighbor of our "opposite" element.
                         *
                         * Notation:
                         *        *
                         *       / \
                         *      /opp_neigh
                         *     *-----*
                         *    / \opp/
                         *   / el\ /
                         *  *-----*
                         *
                         * el: The element currently under consideration in the 
                         *   "big" loop.
                         * opp: The "opposite" element that we just decided to
                         *   activate.
                         * opp_neigh: Neighbor of opp, also part of this 
                         *   advected_particle
                         */

                         // First, tell opp that opp_neigh exists.
                        opp_element.m_connections[opp_fn] = opp_neigh_en;

                        // Next, tell opp_neigh that opp exists.
                        const mesh_data::element_info &opp_neigh_einfo(
                            CONST_PIC_THIS->m_mesh_data.m_element_info[opp_neigh_en]);
                        unsigned opp_index_in_opp_neigh = std::find(
                              opp_neigh_einfo.m_neighbors.begin(),
                              opp_neigh_einfo.m_neighbors.end(),
                              opp_en) - opp_neigh_einfo.m_neighbors.begin();

                        opp_neigh_el->m_connections[opp_index_in_opp_neigh] = opp_en;
                      }

                      ++opp_fn;
                    }

                    p.m_elements.push_back(opp_element);

                    // modification of m_elements might have invalidated el,
                    // refresh it

                    el = &p.m_elements[i_el];

                    active = true;
                  }
                }
                
                // treat fluxes between active elements -----------------------
                if (active)
                {
                  const active_element *opp_el = p.find_element(el->m_connections[fn]);
                  const mesh_data::node_index opp_base_idx = opp_el->m_start_index;

                  if (opp_el == 0)
                  {
                    dump_particle(p);
                    throw std::runtime_error(
                        str(boost::format("opposite element %d of (el:%d,face:%d) for active connection not found")
                        % el->m_connections[fn] % en % fn).c_str());
                  }

                  const unsigned face_length = m_face_mass_matrix.size1();
                  assert(face_length == idx_list->size());
                  assert(face_length == opp_idx_list->size());

                  for (unsigned i = 0; i < face_length; i++)
                  {
                    const int ili = this_base_idx+(*idx_list)[i];

                    hedge::index_list::const_iterator ilj_iterator = idx_list->begin();
                    hedge::index_list::const_iterator oilj_iterator = opp_idx_list->begin();

                    scalar_t res_ili_addition = 0;

                    for (unsigned j = 0; j < face_length; j++)
                    {
                      const scalar_t fmm_entry = m_face_mass_matrix(i, j);

                      const int ilj = this_base_idx+*ilj_iterator++;
                      const int oilj = opp_base_idx+*oilj_iterator++;

                      res_ili_addition += 
                        m_rho[ilj]*int_coeff*fmm_entry
                        +m_rho[oilj]*ext_coeff*fmm_entry;
                    }

                    fluxes[ili] += res_ili_addition;
                  }
                }

                // handle zero inflow from inactive neighbors -----------------
                else if (inflow)
                {
                  const unsigned face_length = m_face_mass_matrix.size1();
                  assert(face_length == idx_list->size());

                  for (unsigned i = 0; i < face_length; i++)
                  {
                    const int ili = this_base_idx+(*idx_list)[i];

                    hedge::index_list::const_iterator ilj_iterator = idx_list->begin();

                    scalar_t res_ili_addition = 0;

                    for (unsigned j = 0; j < face_length; j++)
                      res_ili_addition += m_rho[this_base_idx+*ilj_iterator++]
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




        hedge::vector apply_elementwise_inverse_mass_matrix(hedge::vector const &operand) const
        {
          hedge::vector result(m_rho.size());
          result.clear();

          const unsigned active_contiguous_elements = 
            m_active_elements + m_freelist.size();

          using namespace boost::numeric::bindings;
          using blas::detail::gemm;

          const hedge::matrix &matrix = m_inverse_mass_matrix;
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
          BOOST_FOREACH(const advected_particle &p, m_advected_particles)
            BOOST_FOREACH(const active_element &el, p.m_elements)
            {
              subrange(result, 
                  el.m_start_index, 
                  el.m_start_index+m_dofs_per_element) *= 
              1/el.m_element_info->m_jacobian;
            }

          return result;
        }




        hedge::vector get_advective_particle_rhs(hedge::vector const &velocities)
        {
          return calculate_local_div(velocities) 
          - apply_elementwise_inverse_mass_matrix(calculate_fluxes(velocities));
        }
        



        void apply_advective_particle_rhs(hedge::vector const &rhs)
        {
          m_rho += rhs;
        }




        template <class VectorExpression>
        double element_integral(double jacobian, const VectorExpression &ve)
        {
          return jacobian * inner_prod(m_integral_weights, ve);
        }




        template <class VectorExpression>
        double element_l1(double jacobian, const VectorExpression &ve)
        {
          hedge::vector u(ve.size());
          for (unsigned i = 0; i < ve.size(); ++i)
            u(i) = fabs(ve(i));
          return jacobian * inner_prod(m_integral_weights, u);
        }
    };
  };
}




#endif
