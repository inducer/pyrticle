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
    class type : public shape_element_finder<PICAlgorithm>
    {
      public:
        static const unsigned max_faces = 4;

        // FIXME This should be a safe assumption, but it is nonetheless just
        // that: an assumption.
        static const unsigned dimensions_mesh = PICAlgorithm::dimensions_pos;

        // member types -------------------------------------------------------
        struct active_element
        {
          const mesh_data::element_info *m_element_info;
          boost::array<mesh_data::element_number, type::max_faces> m_connections;
          unsigned m_start_index;

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
        };




        // member data --------------------------------------------------------
        unsigned                        m_faces_per_element;
        unsigned                        m_dofs_per_element;
        unsigned                        m_active_elements;
        std::vector<unsigned>           m_freelist;

        hedge::matrix                   m_mass_matrix;
        hedge::matrix                   m_inverse_mass_matrix;
        hedge::matrix                   m_face_mass_matrix;
        std::vector<hedge::matrix>      m_local_diff_matrices;
        boost::shared_ptr<hedge::face_group> m_face_group;

        boost::unordered_map<mesh_data::el_face, hedge::face_pair const *> m_el_face_to_face_pair;

        std::vector<advected_particle>  m_advected_particles;

        hedge::vector                   m_rho;





        // public interface ---------------------------------------------------
        type()
          : m_faces_per_element(0), m_dofs_per_element(0), m_active_elements(0)
        { }




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
          const unsigned dim = dimensions_mesh;
          particle_number pn = 0;
          BOOST_FOREACH(advected_particle &p, m_advected_particles)
          {
            hedge::vector v = subrange(velocities, 
                PICAlgorithm::dimensions_velocity*pn,
                PICAlgorithm::dimensions_velocity*(pn+1));

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
          BOOST_FOREACH(advected_particle &p, m_advected_particles)
          {
            BOOST_FOREACH(active_element &el, p.m_elements)
            {
              const mesh_data::element_info &einfo = *el.m_element_info;
              noalias(subrange(rho, einfo.m_start, einfo.m_end)) +=
                subrange(m_rho, el.m_start_index, el.m_start_index+m_dofs_per_element);
            }
          }
        }




        hedge::vector rhs_mesh_field(hedge::vector const &velocities)
        {
          hedge::vector result(m_dofs_per_element*PIC_THIS->m_mesh_data.m_element_info.size());
          result.clear();

          hedge::vector rhs = get_advective_particle_rhs(velocities);
          BOOST_FOREACH(advected_particle &p, m_advected_particles)
          {
            BOOST_FOREACH(active_element &el, p.m_elements)
            {
              const mesh_data::element_info &einfo = *el.m_element_info;
              noalias(subrange(result, einfo.m_start, einfo.m_end)) +=
                subrange(rhs, el.m_start_index, el.m_start_index+m_dofs_per_element);
            }
          }
          return result;
        }




        void perform_reconstructor_upkeep()
        {
          // retire empty particle subelements 
        }




        // initialization -----------------------------------------------------
        void setup_advective_reconstructor(
            unsigned faces_per_element, 
            unsigned dofs_per_element,
            const hedge::matrix &mass_matrix,
            const hedge::matrix &inverse_mass_matrix,
            const hedge::matrix &face_mass_matrix,
            boost::shared_ptr<hedge::face_group> fg
            )
        {
          m_faces_per_element = faces_per_element;
          m_dofs_per_element = dofs_per_element;
          resize_state(m_dofs_per_element * 1024);

          m_mass_matrix = mass_matrix;
          m_inverse_mass_matrix = inverse_mass_matrix;
          m_face_mass_matrix = face_mass_matrix;

          m_face_group = fg;

          // build m_el_face_to_face_pair
          BOOST_FOREACH(const hedge::face_pair &fp, fg->face_pairs)
          {
            hedge::fluxes::face &f = fg->flux_faces[fp.flux_face_index];
            m_el_face_to_face_pair[std::make_pair(f.element_id, f.face_id)] = &fp;

            hedge::fluxes::face &opp_f = fg->flux_faces[fp.opp_flux_face_index];
            m_el_face_to_face_pair[std::make_pair(opp_f.element_id, opp_f.face_id)] = &fp;
          }
        }




        void add_local_diff_matrix(unsigned coordinate, const hedge::matrix &dmat)
        {
          if (coordinate != m_local_diff_matrices.size())
            throw std::runtime_error("local diff matrices added out of order");

          m_local_diff_matrices.push_back(dmat);
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

          if (m_freelist.size())
          {
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

          // unless we're deallocating the last element, add it to the freelist.
          if (el_index != m_active_elements+m_freelist.size())
            m_freelist.push_back(el_index);
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

          shape_function sf(new_particle.m_radius, dimensions_mesh);

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
            for (int f = 0; f < m_faces_per_element; f++)
            {
              mesh_data::element_number connected_el = el.m_connections[f];
              if (new_particle.find_element(connected_el))
                el.m_connections[f] = connected_el;
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
        hedge::vector get_advective_particle_rhs(hedge::vector const &velocities)
        {
          const unsigned dofs = m_rho.size();
          const unsigned active_contiguous_elements = 
            m_active_elements + m_freelist.size();
          hedge::vector local_div(dofs);
          local_div.clear();

          // calculate local rst derivatives ----------------------------------
          hedge::vector rst_derivs(dimensions_mesh*dofs);
          using namespace boost::numeric::bindings;
          using blas::detail::gemm;

          for (unsigned loc_axis = 0; loc_axis < dimensions_mesh; ++loc_axis)
          {
            hedge::matrix &matrix = m_local_diff_matrices.at(loc_axis);

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
            BOOST_FOREACH(advected_particle &p, m_advected_particles)
            {
              hedge::vector v = subrange(velocities, 
                  PICAlgorithm::dimensions_velocity*pn,
                  PICAlgorithm::dimensions_velocity*(pn+1));

              BOOST_FOREACH(active_element &el, p.m_elements)
              {
                for (unsigned loc_axis = 0; loc_axis < dimensions_mesh; ++loc_axis)
                {
                  double coeff = 0;
                  for (unsigned glob_axis = 0; glob_axis < dimensions_mesh; ++glob_axis)
                    coeff += v[glob_axis] *
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

          // compute fluxes ---------------------------------------------------
          hedge::vector fluxes(dofs);
          fluxes.clear();
          {
            typedef hedge::vector::value_type scalar_t;
            particle_number pn = 0;
            BOOST_FOREACH(advected_particle &p, m_advected_particles)
            {
              const hedge::vector v = subrange(velocities, 
                  PICAlgorithm::dimensions_velocity*pn,
                  PICAlgorithm::dimensions_velocity*(pn+1));

              scalar_t norm_v = norm_2(v);

              // perform fluxes
              BOOST_FOREACH(active_element &el, p.m_elements)
              {
                /* Here, we look for active inbound and inactive outbound connections.
                 * Thanks to the chosen DG strong form and upwind flux, only these 
                 * require any action on our part.
                 *
                 * Consider:
                 * INFLOW: If there's an active neighbor past that inflow boundary,
                 * then of course we need to advect his data in--active inflow needs
                 * to be treated, hence. If we don't treat inactive inflow elements,
                 * we're pretending u+ == u-; this will probably be unstable, we'll 
                 * see.
                 *
                 * OUTLFOW: For a strong form DG with an upwind flux, and assuming that
                 * u+ == u- along outflow boundaries, outflow fluxes contribute
                 * nothing, and so can be ignored. Inactive outflow boundaries need 
                 * to be checked if we need to activate the past-interface element.
                 */
                for (mesh_data::face_number fn = 0; fn < m_faces_per_element; ++fn)
                {
                  const mesh_data::element_number en = el.m_element_info->m_id;

                  const hedge::face_pair &fp = *m_el_face_to_face_pair[
                    std::make_pair(en, fn)];
                  const hedge::fluxes::face &this_f = m_face_group->flux_faces[
                    fp.flux_face_index];
                  const hedge::fluxes::face &opp_f = m_face_group->flux_faces[
                    fp.opp_flux_face_index];

                  bool is_opposite = en != this_f.element_id;

                  if (is_opposite && en != opp_f.element_id)
                    throw std::runtime_error("el/face lookup failed");

                  const hedge::fluxes::face &flux_face(is_opposite ? opp_f : this_f);

                  const scalar_t n_dot_v = inner_prod(v, flux_face.normal);
                  const bool inflow = n_dot_v <= 0;
                  const bool active = el.m_connections[fn] != mesh_data::INVALID_ELEMENT;

                  if (active && inflow)
                  {
                    hedge::index_list &idx_list(is_opposite ?
                        m_face_group->index_lists[fp.opp_face_index_list_number]
                        : m_face_group->index_lists[fp.face_index_list_number]);
                    hedge::index_list &opp_idx_list(is_opposite ?
                        m_face_group->index_lists[fp.face_index_list_number]
                        : m_face_group->index_lists[fp.opp_face_index_list_number]);

                    const mesh_data::node_index this_base_idx = el.m_start_index;
                    const active_element *opp_el = p.find_element(el.m_connections[fn]);
                    const mesh_data::node_index opp_base_idx = opp_el->m_start_index;
                    
                    if (opp_el == 0)
                      throw std::runtime_error("opposite element for active connection not found");

                    const unsigned face_length = m_face_mass_matrix.size1();
                    assert(face_length == idx_list.size());
                    assert(face_length == opp_idx_list.size());

                    const double int_coeff = 
                      flux_face.face_jacobian*0.5*(n_dot_v + norm_v);
                    const double ext_coeff = 
                      flux_face.face_jacobian*0.5*-(n_dot_v + norm_v);

                    for (unsigned i = 0; i < face_length; i++)
                    {
                      const int ili = this_base_idx+idx_list[i];

                      hedge::index_list::const_iterator ilj_iterator = idx_list.begin();
                      hedge::index_list::const_iterator oilj_iterator = opp_idx_list.begin();

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
                  else if (!active && inflow)
                  {

                  }
                  else if (!active && !inflow)
                  {

                  }
                }
              }
            }
            ++pn;
          }

          // apply inverse mass matrix to fluxes ------------------------------
          hedge::vector minv_fluxes(dofs);
          {
            hedge::matrix &matrix = m_inverse_mass_matrix;
            gemm(
                'T', // "matrix" is row-major
                'N', // a contiguous array of vectors is column-major
                matrix.size1(),
                active_contiguous_elements,
                matrix.size2(),
                /*alpha*/ 1,
                /*a*/ traits::matrix_storage(matrix), 
                /*lda*/ matrix.size2(),
                /*b*/ traits::vector_storage(fluxes), 
                /*ldb*/ m_dofs_per_element,
                /*beta*/ 1,
                /*c*/ traits::vector_storage(minv_fluxes),
                /*ldc*/ m_dofs_per_element
                );

            // perform jacobian scaling
            BOOST_FOREACH(advected_particle &p, m_advected_particles)
              BOOST_FOREACH(active_element &el, p.m_elements)
                subrange(minv_fluxes, 
                    el.m_start_index, 
                    el.m_start_index+m_dofs_per_element) *= 1/el.m_element_info->m_jacobian;
          }

          // end result -------------------------------------------------------
          return local_div - minv_fluxes;
        }
        



        void apply_advective_particle_rhs(hedge::vector const &rhs)
        {
          m_rho += rhs;
        }




        template <class VectorExpression>
        double element_integral(double jacobian, const VectorExpression &ve)
        {
          return jacobian * inner_prod(
              boost::numeric::ublas::scalar_vector<typename VectorExpression::value_type>
              (m_mass_matrix.size1(), 1), prod(m_mass_matrix, ve));
        }
    };
  };
}




#endif
