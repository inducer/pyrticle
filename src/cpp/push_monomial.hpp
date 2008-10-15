// Pyrticle - Particle in Cell in Python
// Particle pusher based on monomial interpolation
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





#ifndef _BADFJAH_PYRTICLE_PUSH_MONOMIAL_HPP_INCLUDED
#define _BADFJAH_PYRTICLE_PUSH_MONOMIAL_HPP_INCLUDED




#include <vector>
#include <boost/foreach.hpp>
#include <boost/assign/list_of.hpp> 
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector2.hpp>
#include "tools.hpp"
#include "bases.hpp"
#include "meshdata.hpp"
#include "pic_algorithm.hpp"




namespace pyrticle
{
  struct monomial_basis_function
  {
    /* exponents in each coordinate direction */
    std::vector<unsigned> m_exponents;

    monomial_basis_function(const std::vector<unsigned> &exponents)
      : m_exponents(exponents)
    { }

    monomial_basis_function(unsigned i, unsigned j)
    { m_exponents = boost::assign::list_of(i)(j); }

    monomial_basis_function(unsigned i, unsigned j, unsigned k)
    { m_exponents = boost::assign::list_of(i)(j)(k); }

    template <class VecType>
    const double operator()(const VecType &v) const
    {
      double result = 1;
      unsigned i = 0;
      BOOST_FOREACH(unsigned exp, m_exponents)
        result *= pow(v[i++], exp);

      return result;
    }
  };




  struct local_monomial_discretization
  {
    std::vector<monomial_basis_function> m_basis;
    py_fortran_matrix m_lu_vandermonde_t;
    py_int_vector m_lu_piv_vandermonde_t;
  };




  class interpolator
  {
    public:
      const local_monomial_discretization &m_ldis;
      dyn_vector m_interpolation_coefficients;

      interpolator(
          const local_monomial_discretization &ldis,
          unsigned particle_count)
        : m_ldis(ldis),
        m_interpolation_coefficients(ldis.m_basis.size()*particle_count)
      { }

      template <class VecType>
      const double operator()(
          const particle_number pn,
          const mesh_data::mesh_data::element_number en,
          const VecType &data) const
      {
        const unsigned basis_size = m_ldis.m_basis.size();
        return inner_prod(
            subrange(m_interpolation_coefficients,
              pn*basis_size, (pn+1)*basis_size),
            subrange(data, 
              en*basis_size, (en+1)*basis_size)
            );
      }

      const double operator()(
          const particle_number pn,
          const mesh_data::mesh_data::element_number en,
          zero_vector) const
      { return 0; }
  };




  template <class ParticleState>
  struct monomial_particle_pusher
  {
    public:
      typedef ParticleState particle_state;

      const mesh_data &m_mesh_data;
      std::vector<local_monomial_discretization> 
        m_local_discretizations;
      std::vector<unsigned> 
        m_ldis_indices;


      monomial_particle_pusher(const mesh_data &md)
        : m_mesh_data(md)
      { }

      interpolator make_interpolator(const ParticleState &ps) const
      {
        const unsigned xdim = ps.xdim();
        
        interpolator result(
            m_local_discretizations[0], ps.particle_count);

        for (particle_number pn = 0; pn < ps.particle_count; pn++)
        {
          mesh_data::mesh_data::element_number in_el = 
            ps.containing_elements[pn];
          const mesh_data::element_info &el_inf = 
            m_mesh_data.m_element_info[in_el];
        
          if (m_ldis_indices[in_el] != 0)
            throw std::runtime_error("more than one "
                "local discretization is currently not "
                "supported");

          bounded_vector unit_pt = el_inf.m_inverse_map
            .operator()<bounded_vector>(
                subrange(ps.positions, xdim*pn, xdim*(pn+1))
                );
          unsigned base_idx = result.m_ldis.m_basis.size()*pn;

          for (unsigned i = 0; i < result.m_ldis.m_basis.size(); i++)
            result.m_interpolation_coefficients[base_idx+i] 
              = result.m_ldis.m_basis[i](unit_pt);
        }

        {
          using namespace boost::numeric::bindings;
          
          const py_fortran_matrix &matrix = 
            result.m_ldis.m_lu_vandermonde_t;

          int info;
          lapack::detail::getrs(
              'N', 
              /*n*/ matrix.size1(),
              /*nrhs*/ ps.particle_count,
              traits::matrix_storage(matrix.as_ublas()),
              /*lda*/ matrix.size1(),
              traits::vector_storage(
                result.m_ldis.m_lu_piv_vandermonde_t),
              traits::vector_storage(
                result.m_interpolation_coefficients),
              /*ldb*/ matrix.size1(),
              &info);

          if (info < 0)
            throw std::runtime_error("invalid argument to getrs");
        }

        return result;
      }




      // why all these template arguments? In 2D and 1D,
      // instead of passing a hedge::vector, you may simply
      // pass a zero_vector, and interpolation will know to
      // not even compute anything, but just return zero.
      template <class EX, class EY, class EZ, 
               class BX, class BY, class BZ>
      py_vector forces(
          const EX &ex, const EY &ey, const EZ &ez,
          const BX &bx, const BY &by, const BZ &bz,
          ParticleState &ps,
          const py_vector &velocities,
          visualization_listener *vis_listener
          )
      {
        const unsigned vdim = ps.vdim();

        npy_intp res_dims[] = { ps.particle_count, vdim };
        py_vector result(2, res_dims);

        py_vector
          vis_e, vis_b, vis_el_force, vis_mag_force;

        if (vis_listener)
        {
          npy_intp dims[] = { ps.particle_count, 3 };
          vis_e = py_vector(2, dims);
          vis_b = py_vector(2, dims);
          vis_el_force = py_vector(2, dims);
          vis_mag_force = py_vector(2, dims);
        }

        interpolator interp = make_interpolator(ps);

        for (particle_number pn = 0; pn < ps.particle_count; pn++)
        {
          const unsigned v_pstart = vdim*pn;
          const unsigned v_pend = vdim*(pn+1);

          mesh_data::mesh_data::element_number in_el = ps.containing_elements[pn];

          bounded_vector e(3);
          e[0] = interp(pn, in_el, ex);
          e[1] = interp(pn, in_el, ey);
          e[2] = interp(pn, in_el, ez);

          bounded_vector b(3);
          b[0] = interp(pn, in_el, bx);
          b[1] = interp(pn, in_el, by);
          b[2] = interp(pn, in_el, bz);

          const double charge = ps.charges[pn];

          bounded_vector el_force(3);
          el_force[0] = charge*e[0];
          el_force[1] = charge*e[1];
          el_force[2] = charge*e[2];

          const bounded_vector v = subrange(velocities, v_pstart, v_pend);
          bounded_vector mag_force = cross(v, charge*b);

          // truncate forces to dimensions_velocity entries
          subrange(result, v_pstart, v_pend) = subrange(
              el_force + mag_force, 0, ps.vdim());

          if (vis_listener)
          {
            subrange(vis_e, 3*pn, 3*(pn+1)) = e;
            subrange(vis_b, 3*pn, 3*(pn+1)) = b;
            subrange(vis_el_force, 3*pn, 3*(pn+1)) = el_force;
            subrange(vis_mag_force, 3*pn, 3*(pn+1)) = mag_force;
          }
        }

        if (vis_listener)
        {
          vis_listener->store_particle_vis_vector("pt_e", vis_e);
          vis_listener->store_particle_vis_vector("pt_b", vis_b);
          vis_listener->store_particle_vis_vector("el_force", vis_el_force);
          vis_listener->store_particle_vis_vector("mag_force", vis_mag_force);
        }

        return result;
      }
  };
}




#endif
