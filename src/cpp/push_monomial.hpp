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
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include "tools.hpp"
#include "meshdata.hpp"




namespace pyrticle
{
  struct monomial_basis_function
  {
    std::vector<unsigned> m_exponents;

    monomial_basis_function(const std::vector<unsigned> &exponents)
      : m_exponents(exponents)
    { }

    monomial_basis_function(unsigned i, unsigned j)
    { m_exponents = boost::assign::list_of(i)(j); }

    monomial_basis_function(unsigned i, unsigned j, unsigned k)
    { m_exponents = boost::assign::list_of(i)(j)(k); }

    const double operator()(const hedge::vector &v) const
    {
      double result = 1;
      unsigned i = 0;
      BOOST_FOREACH(unsigned exp, m_exponents)
        result *= pow(v[i++], exp);

      return result;
    }
  };




  class interpolator
  {
    private:
      unsigned m_el_start, m_el_end;
      hedge::vector m_interpolation_coefficients;

    public:
      interpolator(unsigned el_start, unsigned el_end, 
          const hedge::vector &intp_coeff)
        : m_el_start(el_start), m_el_end(el_end),
        m_interpolation_coefficients(intp_coeff)
      { }

      const double operator()(const hedge::vector &data) const
      {
        return inner_prod(m_interpolation_coefficients,
            subrange(data, m_el_start, m_el_end));
      }

      const double operator()(zero_vector) const
      { return 0; }
  };




  struct local_monomial_discretization
  {
    std::vector<monomial_basis_function> m_basis;
    hedge::matrix m_l_vandermonde_t;
    hedge::matrix m_u_vandermonde_t;
    csr_matrix m_p_vandermonde_t;
  };




  struct monomial_particle_pusher
  {
    template <class PICAlgorithm>
    class type
    {
      public:
        static const char *get_name()
        { return "Monomial"; }

        std::vector<local_monomial_discretization> 
          m_local_discretizations;
        std::vector<unsigned> 
          m_ldis_indices;

        const interpolator make_interpolator(
            const hedge::vector &pt, 
            mesh_data::element_number in_element) const
        {
          const mesh_data::element_info &el_inf = 
            CONST_PIC_THIS->m_mesh_data.m_element_info[in_element];
          const local_monomial_discretization &ldis = 
            m_local_discretizations[m_ldis_indices[in_element]];

          unsigned basis_length = ldis.m_basis.size();
          hedge::vector mon_basis_values_at_pt(basis_length);

          hedge::vector unit_pt = el_inf.m_inverse_map(pt);

          for (unsigned i = 0; i < basis_length; i++)
            mon_basis_values_at_pt[i] = ldis.m_basis[i](unit_pt);

          hedge::vector permuted_basis_values = prod(
                    ldis.m_p_vandermonde_t,
                    mon_basis_values_at_pt);

          hedge::vector coeff = 
            solve(
                ldis.m_u_vandermonde_t,
                solve(
                  ldis.m_l_vandermonde_t,
                  permuted_basis_values,
                  ublas::lower_tag()),
                ublas::upper_tag());

          return interpolator(el_inf.m_start, el_inf.m_end, coeff);
        }




        // why all these template arguments? In 2D and 1D,
        // instead of passing a hedge::vector, you may simply
        // pass a zero_vector, and interpolation will know to
        // not even compute anything, but just return zero.
        template <class EX, class EY, class EZ, 
                 class BX, class BY, class BZ>
        hedge::vector forces(
            const EX &ex, const EY &ey, const EZ &ez,
            const BX &bx, const BY &by, const BZ &bz,
            const hedge::vector &velocities,
            bool verbose_vis
            )
        {
          hedge::vector result(PIC_THIS->m_momenta.size());
          std::auto_ptr<hedge::vector> 
            vis_e, vis_b, vis_el_force, vis_lorentz_force;

          if (verbose_vis)
          {
            vis_e = std::auto_ptr<hedge::vector>(
                new hedge::vector(3*PIC_THIS->m_containing_elements.size()));
            vis_b = std::auto_ptr<hedge::vector>(
                new hedge::vector(3*PIC_THIS->m_containing_elements.size()));
            vis_el_force = std::auto_ptr<hedge::vector>(
                new hedge::vector(3*PIC_THIS->m_containing_elements.size()));
            vis_lorentz_force = std::auto_ptr<hedge::vector>(
                new hedge::vector(3*PIC_THIS->m_containing_elements.size()));
          }

          for (particle_number i = 0; i < PIC_THIS->m_containing_elements.size(); i++)
          {
            unsigned x_pstart = PIC_THIS->dimensions_pos*i;
            unsigned x_pend = PIC_THIS->dimensions_pos*(i+1);
            unsigned v_pstart = PIC_THIS->dimensions_velocity*i;
            unsigned v_pend = PIC_THIS->dimensions_velocity*(i+1);

            mesh_data::mesh_data::element_number in_el = 
              PIC_THIS->m_containing_elements[i];
            if (in_el == mesh_data::INVALID_ELEMENT)
            {
              subrange(result, v_pstart, v_pend) = zero_vector(
                  PIC_THIS->dimensions_pos);
              continue;
            }

            interpolator interp = make_interpolator(
                subrange(PIC_THIS->m_positions, x_pstart, x_pend), in_el);

            hedge::vector e(3);
            e[0] = interp(ex);
            e[1] = interp(ey);
            e[2] = interp(ez);

            hedge::vector b(3);
            b[0] = interp(bx);
            b[1] = interp(by);
            b[2] = interp(bz);

            const double charge = PIC_THIS->m_charges[i];

            hedge::vector el_force(3);
            el_force[0] = charge*e[0];
            el_force[1] = charge*e[1];
            el_force[2] = charge*e[2];

            const hedge::vector v = subrange(velocities, v_pstart, v_pend);
            hedge::vector lorentz_force = cross(v, charge*b);

            // truncate forces to dimensions_velocity entries
            subrange(result, v_pstart, v_pend) = subrange(
                el_force + lorentz_force, 0, PIC_THIS->dimensions_velocity);

            if (verbose_vis)
            {
              subrange(*vis_e, 3*i, 3*(i+1)) = e;
              subrange(*vis_b, 3*i, 3*(i+1)) = b;
              subrange(*vis_el_force, 3*i, 3*(i+1)) = el_force;
              subrange(*vis_lorentz_force, 3*i, 3*(i+1)) = lorentz_force;
            }
          }

          if (verbose_vis)
          {
            PIC_THIS->store_vis_vector("pt_e", *vis_e);
            PIC_THIS->store_vis_vector("pt_b", *vis_b);
            PIC_THIS->store_vis_vector("el_force", *vis_el_force);
            PIC_THIS->store_vis_vector("lorentz_force", *vis_lorentz_force);
          }

          return result;
        }
    };
  };
}




#endif
