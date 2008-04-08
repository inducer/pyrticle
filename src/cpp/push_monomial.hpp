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




  class interpolator
  {
    private:
      unsigned m_el_start, m_el_end;
      py_vector m_interpolation_coefficients;

    public:
      interpolator(unsigned el_start, unsigned el_end, 
          const py_vector &intp_coeff)
        : m_el_start(el_start), m_el_end(el_end),
        m_interpolation_coefficients(intp_coeff)
      { }

      template <class VecType>
      const double operator()(const VecType &data) const
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
    py_matrix m_l_vandermonde_t;
    py_matrix m_u_vandermonde_t;
    csr_matrix m_p_vandermonde_t;
  };




  struct monomial_particle_pusher
  {
    template <class PICAlgorithm>
    class type : public pusher_base
    {
      public:
        static const char *get_name()
        { return "Monomial"; }

        std::vector<local_monomial_discretization> 
          m_local_discretizations;
        std::vector<unsigned> 
          m_ldis_indices;




        const interpolator make_interpolator(
            const bounded_vector &pt, 
            mesh_data::element_number in_element) const
        {
          const mesh_data::element_info &el_inf = 
            CONST_PIC_THIS->m_mesh_data.m_element_info[in_element];
          const local_monomial_discretization &ldis = 
            m_local_discretizations[m_ldis_indices[in_element]];

          unsigned basis_length = ldis.m_basis.size();
          dyn_vector mon_basis_values_at_pt(basis_length);

          dyn_vector unit_pt = el_inf.m_inverse_map(pt);

          for (unsigned i = 0; i < basis_length; i++)
            mon_basis_values_at_pt[i] = ldis.m_basis[i](unit_pt);

          dyn_vector permuted_basis_values = prod(
                    ldis.m_p_vandermonde_t,
                    mon_basis_values_at_pt);

          py_vector coeff = 
            solve(
                ldis.m_u_vandermonde_t,
                solve(
                  ldis.m_l_vandermonde_t,
                  permuted_basis_values,
                  boost::numeric::ublas::lower_tag()),
                boost::numeric::ublas::upper_tag());

          return interpolator(el_inf.m_start, el_inf.m_end, coeff);
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
            const py_vector &velocities,
            bool verbose_vis
            )
        {
          const unsigned xdim = CONST_PIC_THIS->get_dimensions_pos();
          const unsigned vdim = CONST_PIC_THIS->get_dimensions_velocity();

          py_vector result(CONST_PIC_THIS->m_particle_count * vdim);
          std::auto_ptr<py_vector> 
            vis_e, vis_b, vis_el_force, vis_mag_force;

          if (verbose_vis)
          {
            vis_e = std::auto_ptr<py_vector>(
                new py_vector(3*PIC_THIS->m_particle_count));
            vis_b = std::auto_ptr<py_vector>(
                new py_vector(3*PIC_THIS->m_particle_count));
            vis_el_force = std::auto_ptr<py_vector>(
                new py_vector(3*PIC_THIS->m_particle_count));
            vis_mag_force = std::auto_ptr<py_vector>(
                new py_vector(3*PIC_THIS->m_particle_count));
          }

          for (particle_number i = 0; i < PIC_THIS->m_particle_count; i++)
          {
            const unsigned v_pstart = vdim*i;
            const unsigned v_pend = vdim*(i+1);

            mesh_data::mesh_data::element_number in_el = 
              PIC_THIS->m_containing_elements[i];

            interpolator interp = make_interpolator(
                subrange(PIC_THIS->m_positions, xdim*i, xdim*(i+1)), in_el);

            bounded_vector e(3);
            e[0] = interp(ex);
            e[1] = interp(ey);
            e[2] = interp(ez);

            bounded_vector b(3);
            b[0] = interp(bx);
            b[1] = interp(by);
            b[2] = interp(bz);

            const double charge = PIC_THIS->m_charges[i];

            bounded_vector el_force(3);
            el_force[0] = charge*e[0];
            el_force[1] = charge*e[1];
            el_force[2] = charge*e[2];

            const bounded_vector v = subrange(velocities, v_pstart, v_pend);
            bounded_vector mag_force = cross(v, charge*b);

            // truncate forces to dimensions_velocity entries
            subrange(result, v_pstart, v_pend) = subrange(
                el_force + mag_force, 0, PIC_THIS->get_dimensions_velocity());

            if (verbose_vis)
            {
              subrange(*vis_e, 3*i, 3*(i+1)) = e;
              subrange(*vis_b, 3*i, 3*(i+1)) = b;
              subrange(*vis_el_force, 3*i, 3*(i+1)) = el_force;
              subrange(*vis_mag_force, 3*i, 3*(i+1)) = mag_force;
            }
          }

          if (verbose_vis)
          {
            PIC_THIS->store_particle_vis_vector("pt_e", *vis_e, 3);
            PIC_THIS->store_particle_vis_vector("pt_b", *vis_b, 3);
            PIC_THIS->store_particle_vis_vector("el_force", *vis_el_force, 3);
            PIC_THIS->store_particle_vis_vector("mag_force", *vis_mag_force, 3);
          }

          return result;
        }
    };
  };
}




#endif
