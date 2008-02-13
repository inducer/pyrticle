// Pyrticle - Particle in Cell in Python
// Particle pusher based on weighted averaging
// Copyright (C) 2008 Andreas Kloeckner
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




#ifndef _AFYAFDFYA_PYRTICLE_PUSH_AVERAGE_HPP_INCLUDED
#define _AFYAFDFYA_PYRTICLE_PUSH_AVERAGE_HPP_INCLUDED




#include "meshdata.hpp"
#include "bases.hpp"




namespace pyrticle
{
  template <unsigned DimensionsVelocity, class FX, class FY, class FZ>
  class force_averaging_target
  {
    protected:
      particle_number m_current_particle;
      const FX       &m_fx;
      const FY       &m_fy;
      const FZ       &m_fz;
      hedge::vector  &m_result;

    public:
      force_averaging_target(
          const FX &fx, const FY &fy, const FZ &fz,
          hedge::vector &result)
        : 
          m_current_particle(INVALID_PARTICLE), 
          m_fx(fx), m_fy(fy), m_fz(fz),
          m_result(result)
      { }

      void begin_particle(particle_number pn)
      {
        m_current_particle = pn;
      }

      template <class RhoExpression>
      double weigh_force_along_axis(
          const unsigned axis,
          const mesh_data::node_number start_idx, 
          const hedge::vector &field, 
          const RhoExpression &rho_contrib)
      {
        std::cout 
          << "WEIGH"
          << ' ' << start_idx 
          << ' ' << field.size()
          << ' ' << rho_contrib.size()
          << std::endl;
        return inner_prod(
              subrange(field, start_idx, rho_contrib.size()),
              rho_contrib);
      }

      template <class RhoExpression>
      double weigh_force_along_axis(
          const unsigned axis,
          const mesh_data::node_number start_idx, 
          const zero_vector &field, 
          const RhoExpression &rho_contrib)
      {
        return 0;
      }

      void end_particle(particle_number pn)
      { }
  };




  template <unsigned DimensionsVelocity, class FX, class FY, class FZ>
  class el_force_averaging_target : 
    public force_averaging_target<DimensionsVelocity, FX, FY, FZ>
  {
    private:
      typedef force_averaging_target<DimensionsVelocity, FX, FY, FZ> super;

    public:
      el_force_averaging_target(
          const FX &fx, const FY &fy, const FZ &fz,
          hedge::vector &result)
        : 
          super(fx, fy, fz, result)
      { }

      template <class RhoExpression>
      void add_shape_on_element(const mesh_data::element_number en, 
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib)
      {
        std::cout 
          << "EENTER " 
          << this->m_current_particle << std::endl;
        unsigned base_idx = this->m_current_particle*DimensionsVelocity;

        if (DimensionsVelocity >= 1)
          this->m_result[base_idx+0] += 
            this->weigh_force_along_axis(0, start_idx, this->m_fx, rho_contrib);
        if (DimensionsVelocity >= 2)
          this->m_result[base_idx+1] += 
            this->weigh_force_along_axis(1, start_idx, this->m_fy, rho_contrib);
        if (DimensionsVelocity >= 3)
          this->m_result[base_idx+2] += 
            this->weigh_force_along_axis(2, start_idx, this->m_fz, rho_contrib);
        std::cout << "ELEAVE" << std::endl;
      }
  };




  template <unsigned DimensionsVelocity, class FX, class FY, class FZ>
  class mag_force_averaging_target : 
    public force_averaging_target<DimensionsVelocity, FX, FY, FZ>
  {
    private:
      typedef force_averaging_target<DimensionsVelocity, FX, FY, FZ> super;
      const hedge::vector               &m_velocities;

    public:
      mag_force_averaging_target(
          const FX &fx, const FY &fy, const FZ &fz,
          const hedge::vector &velocities,
          hedge::vector &result
          )
        : 
          super(fx, fy, fz, result), m_velocities(velocities)
      { }

      template <class RhoExpression>
      void add_shape_on_element(const mesh_data::element_number en, 
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib)
      {
        std::cout << "MENTER" << std::endl;
        double qB[3] = {
          this->weigh_force_along_axis(0, start_idx, this->m_fx, rho_contrib),
          this->weigh_force_along_axis(1, start_idx, this->m_fy, rho_contrib),
          this->weigh_force_along_axis(2, start_idx, this->m_fz, rho_contrib)
        };

        unsigned 
          pstart = this->m_current_particle*DimensionsVelocity,
          pend = (this->m_current_particle+1)*DimensionsVelocity;

        noalias(subrange(this->m_result, pstart, pend)) += 
          cross(subrange(m_velocities, pstart, pend), qB);
        std::cout << "MLEAVE" << std::endl;
      }
  };




  struct averaging_particle_pusher
  {
    template <class PICAlgorithm>
    class type : public pusher_base
    {
      public:
        static const char *get_name()
        { return "Average"; }

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
          const unsigned vdim = CONST_PIC_THIS->get_dimensions_velocity();

          hedge::vector e_result = 
            zero_vector(CONST_PIC_THIS->m_particle_count * vdim);
          hedge::vector m_result = 
            zero_vector(CONST_PIC_THIS->m_particle_count * vdim);

          typedef el_force_averaging_target<
            PICAlgorithm::dimensions_velocity, EX, EY, EZ> el_tgt_t;
          typedef mag_force_averaging_target
            <PICAlgorithm::dimensions_velocity, BX, BY, BZ> mag_tgt_t;

          el_tgt_t el_tgt(ex, ey, ez, e_result);
          mag_tgt_t mag_tgt(bx, by, bz, velocities, m_result);

          chained_reconstruction_target<el_tgt_t, mag_tgt_t> force_tgt(el_tgt, mag_tgt);
          CONST_PIC_THIS->reconstruct_densities_on_target(force_tgt);

          if (verbose_vis)
          {
            PIC_THIS->store_vis_vector("el_force", e_result);
            PIC_THIS->store_vis_vector("mag_force", m_result);
          }

          return e_result+m_result;
        }
    };
  };
}




#endif
