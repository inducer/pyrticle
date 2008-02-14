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
      const hedge::vector &m_charges;
      const FX       &m_fx;
      const FY       &m_fy;
      const FZ       &m_fz;
      hedge::vector  *m_particlewise_field;
      hedge::vector  *m_field_variance;
      double         m_qfield_square_accumulator;

    public:
      static const unsigned field_components = 3;
      typedef boost::numeric::ublas::vector<double,
              boost::numeric::ublas::bounded_array<double, field_components> >
                field_vector_t;

      force_averaging_target(
          const hedge::vector &charges,
          const FX &fx, const FY &fy, const FZ &fz,
          hedge::vector *particlewise_field,
          hedge::vector *field_variance
          )
        : 
          m_current_particle(INVALID_PARTICLE), 
          m_charges(charges),
          m_fx(fx), m_fy(fy), m_fz(fz),
          m_particlewise_field(particlewise_field),
          m_field_variance(field_variance),
          m_qfield_square_accumulator(0)
      { 
        if (m_field_variance && !m_particlewise_field)
          throw std::runtime_error("cannot acquire field variance without acquiring mean field");
      }

      void begin_particle(particle_number pn)
      {
        m_current_particle = pn;
        m_qfield_square_accumulator = 0;
      }

      template <class RhoExpression>
      double weigh_field_component(
          const hedge::vector &field, 
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib)
      {
        return inner_prod(
              subrange(field, start_idx, start_idx+rho_contrib.size()),
              rho_contrib);
      }

      template <class RhoExpression>
      double weigh_field_component(
          const zero_vector &field, 
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib)
      {
        return 0;
      }

      template <class RhoExpression>
      field_vector_t get_qavg_field(
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib
          )
      {
        // Important: recall that these result in an average normalized
        // to the particle's charge.
        field_vector_t qfield(3);
        qfield[0] = weigh_field_component(m_fx, start_idx, rho_contrib);
        qfield[1] = weigh_field_component(m_fy, start_idx, rho_contrib);
        qfield[2] = weigh_field_component(m_fz, start_idx, rho_contrib);

        if (m_particlewise_field)
        {
          const unsigned fld_base = m_current_particle*3;
          for (unsigned i = 0; i < field_components; ++i)
            (*m_particlewise_field)[fld_base+i] += qfield[i];
        }
        if (m_field_variance)
        {
          for (unsigned i = 0; i < field_components; ++i)
            m_qfield_square_accumulator += qfield[i]*qfield[i];
        }

        return qfield;
      }

      void end_particle(particle_number pn)
      { 
        if (m_particlewise_field)
        {
          double charge = m_charges[pn];

          boost::numeric::ublas::vector_range<hedge::vector>
            particle_field(*m_particlewise_field, 
                boost::numeric::ublas::range(
                  pn*field_components, 
                  (pn+1)*field_components));

          particle_field /= charge;

          if (m_field_variance)
            (*m_field_variance)[pn] = 
              inner_prod(particle_field, particle_field)
              - m_qfield_square_accumulator/square(charge);
        }
      }

  };




  template <unsigned DimensionsVelocity, class FX, class FY, class FZ>
  class el_force_averaging_target : 
    public force_averaging_target<DimensionsVelocity, FX, FY, FZ>
  {
    private:
      typedef force_averaging_target<DimensionsVelocity, FX, FY, FZ> super;
      hedge::vector  &m_result;

    public:
      el_force_averaging_target(
          const hedge::vector &charges,
          const FX &fx, const FY &fy, const FZ &fz,
          hedge::vector *particlewise_field,
          hedge::vector *field_variance,
          hedge::vector &result
          )
        : 
          super(charges, fx, fy, fz, particlewise_field, field_variance),
          m_result(result)
      { }

      template <class RhoExpression>
      void add_shape_on_element(const mesh_data::element_number en, 
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib)
      {
        subrange(m_result, 
            this->m_current_particle*DimensionsVelocity,
            (1+this->m_current_particle)*DimensionsVelocity)
          += subrange(
              this->get_qavg_field(start_idx, rho_contrib),
              0, DimensionsVelocity);
      }
  };




  template <unsigned DimensionsVelocity, class FX, class FY, class FZ>
  class mag_force_averaging_target : 
    public force_averaging_target<DimensionsVelocity, FX, FY, FZ>
  {
    private:
      typedef force_averaging_target<DimensionsVelocity, FX, FY, FZ> super;
      const hedge::vector               &m_velocities;
      
      hedge::vector  &m_result;
      hedge::vector m_particle_qB_accumulator;

    public:
      mag_force_averaging_target(
          const hedge::vector &charges,
          const FX &fx, const FY &fy, const FZ &fz,
          const hedge::vector &velocities,
          hedge::vector *particlewise_field,
          hedge::vector *field_variance,
          hedge::vector &result
          )
        : 
          super(charges, fx, fy, fz, particlewise_field, field_variance), 
          m_velocities(velocities),
          m_result(result),
          m_particle_qB_accumulator(super::field_components)
      { }

      void begin_particle(particle_number pn)
      {
        super::begin_particle(pn);
        m_particle_qB_accumulator.clear();
      }

      template <class RhoExpression>
      void add_shape_on_element(const mesh_data::element_number en, 
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib)
      {
        m_particle_qB_accumulator +=
          this->get_qavg_field(start_idx, rho_contrib);
      }

      void end_particle(particle_number pn)
      {
        super::end_particle(pn);

        unsigned 
          pstart = pn*DimensionsVelocity,
          pend = (pn+1)*DimensionsVelocity;

        noalias(subrange(this->m_result, pstart, pend)) += 
          subrange(
              cross(
                subrange(m_velocities, pstart, pend), 
                m_particle_qB_accumulator),
              0, DimensionsVelocity);
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
        // pass a zero_vector, and the field averaging machinery
        // will statically know to not even compute anything, 
        // and just return zero.
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

          typedef el_force_averaging_target<
            PICAlgorithm::dimensions_velocity, EX, EY, EZ> el_tgt_t;
          typedef mag_force_averaging_target
            <PICAlgorithm::dimensions_velocity, BX, BY, BZ> mag_tgt_t;

          static const unsigned field_components = el_tgt_t::field_components;
          static const unsigned pcount = PIC_THIS->m_particle_count;


          hedge::vector el_force = 
            zero_vector(CONST_PIC_THIS->m_particle_count * vdim);
          hedge::vector mag_force = 
            zero_vector(CONST_PIC_THIS->m_particle_count * vdim);

          std::auto_ptr<hedge::vector> 
            vis_e, vis_b, vis_e_var, vis_b_var;

          if (verbose_vis)
          {
            vis_e = std::auto_ptr<hedge::vector>(
                new hedge::vector(field_components*pcount));
            vis_b = std::auto_ptr<hedge::vector>(
                new hedge::vector(field_components*pcount));
            vis_e_var = std::auto_ptr<hedge::vector>(
                new hedge::vector(pcount));
            vis_b_var = std::auto_ptr<hedge::vector>(
                new hedge::vector(pcount));
          }

          el_tgt_t el_tgt(CONST_PIC_THIS->m_charges,
              ex, ey, ez, vis_e.get(), vis_e_var.get(), el_force);
          mag_tgt_t mag_tgt(CONST_PIC_THIS->m_charges,
              bx, by, bz, velocities, vis_b.get(), vis_b_var.get(), mag_force);

          chained_reconstruction_target<el_tgt_t, mag_tgt_t> force_tgt(el_tgt, mag_tgt);
          CONST_PIC_THIS->reconstruct_densities_on_target(force_tgt);

          if (verbose_vis)
          {
            PIC_THIS->store_vis_vector("el_force", el_force);
            PIC_THIS->store_vis_vector("mag_force", mag_force);
            PIC_THIS->store_vis_vector("pt_e", *vis_e);
            PIC_THIS->store_vis_vector("pt_b", *vis_b);
            PIC_THIS->store_vis_vector("pt_e_var", *vis_e_var);
            PIC_THIS->store_vis_vector("pt_b_var", *vis_b_var);
          }

          return el_force+mag_force;
        }
    };
  };
}




#endif
