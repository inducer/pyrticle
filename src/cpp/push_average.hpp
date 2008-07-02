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
#include "pic_algorithm.hpp"




namespace pyrticle
{
  template <unsigned DimensionsVelocity, class FX, class FY, class FZ>
  class force_averaging_target
  {
    public:
      static const unsigned field_components = 3;
      typedef boost::numeric::ublas::vector<double,
              boost::numeric::ublas::bounded_array<double, field_components> >
                field_vector_t;

    protected:
      particle_number m_current_particle;
      const mesh_data &m_mesh_data;
      const dyn_vector &m_integral_weights;
      const FX       &m_fx;
      const FY       &m_fy;
      const FZ       &m_fz;
      py_vector      m_particlewise_field;
      py_vector      m_field_stddev;
      field_vector_t m_qfield_accumulator;
      double         m_qfield_square_accumulator;
      double         m_particle_charge;

    public:
      force_averaging_target(
          const mesh_data &md,
          const dyn_vector &integral_weights,
          const FX &fx, const FY &fy, const FZ &fz,
          py_vector particlewise_field,
          py_vector field_stddev
          )
        : 
          m_current_particle(INVALID_PARTICLE), 
          m_mesh_data(md),
          m_integral_weights(integral_weights),
          m_fx(fx), m_fy(fy), m_fz(fz),
          m_particlewise_field(particlewise_field),
          m_field_stddev(field_stddev),
          m_qfield_accumulator(field_components),
          m_qfield_square_accumulator(0),
          m_particle_charge(0)
      { 
        if (m_field_stddev.is_valid() && !m_particlewise_field.is_valid())
          throw std::runtime_error("cannot acquire field stddev without acquiring mean field");
      }

      void begin_particle(particle_number pn)
      {
        m_current_particle = pn;
        m_qfield_accumulator.clear();
        m_qfield_square_accumulator = 0;
        m_particle_charge = 0;
      }

      template <class RhoExpression>
      double weigh_field_component(
          const py_vector &field, 
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib)
      {
        return inner_prod(
              subrange(field, start_idx, start_idx+rho_contrib.size()),
              element_prod(rho_contrib, m_integral_weights));
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
      double weigh_field_component_squared(
          const py_vector &field, 
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib)
      {
        boost::numeric::ublas::vector_range<const py_vector>
          field_range(field, 
              boost::numeric::ublas::range(start_idx, start_idx+rho_contrib.size()));

        return inner_prod(
            element_prod(field_range, field_range),
            element_prod(rho_contrib, m_integral_weights));
      }

      template <class RhoExpression>
      double weigh_field_component_squared(
          const zero_vector &field, 
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib)
      {
        return 0;
      }

      template <class RhoExpression>
      void add_shape_on_element(
          const mesh_data::element_number en, 
          const mesh_data::node_number start_idx, 
          const RhoExpression &rho_contrib
          )
      {
        const double jacobian = m_mesh_data.m_element_info[en].m_jacobian;

        // Important: recall that these result in an average normalized
        // to the particle's charge.
        field_vector_t qfield(field_components);
        qfield[0] = jacobian*weigh_field_component(m_fx, start_idx, rho_contrib);
        qfield[1] = jacobian*weigh_field_component(m_fy, start_idx, rho_contrib);
        qfield[2] = jacobian*weigh_field_component(m_fz, start_idx, rho_contrib);

        m_qfield_accumulator += qfield;

        m_particle_charge += jacobian * inner_prod(rho_contrib, m_integral_weights);

        if (m_particlewise_field.is_valid())
        {
          subrange(m_particlewise_field, 
              m_current_particle*field_components, 
              (m_current_particle+1)*field_components) += qfield;
          if (m_field_stddev.is_valid())
          {
            m_qfield_square_accumulator += 
              jacobian*(
                  weigh_field_component_squared(m_fx, start_idx, rho_contrib)
                  + weigh_field_component_squared(m_fy, start_idx, rho_contrib)
                  + weigh_field_component_squared(m_fz, start_idx, rho_contrib));
          }
        }
      }

      void end_particle(particle_number pn)
      { 
        if (m_particlewise_field.is_valid())
        {
          boost::numeric::ublas::vector_range<py_vector>
            particle_field(m_particlewise_field, 
                boost::numeric::ublas::range(
                  pn*field_components, 
                  (pn+1)*field_components));

          if (m_particle_charge == 0)
          {
            WARN(str(boost::format(
                    "average pusher: particle %d had zero reconstructed charge") % pn
                  ));
            m_field_stddev[pn] = 0;
            return;
          }

          particle_field /= m_particle_charge;

          if (m_field_stddev.is_valid())
          {
            const double squared_mean = m_qfield_square_accumulator/m_particle_charge;
            const double mean_squared = inner_prod(particle_field, particle_field);

            double variance = squared_mean - mean_squared;
            if (variance < 0)
            {
              if (fabs(variance) > squared_mean * 1e-3)
              {
                std::cout 
                 << boost::format("squared_mean (%g) > mean_squared (%g) in calculating stddev, variance=%g")
                 % squared_mean % mean_squared % variance
                 << std::endl;
              }
              m_field_stddev[pn] = 0;
            }
            else
              m_field_stddev[pn] = sqrt(variance);
          }
        }
      }

  };




  template <unsigned DimensionsVelocity, class FX, class FY, class FZ>
  class el_force_averaging_target : 
    public force_averaging_target<DimensionsVelocity, FX, FY, FZ>
  {
    private:
      typedef force_averaging_target<DimensionsVelocity, FX, FY, FZ> super;
      stats_gatherer<double> *m_normalization_stats;
      py_vector::const_iterator m_charges;
      py_vector  &m_result;

    public:
      el_force_averaging_target(
          const mesh_data &md,
          const dyn_vector &integral_weights,
          const FX &fx, const FY &fy, const FZ &fz,
          py_vector particlewise_field,
          py_vector field_stddev,
          stats_gatherer<double> *normalization_stats,
          const py_vector &charges,
          py_vector &result
          )
        : 
          super(md, integral_weights, 
              fx, fy, fz, 
              particlewise_field, field_stddev),
          m_normalization_stats(normalization_stats),
          m_charges(charges.begin()),
          m_result(result)
      { }

      void end_particle(particle_number pn)
      {
        super::end_particle(pn);

        unsigned 
          pstart = pn*DimensionsVelocity,
          pend = (pn+1)*DimensionsVelocity;

        if (this->m_particle_charge == 0)
          return;

        const double scale = m_charges[pn]/this->m_particle_charge;

        if (m_normalization_stats)
          m_normalization_stats->add(scale);

        noalias(subrange(m_result, pstart, pend)) += 
          scale*subrange(
              this->m_qfield_accumulator,
              0, DimensionsVelocity);
      }
  };




  template <unsigned DimensionsVelocity, class FX, class FY, class FZ>
  class mag_force_averaging_target : 
    public force_averaging_target<DimensionsVelocity, FX, FY, FZ>
  {
    private:
      typedef force_averaging_target<DimensionsVelocity, FX, FY, FZ> super;
      const py_vector &m_velocities;
      
      stats_gatherer<double> *m_normalization_stats;
      py_vector::const_iterator m_charges;
      py_vector &m_result;

    public:
      mag_force_averaging_target(
          const mesh_data &md,
          const dyn_vector &integral_weights,
          const FX &fx, const FY &fy, const FZ &fz,
          const py_vector &velocities,
          py_vector particlewise_field,
          py_vector field_stddev,
          stats_gatherer<double> *normalization_stats,
          const py_vector &charges,
          py_vector &result
          )
        : 
          super(md, integral_weights, 
              fx, fy, fz, 
              particlewise_field, field_stddev), 
          m_velocities(velocities),
          m_normalization_stats(normalization_stats),
          m_charges(charges.begin()),
          m_result(result)
      { }

      void end_particle(particle_number pn)
      {
        super::end_particle(pn);

        unsigned 
          pstart = pn*DimensionsVelocity,
          pend = (pn+1)*DimensionsVelocity;

        if (this->m_particle_charge == 0)
          return;

        const double scale = m_charges[pn]/this->m_particle_charge;

        if (m_normalization_stats)
          m_normalization_stats->add(scale);

        noalias(subrange(m_result, pstart, pend)) += 
          subrange(
              scale*cross<bounded_vector>(
                subrange(m_velocities, pstart, pend), 
                this->m_qfield_accumulator),
              0, DimensionsVelocity);
      }
  };




  struct averaging_particle_pusher
  {
    template <class PICAlgorithm>
    class type : public pusher_base
    {
      private:
        dyn_vector   m_integral_weights;

      public:
        stats_gatherer<double> m_e_normalization_stats;
        stats_gatherer<double> m_b_normalization_stats;

        static const char *get_name()
        { return "Average"; }

        // initialization -----------------------------------------------------
        void setup_averaging_particle_pusher(const py_matrix &mass_matrix)
        {
          m_integral_weights = prod(mass_matrix, 
              boost::numeric::ublas::scalar_vector<double>
              (mass_matrix.size1(), 1));
        }




        // force calculation --------------------------------------------------
        // why all these template arguments? In 2D and 1D,
        // instead of passing a py_vector, you may simply
        // pass a zero_vector, and the field averaging machinery
        // will statically know to not even compute anything, 
        // and just return zero.
        template <class EX, class EY, class EZ, 
                 class BX, class BY, class BZ>
        py_vector forces(
            const EX &ex, const EY &ey, const EZ &ez,
            const BX &bx, const BY &by, const BZ &bz,
            const py_vector &velocities,
            bool verbose_vis
            )
        {
          const unsigned vdim = CONST_PIC_THIS->get_dimensions_velocity();

          typedef el_force_averaging_target<
            PICAlgorithm::dimensions_velocity, EX, EY, EZ> el_tgt_t;
          typedef mag_force_averaging_target
            <PICAlgorithm::dimensions_velocity, BX, BY, BZ> mag_tgt_t;

          const unsigned field_components = el_tgt_t::field_components;
          const unsigned pcount = PIC_THIS->m_particle_count;

          npy_intp res_dims[] = { PIC_THIS->m_particle_count, vdim };

          py_vector el_force(2, res_dims);
          py_vector mag_force(2, res_dims);
          el_force.clear();
          mag_force.clear();

          py_vector vis_e, vis_b, vis_e_stddev, vis_b_stddev;

          if (verbose_vis)
          {
            npy_intp dims[] = { pcount, field_components };
            vis_e = py_vector(2, dims);
            vis_b = py_vector(2, dims);
            vis_e.clear();
            vis_b.clear();

            vis_e_stddev = py_vector(zero_vector(pcount));
            vis_b_stddev = py_vector(zero_vector(pcount));
          }

          el_tgt_t el_tgt(CONST_PIC_THIS->m_mesh_data, m_integral_weights,
              ex, ey, ez, vis_e, vis_e_stddev, 
              &m_e_normalization_stats,
              CONST_PIC_THIS->m_charges, el_force);
          mag_tgt_t mag_tgt(CONST_PIC_THIS->m_mesh_data, m_integral_weights,
              bx, by, bz, velocities, vis_b, vis_b_stddev, 
              &m_b_normalization_stats,
              CONST_PIC_THIS->m_charges, mag_force);

          chained_reconstruction_target<el_tgt_t, mag_tgt_t> force_tgt(el_tgt, mag_tgt);
          CONST_PIC_THIS->reconstruct_densities_on_target(force_tgt,
              boost::python::slice());

          if (verbose_vis)
          {
            PIC_THIS->store_particle_vis_vector("el_force", el_force);
            PIC_THIS->store_particle_vis_vector("mag_force", mag_force);
            PIC_THIS->store_particle_vis_vector("pt_e", vis_e);
            PIC_THIS->store_particle_vis_vector("pt_b", vis_b);
            PIC_THIS->store_particle_vis_vector("pt_e_stddev", vis_e_stddev);
            PIC_THIS->store_particle_vis_vector("pt_b_stddev", vis_b_stddev);
          }

          return el_force+mag_force;
        }
    };
  };
}




#endif
