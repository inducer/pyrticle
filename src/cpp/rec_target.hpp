// Pyrticle - Particle in Cell in Python
// Generic reconstruction target interface
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





#ifndef _ATHCNHF_PYRTICLE_REC_TARGET_HPP_INCLUDED
#define _ATHCNHF_PYRTICLE_REC_TARGET_HPP_INCLUDED




#include "meshdata.hpp"
#include "bases.hpp"




namespace pyrticle 
{
  /** The ReconstructionTarget protocol:
   *
   * template <class Scaler>
   * class reconstruction_target
   * {
   *   void begin_particle(particle_number pn);
   *
   *   template <class VectorExpression>
   *   void add_shape_on_element(mesh_data::element_number en, 
   *     mesh_data::node_number start_idx, 
   *     VectorExpression const &rho_contrib)
   *
   *   void end_particle(particle_number pn)
   * };
   *
   * Note: this is a stateful protocol.
   */

  class rho_reconstruction_target
  {
    private:
      py_vector &m_target_vector;

    public:
      rho_reconstruction_target(py_vector &target_vector)
        : m_target_vector(target_vector)
      { 
        m_target_vector.clear();
      }

      void begin_particle(particle_number pn)
      { }

      template <class VectorExpression>
      void add_shape_on_element(const mesh_data::element_number en, 
          const mesh_data::node_number start_idx, 
          VectorExpression const &rho_contrib)
      {
        noalias(subrange(m_target_vector, start_idx, start_idx+rho_contrib.size()))
          += rho_contrib;
      }

      void end_particle(particle_number pn)
      { }

      const py_vector &result() const
      {
        return m_target_vector;
      }
  };




  /** Reconstruction Target for the current density.
   */
  template<unsigned DimensionsVelocity>
  class j_reconstruction_target
  {
    private:
      py_vector &m_target_vector;
      const py_vector &m_velocities;
      double m_scale_factors[DimensionsVelocity];

    public:
      j_reconstruction_target(
          py_vector &target_vector, 
          const py_vector &velocities)
        : m_target_vector(target_vector), m_velocities(velocities)
      { 
        m_target_vector.clear();
        for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
          m_scale_factors[axis] = 0;
      }

      void begin_particle(particle_number pn)
      {
        for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
          m_scale_factors[axis] = m_velocities[pn*DimensionsVelocity+axis];
      }

      template <class VectorExpression>
      void add_shape_on_element(const mesh_data::element_number en, 
          const mesh_data::node_number start_idx, 
          VectorExpression const &rho_contrib)
      {
        for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
          noalias(subslice(m_target_vector, 
                start_idx*DimensionsVelocity+axis, 
                DimensionsVelocity, 
                rho_contrib.size()))
              += m_scale_factors[axis] * rho_contrib;
      }

      void end_particle(particle_number pn)
      { }

      const py_vector &result() const
      {
        return m_target_vector;
      }
  };





  template <class T1, class T2>
  class chained_reconstruction_target
  {
    private:
      T1 m_target1;
      T2 m_target2;

    public:
      chained_reconstruction_target(T1 &target1, T2 &target2)
        : m_target1(target1), m_target2(target2)
      { }

      void begin_particle(const particle_number pn)
      {
        m_target1.begin_particle(pn);
        m_target2.begin_particle(pn);
      }

      template <class VectorExpression>
      void add_shape_on_element(const mesh_data::element_number en, 
          const mesh_data::node_number start_idx, 
          VectorExpression const &rho_contrib)
      {
        m_target1.add_shape_on_element(en, start_idx, rho_contrib);
        m_target2.add_shape_on_element(en, start_idx, rho_contrib);
      }

      void end_particle(const particle_number pn)
      {
        m_target1.end_particle(pn);
        m_target2.end_particle(pn);
      }
  };




  template <class T1, class T2>
  inline
  chained_reconstruction_target<T1, T2> 
  make_chained_reconstruction_target(T1 &target1, T2 &target2)
  {
    return chained_reconstruction_target<T1, T2>(target1, target2);
  }




  // reconstructor base class -------------------------------------------------
  template <class PICAlgorithm>
  class target_reconstructor_mixin
  {
    public:
      boost::tuple<py_vector, py_vector> 
        reconstruct_densities(const py_vector &velocities, 
            boost::python::slice const &pslice)
      {
        py_vector rho(CONST_PIC_THIS->m_mesh_data.node_count());
        npy_intp dims[] = {
          CONST_PIC_THIS->m_mesh_data.node_count(),
          CONST_PIC_THIS->get_dimensions_velocity()
        };
        py_vector j(2, dims);

        rho_reconstruction_target rho_tgt(rho);
        typedef j_reconstruction_target<PICAlgorithm::dimensions_velocity> j_tgt_t;
        j_tgt_t j_tgt(j, velocities);

        chained_reconstruction_target<rho_reconstruction_target, j_tgt_t>
            tgt(rho_tgt, j_tgt);
        PIC_THIS->reconstruct_densities_on_target(tgt, pslice);

        rho = rho_tgt.result();

        return boost::make_tuple(rho, j);
      }




      py_vector reconstruct_j(const py_vector &velocities,
          boost::python::slice const &pslice)
      {
        npy_intp dims[] = {
          CONST_PIC_THIS->m_mesh_data.node_count(),
          CONST_PIC_THIS->get_dimensions_velocity()
        };

        py_vector j(2, dims);
        if (j.size() != PIC_THIS->m_mesh_data.node_count() *
            PIC_THIS->get_dimensions_velocity())
          throw std::runtime_error("j field does not have the correct size");

        j_reconstruction_target<PICAlgorithm::dimensions_velocity> j_tgt(
            j, velocities);

        PIC_THIS->reconstruct_densities_on_target(j_tgt, pslice);
        return j;
      }




      py_vector reconstruct_rho(boost::python::slice const &pslice)
      {
        py_vector rho(CONST_PIC_THIS->m_mesh_data.node_count());

        rho_reconstruction_target rho_tgt(rho);
        PIC_THIS->reconstruct_densities_on_target(rho_tgt, pslice);
        return rho;
      }
  };




  template <class PICAlgorithm>
  class target_reconstructor_base 
  : public reconstructor_base,
  public target_reconstructor_mixin<PICAlgorithm>
  { };
}




#endif
