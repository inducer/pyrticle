// Pyrticle - Particle in Cell in Python
// Grid-based reconstruction
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




#ifndef _AKFAHFHASHDFH_PYRTICLE_DEP_GRID_BASE_HPP_INCLUDED
#define _AKFAHFHASHDFH_PYRTICLE_DEP_GRID_BASE_HPP_INCLUDED




#include "bases.hpp"
#include "meshdata.hpp"
#include "tools.hpp"
#include "grid.hpp"




namespace pyrticle
{
  struct grid_targets
  {
    template <class VecType>
    class rho_target
    {
      private:
        typename VecType::iterator m_target;

      public:
        rho_target(VecType &target_vector)
          : m_target(target_vector.begin())
        { target_vector.clear(); }

        void begin_particle(particle_number pn)
        { }

        void add_shape_value(unsigned vec_idx, double q_shapeval)
        { m_target[vec_idx] += q_shapeval; }

        void end_particle(particle_number pn)
        { }
    };




    /** Reconstruction Target for the current density.
     */
    template<unsigned DimensionsVelocity, class JVecType, class VVecType>
    class j_target
    {
      private:
        typename JVecType::iterator m_target;
        typename VVecType::const_iterator m_velocities;
        double m_scale_factors[DimensionsVelocity];

      public:
        j_target(JVecType &target_vector, const VVecType &velocities)
          : m_target(target_vector.begin()), m_velocities(velocities.begin())
        { 
          target_vector.clear();
          for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
            m_scale_factors[axis] = 0;
        }

        void begin_particle(particle_number pn)
        {
          for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
            m_scale_factors[axis] = m_velocities[pn*DimensionsVelocity+axis];
        }

        void add_shape_value(const unsigned vec_idx, const double q_shapeval)
        { 
          unsigned const base = vec_idx*DimensionsVelocity;
          for (unsigned axis = 0; axis < DimensionsVelocity; axis++)
            m_target[base+axis] += m_scale_factors[axis]*q_shapeval;
        }

        void end_particle(particle_number pn)
        { }
    };





    template <class T1, class T2>
    class chained_target
    {
      private:
        T1 m_target1;
        T2 m_target2;

      public:
        chained_target(T1 &target1, T2 &target2)
          : m_target1(target1), m_target2(target2)
        { }

        void begin_particle(const particle_number pn)
        {
          m_target1.begin_particle(pn);
          m_target2.begin_particle(pn);
        }

        void add_shape_value(unsigned vec_idx, double q_shapeval)
        { 
          m_target1.add_shape_value(vec_idx, q_shapeval);
          m_target2.add_shape_value(vec_idx, q_shapeval);
        }

        void end_particle(const particle_number pn)
        {
          m_target1.end_particle(pn);
          m_target2.end_particle(pn);
        }
    };




    template <class T1, class T2>
    static 
    chained_target<T1, T2> 
    make_chained_target(T1 &target1, T2 &target2)
    {
      return chained_target<T1, T2>(target1, target2);
    }
  };
}




#endif
