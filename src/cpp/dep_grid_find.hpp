// Pyrticle - Particle in Cell in Python
// Grid-based node finding deposition
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




#ifndef _ADFYUYHQHF_PYRTICLE_DEP_GRID_FIND_HPP_INCLUDED
#define _ADFYUYHQHF_PYRTICLE_DEP_GRID_FIND_HPP_INCLUDED




#include "grid.hpp"
#include "dep_grid_base.hpp"




namespace pyrticle
{
  template <class ParticleState, class ShapeFunction, class Brick>
  class grid_find_depositor : 
    public grid_targets, 
    public grid_depositor_base<grid_find_depositor<ParticleState, ShapeFunction, Brick>,
      ParticleState, ShapeFunction, Brick, grid_depositor_base_state>
    // public target_depositor_mixin<PICAlgorithm>,
  {
    private:
      typedef grid_depositor_base<
        grid_find_depositor<ParticleState, ShapeFunction, Brick>, 
        ParticleState, ShapeFunction, Brick, grid_depositor_base_state>
          base_type;

    public:
      typedef typename base_type::particle_state particle_state;
      typedef typename base_type::brick_type brick_type;
      typedef typename base_type::depositor_state depositor_state;

      typedef std::vector<unsigned> node_number_list_starts_t;
      typedef std::vector<mesh_data::node_number> node_number_lists_t;

      // member data ----------------------------------------------------------
      node_number_list_starts_t m_node_number_list_starts;
      node_number_lists_t m_node_number_lists;




      grid_find_depositor(const mesh_data &md)
        : base_type(md)
      { }




      template <class Target>
      bool deposit_particle_on_one_brick(Target tgt, 
          const brick &brk, 
          const bounded_vector &center,
          const bounded_box &particle_box,
          double charge,
          bool *complete = 0) const
      {
        const bounded_box intersect_box = 
          brk.bounding_box().intersect(particle_box);

        if (complete)
          *complete = intersect_box == particle_box;

        if (intersect_box.is_empty())
          return false;

        const bounded_int_box particle_brick_index_box = 
          brk.index_range(intersect_box);

        brick::iterator it(brk, particle_brick_index_box);

        while (!it.at_end())
        {
          BOOST_FOREACH(mesh_data::node_number nn, std::make_pair(
            m_node_number_lists.begin() 
            + m_node_number_list_starts[it.index()],
            m_node_number_lists.begin() 
            + m_node_number_list_starts[it.index()+1])
              )
          {
            tgt.add_shape_value(nn, 
                charge*this->m_shape_function(center
                  - this->m_mesh_data.mesh_node(nn))
                );
          }

          ++it;
        }
        
        return true;
      }




      // deposition entry points ----------------------------------------
      boost::tuple<py_vector, py_vector> 
      deposit_densities(
          depositor_state &ds,
          const particle_state &ps,
          const py_vector &velocities,
          boost::python::slice const &pslice)
      {
        py_vector rho(this->m_mesh_data.node_count());
        npy_intp dims[] = {
          this->m_mesh_data.node_count(),
          particle_state::vdim()
        };
        py_vector j(2, dims);

        rho_target<py_vector> rho_tgt(rho);
        typedef j_target<particle_state::m_vdim, py_vector, py_vector> 
          j_tgt_t;
        j_tgt_t j_tgt(j, velocities);
        chained_target<rho_target<py_vector>, j_tgt_t>
            tgt(rho_tgt, j_tgt);

        this->deposit_densities_on_grid_target(ds, ps, tgt, pslice);

        return boost::make_tuple(rho, j);
      }




      py_vector deposit_j(
          depositor_state &ds,
          const particle_state &ps,
          py_vector const &velocities,
          boost::python::slice const &pslice)
      {
        npy_intp dims[] = {
          this->m_mesh_data.node_count(),
          particle_state::vdim()
        };
        py_vector j(2, dims);

        j_target<particle_state::m_vdim, py_vector, py_vector> 
          j_tgt(j, velocities);
        this->deposit_densities_on_grid_target(ds, ps, j_tgt, pslice);
        return j;
      }




      py_vector deposit_rho(
          depositor_state &ds,
          const particle_state &ps,
          boost::python::slice const &pslice)
      {
        py_vector rho(this->m_mesh_data.node_count());

        rho_target<py_vector> rho_tgt(rho);
        this->deposit_densities_on_grid_target(ds, ps, rho_tgt, pslice);
        return rho;
      }
  };
}




#endif
