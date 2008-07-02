// Pyrticle - Particle in Cell in Python
// Grid-based node finding reconstruction
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




#ifndef _ADFYUYHQHF_PYRTICLE_REC_GRID_FIND_HPP_INCLUDED
#define _ADFYUYHQHF_PYRTICLE_REC_GRID_FIND_HPP_INCLUDED




#include "grid.hpp"




namespace pyrticle
{
  struct grid_find_reconstructor : public grid_targets
  {
    template <class PICAlgorithm>
    class type : public grid_reconstructor_base<PICAlgorithm, brick>,
      public target_reconstructor_mixin<PICAlgorithm>
    {
      public:
        typedef std::vector<unsigned> node_number_list_starts_t;
        node_number_list_starts_t m_node_number_list_starts;
        typedef std::vector<mesh_data::node_number> node_number_lists_t;
        node_number_lists_t m_node_number_lists;




        template <class Target>
        bool reconstruct_particle_on_one_brick(Target tgt, 
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
                    - CONST_PIC_THIS->m_mesh_data.mesh_node(nn))
                  );
            }

            ++it;
          }
          
          return true;
        }




        // reconstruction entry points ----------------------------------------
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

          rho_target<py_vector> rho_tgt(rho);
          typedef j_target<PICAlgorithm::dimensions_velocity, py_vector, py_vector> 
            j_tgt_t;
          j_tgt_t j_tgt(j, velocities);
          chained_target<rho_target<py_vector>, j_tgt_t>
              tgt(rho_tgt, j_tgt);

          this->reconstruct_densities_on_grid_target(tgt, pslice);

          return boost::make_tuple(rho, j);
        }




        py_vector reconstruct_j(py_vector const &velocities,
            boost::python::slice const &pslice)
        {
          npy_intp dims[] = {
            CONST_PIC_THIS->m_mesh_data.node_count(),
            CONST_PIC_THIS->get_dimensions_velocity()
          };
          py_vector j(2, dims);

          j_target<PICAlgorithm::dimensions_velocity, py_vector, py_vector> 
            j_tgt(j, velocities);
          this->reconstruct_densities_on_grid_target(j_tgt, pslice);
          return j;
        }




        py_vector reconstruct_rho(boost::python::slice const &pslice)
        {
          py_vector rho(CONST_PIC_THIS->m_mesh_data.node_count());

          rho_target<py_vector> rho_tgt(rho);
          this->reconstruct_densities_on_grid_target(rho_tgt, pslice);
          return rho;
        }




        // public interface -----------------------------------------------------
        static const std::string get_name()
        { return "GridFind"; }
    };
  };
}




#endif
