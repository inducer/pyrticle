// Pyrticle - Particle in Cell in Python
// Python wrapper for PIC algorithm
// Copyright (C) 2007 Andreas Kloeckner
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or // (at your option) any later version.  // 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.





#include "wrap_pic.hpp"
#include "rec_shape.hpp"
#include "rec_target.hpp"




using namespace pyrticle;
using namespace boost::python;




namespace
{
  template <class Depositor>
  void expose_deposition_functions()
  {
    def("deposit_densities", reconstruct_densities<Depositor>);
    def("deposit_j", reconstruct_j<Depositor>);
    def("deposit_rho", reconstruct_rho<Depositor>);
  }




  template <class ParticleState>
  void expose_depositors_for_pstate()
  {
    {
      typedef shape_function_reconstructor<ParticleState, polynomial_shape_function> cl;
      class_<cl>(
        ("InterpolatingDepositor"+get_state_class_suffix<ParticleState>()).c_str(), 
        init<const mesh_data &>())
        .DEF_RW_MEMBER(shape_function)
        ;
    }
    
    EXPOSE_FOR_ALL_TARGET_RECONSTRUCTORS(expose_deposition_functions, ());
  }
}




void expose_deposition()
{
  EXPOSE_FOR_ALL_STATE_TYPES(expose_depositors_for_pstate, ());
}
