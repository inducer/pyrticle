// Pyrticle - Particle in Cell in Python
// Main Module of the Python wrapper
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




#include <boost/python.hpp>
#include "pic_algorithm.hpp"
#include "wrap_pic.hpp"




using namespace pyrticle;
using namespace boost::python;




namespace
{
  template <class ParticleState>
  void expose_pic_basics()
  {
    {
      typedef ParticleState cl;
      class_<ParticleState>(
          ("ParticleState"+get_state_class_suffix<ParticleState>()).c_str())
        .SDEF_RW_MEMBER(particle_count)

        .SDEF_BYVAL_RW_MEMBER(containing_elements)
        .SDEF_BYVAL_RW_MEMBER(positions)
        .SDEF_BYVAL_RW_MEMBER(momenta)
        .SDEF_BYVAL_RW_MEMBER(charges)
        .SDEF_BYVAL_RW_MEMBER(masses)
        ;
    }

    def("get_velocities", get_velocities<ParticleState>);
    def("find_new_containing_element", find_new_containing_element<ParticleState>);
    def("update_containing_elements", update_containing_elements<ParticleState>);
  }
}




void expose_pic()
{
  EXPOSE_FOR_ALL_STATE_TYPES(expose_pic_basics, ());
}
