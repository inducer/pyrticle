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
#include "diagnostics.hpp"




using namespace pyrticle;
using namespace boost::python;




namespace
{
  template <class ParticleState>
  void expose_diagnostics()
  {
    python::def("kinetic_energies", kinetic_energies<ParticleState>);
    python::def("total_charge", total_charge<ParticleState>);
    python::def("particle_momentum", particle_momentum<ParticleState>);
    python::def("particle_current", particle_current<ParticleState>);
    python::def("rms_beam_size", rms_beam_size<ParticleState>);
    python::def("rms_beam_emittance", rms_beam_emittance<ParticleState>);
    python::def("rms_energy_spread", rms_energy_spread<ParticleState>);
  }




  template <class ParticleState>
  void expose_pic_basics()
  {
    {
      typedef ParticleState cl;
      class_<ParticleState>(
          ("ParticleState"+get_state_class_suffix<ParticleState>()).c_str())
        .add_static_property("xdim", &cl::xdim)
        .add_static_property("vdim", &cl::vdim)

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
  EXPOSE_FOR_ALL_STATE_TYPES(expose_diagnostics, ());

  {
    typedef find_event_counters cl;
    class_<find_event_counters>("FindEventCounters")
      .SDEF_RW_MEMBER(find_same)
      .SDEF_RW_MEMBER(find_by_neighbor)
      .SDEF_RW_MEMBER(find_by_vertex)
      .SDEF_RW_MEMBER(find_global)
      ;
  }
}
