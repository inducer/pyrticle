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




#ifndef _AFAVCZ_PYRTICLE_WRAP_PIC_HPP_INCLUDED
#define _AFAVCZ_PYRTICLE_WRAP_PIC_HPP_INCLUDED




#include <boost/lexical_cast.hpp>
#include "pic_algorithm.hpp"
#include "diagnostics.hpp"
#include "wrap_reconstructor.hpp"
#include "wrap_pusher.hpp"
#include "wrap_helpers.hpp"




namespace python = boost::python;
using namespace pyrticle;




namespace
{
  template <class PICAlgorithm>
  void expose_diagnostics()
  {
    python::def("kinetic_energies", kinetic_energies<PICAlgorithm>);
    python::def("total_charge", total_charge<PICAlgorithm>);
    python::def("particle_momentum", particle_momentum<PICAlgorithm>);
    python::def("particle_current", particle_current<PICAlgorithm>);
    python::def("rms_beam_size", rms_beam_size<PICAlgorithm>);
    python::def("rms_beam_emittance", rms_beam_emittance<PICAlgorithm>);
    python::def("rms_energy_spread", rms_energy_spread<PICAlgorithm>);
  }




  template <class PICAlgorithm>
  void expose_pic_algorithm()
  {
    std::string name = "PIC";
    name += PICAlgorithm::reconstructor::get_name();
    name += PICAlgorithm::particle_pusher::get_name();
    name += PICAlgorithm::element_finder::get_name();
    name += boost::lexical_cast<std::string>(PICAlgorithm::get_dimensions_pos());
    name += boost::lexical_cast<std::string>(PICAlgorithm::get_dimensions_velocity());

    using python::arg;

    typedef PICAlgorithm cl;

    python::class_<cl, boost::noncopyable> 
      pic_wrap(name.c_str(), python::init<unsigned, double>());

    pic_wrap
      .add_static_property("dimensions_pos", &cl::get_dimensions_pos)
      .add_static_property("dimensions_velocity", &cl::get_dimensions_velocity)

      .DEF_RO_MEMBER(mesh_data)

      .DEF_RW_MEMBER(particle_count)

      .add_property("containing_elements", 
          &cl::containing_elements, 
          &cl::set_containing_elements)
      .add_property("positions", &cl::positions, &cl::set_positions)
      .add_property("momenta", &cl::momenta, &cl::set_momenta)
      .add_property("charges", &cl::charges, &cl::set_charges)
      .add_property("masses", &cl::masses, &cl::set_masses)

      .DEF_RO_MEMBER(vacuum_c)

      .DEF_RO_MEMBER(find_same)
      .DEF_RO_MEMBER(find_by_neighbor)
      .DEF_RO_MEMBER(find_by_vertex)
      .DEF_RO_MEMBER(find_global)

      .DEF_RW_MEMBER(particle_number_shift_listener)
      .DEF_RW_MEMBER(boundary_hit_listener)

      .DEF_SIMPLE_METHOD(velocities)

      .DEF_SIMPLE_METHOD(move_particle)
      .DEF_SIMPLE_METHOD(note_change_particle_count)
      .DEF_SIMPLE_METHOD(kill_particle)

      .DEF_SIMPLE_METHOD(add_rhs)
      .DEF_SIMPLE_METHOD(find_new_containing_element)
      .DEF_SIMPLE_METHOD(update_containing_elements)

      ;

    if (PICAlgorithm::get_dimensions_velocity() == 3)
    {
      pic_wrap
        .def("forces", &cl::template forces< // full-field case
            py_vector, py_vector, py_vector,
            py_vector, py_vector, py_vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("bx"), arg("by"), arg("bz"),
             arg("velocities"), arg("verbose_vis")))
        ;
    }
    else if (PICAlgorithm::get_dimensions_velocity() == 2)
    {
      pic_wrap
        .def("forces", &cl::template forces< // TM case
            zero_vector, zero_vector, py_vector,
            py_vector, py_vector, zero_vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("bx"), arg("by"), arg("bz"),
             arg("velocities"), arg("verbose_vis")))
        .def("forces", &cl::template forces< // TE case
            py_vector, py_vector, zero_vector,
            zero_vector, zero_vector, py_vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("bx"), arg("by"), arg("bz"),
             arg("velocities"), arg("verbose_vis")))
        ;
    }

    pic_wrap
      .DEF_SIMPLE_METHOD(store_particle_vis_vector)
      .DEF_RW_MEMBER(vis_listener)
      ;

    expose_typed_pusher(
        pic_wrap, 
        (typename PICAlgorithm::particle_pusher *) 0
        );
    expose_typed_reconstructor(
        pic_wrap, 
        (typename PICAlgorithm::reconstructor *) 0
        );

    expose_diagnostics<PICAlgorithm>();
  }




  template <class Reconstructor, class Pusher>
  inline
  void expose_pic_all_dim()
  {

    expose_pic_algorithm<
        pic<
          pic_data<2,2>,
          Reconstructor,
          Pusher,
          face_based_element_finder
          >
        >();

    expose_pic_algorithm<
        pic<
          pic_data<3,3>,
          Reconstructor,
          Pusher,
          face_based_element_finder
          >
        >();

    /*
    expose_pic_algorithm<
        pic<
          pic_data<2,2>,
          Reconstructor,
          Pusher,
          heuristic_element_finder
          >
        >();
    expose_pic_algorithm<
        pic<
          pic_data<3,3>,
          Reconstructor,
          Pusher,
          heuristic_element_finder
          >
        >();
        */
  }




  template <class Reconstructor>
  inline
  void expose_pic_all_pushers_all_dim()
  {
    expose_pic_all_dim<Reconstructor, monomial_particle_pusher>();
    expose_pic_all_dim<Reconstructor, averaging_particle_pusher>();
  }




  template <class Reconstructor>
  inline
  void expose_pic_nontarget_pushers_all_dim()
  {
    expose_pic_all_dim<Reconstructor, monomial_particle_pusher>();
  }
}




#endif
