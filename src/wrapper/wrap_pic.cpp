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





#include "pic_algorithm.hpp"
#include "push_monomial.hpp"
#include "rec_shape.hpp"
#include "wrap_helpers.hpp"




namespace python = boost::python;
using namespace pyrticle;




namespace
{
  template <class PICAlgorithm>
  void expose_pic_algorithm(const std::string &name)
  {
    using python::arg;

    typedef PICAlgorithm cl;

    python::class_<cl, boost::noncopyable> 
      pic_wrap(name.c_str(), 
          python::init<mesh_data &, double>());
    pic_wrap
      .add_static_property("dimensions_pos", &cl::get_dimensions_pos)
      .add_static_property("dimensions_velocity", &cl::get_dimensions_velocity)

      .add_property("mesh_data", 
          python::make_function(
            &cl::get_mesh_data,
            python::return_internal_reference<>()
            ))

      .DEF_RO_MEMBER(containing_elements)
      .DEF_RW_MEMBER(positions)
      .DEF_RW_MEMBER(momenta)
      .DEF_RW_MEMBER(charges)
      .DEF_RW_MEMBER(masses)

      .DEF_RO_MEMBER(deadlist)

      // .DEF_RO_MEMBER(vis_info)

      /*
      .DEF_RO_MEMBER(same_searches)
      .DEF_RO_MEMBER(normal_searches)
      .DEF_RO_MEMBER(vertex_searches)
      .DEF_RO_MEMBER(global_searches)
      .DEF_RO_MEMBER(vertex_shape_adds)
      .DEF_RO_MEMBER(neighbor_shape_adds)
      .DEF_RO_MEMBER(periodic_hits)
      */

      .DEF_RO_MEMBER(vacuum_c)

      .DEF_SIMPLE_METHOD(velocities)
      .DEF_SIMPLE_METHOD(find_new_containing_element)
      .DEF_SIMPLE_METHOD(update_containing_elements)
      ;

    if (PICAlgorithm::get_dimensions_velocity() == 3)
    {
      pic_wrap
        .def("forces", &cl::template forces< // full-field case
            hedge::vector, hedge::vector, hedge::vector,
            hedge::vector, hedge::vector, hedge::vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("bx"), arg("by"), arg("bz"),
             arg("velocities"), arg("verbose_vis")))
        ;
    }
    else if (PICAlgorithm::get_dimensions_velocity() == 2)
    {
      pic_wrap
        .def("forces", &cl::template forces< // TM case
            zero_vector, zero_vector, hedge::vector,
            hedge::vector, hedge::vector, zero_vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("bx"), arg("by"), arg("bz"),
             arg("velocities"), arg("verbose_vis")))
        .def("forces", &cl::template forces< // TE case
            hedge::vector, hedge::vector, zero_vector,
            zero_vector, zero_vector, hedge::vector>,
            (arg("ex"), arg("ey"), arg("ez"), 
             arg("bx"), arg("by"), arg("bz"),
             arg("velocities"), arg("verbose_vis")))
        ;
    }

    pic_wrap
      .DEF_SIMPLE_METHOD(reconstruct_densities)
      .DEF_SIMPLE_METHOD(reconstruct_rho)
      .DEF_SIMPLE_METHOD(reconstruct_j)
      ;
  }
}




void expose_pic()
{
  expose_pic_algorithm<
      pic<
        pic_data<3,3>,
        shape_function_reconstructor,
        monomial_particle_pusher
        >
      >
      ("ParticleCloud33");
}
