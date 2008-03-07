// Pyrticle - Particle in Cell in Python
// Python wrapper for little bits of usefulness
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




#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include "tools.hpp"
#include "wrap_helpers.hpp"




namespace python = boost::python;
using namespace pyrticle;




namespace
{
  struct visualization_listener_wrap : 
    visualization_listener,
    python::wrapper<visualization_listener>
  {
    void store_particle_vis_vector(
        const char *name,
        const hedge::vector &vec,
        unsigned entries_per_particle) const
    {
      this->get_override("store_particle_vis_vector")(name, vec, entries_per_particle);
    }
  };




  struct number_shift_listener_wrap : 
    number_shift_listener,
    python::wrapper<number_shift_listener>
  {
    void note_change_size(unsigned new_size) const
    {
      if (python::override f = this->get_override("note_change_size"))
        f(new_size);
      else
        number_shift_listener::note_change_size(new_size);
    }

    void note_move(unsigned orig, unsigned dest, unsigned size) const
    {
      if (python::override f = this->get_override("note_move"))
        f(orig, dest, size);
      else
        number_shift_listener::note_move(orig, dest, size);
    }

    void note_reset(unsigned start, unsigned size) const
    {
      if (python::override f = this->get_override("note_reset"))
        f(start, size);
      else
        number_shift_listener::note_reset(start, size);
    }
  };
}




void expose_tools()
{
  python::def("asinh", (double (*)(double)) boost::math::asinh);
  python::def("acosh", (double (*)(double)) boost::math::acosh);
  python::def("gamma", (double (*)(double)) boost::math::tgamma);

  {
    typedef event_counter cl;
    python::class_<cl>("EventCounter")
      .DEF_SIMPLE_METHOD(get)
      .DEF_SIMPLE_METHOD(pop)
      .DEF_SIMPLE_METHOD(tick)
      ;
  }

  {
    typedef zero_vector cl;
    python::class_<cl>("ZeroVector");
  }

  {
    typedef std::vector<unsigned> cl;
    python::class_<cl>("UnsignedVector")
      .def(python::vector_indexing_suite<cl>())
      .DEF_SIMPLE_METHOD(clear)
      .DEF_SIMPLE_METHOD(reserve)
      ;
  }

  {
    typedef visualization_listener cl;
    python::class_<visualization_listener_wrap, 
      boost::shared_ptr<visualization_listener_wrap>,
      boost::noncopyable>
      ("VisualizationListener")
      .DEF_PURE_VIRTUAL_METHOD(store_particle_vis_vector)
      ;
  }

  {
    typedef number_shift_listener cl;
    typedef number_shift_listener_wrap wrp;
    python::class_<number_shift_listener_wrap, 
      boost::shared_ptr<number_shift_listener_wrap>,
      boost::noncopyable>
      ("NumberShiftListener")
      .def("note_change_size", &cl::note_change_size, &wrp::note_change_size)
      .def("note_move", &cl::note_move, &wrp::note_move)
      .def("note_reset", &cl::note_reset, &wrp::note_reset)
      ;
  }

  {
    typedef stats_gatherer<double> cl;
    python::class_<cl>("StatsGatherer")
      .DEF_SIMPLE_METHOD(add)
      .DEF_SIMPLE_METHOD(reset)
      .DEF_SIMPLE_METHOD(count)
      .DEF_SIMPLE_METHOD(minimum)
      .DEF_SIMPLE_METHOD(maximum)
      .DEF_SIMPLE_METHOD(mean)
      .DEF_SIMPLE_METHOD(variance)
      .DEF_SIMPLE_METHOD(standard_deviation)
      ;
  }

}
