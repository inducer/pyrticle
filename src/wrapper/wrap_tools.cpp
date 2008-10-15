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




#include <pyublas/numpy.hpp>
#include <pyublas/python_helpers.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/numeric/bindings/lapack/gesv.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include "wrap_tuples.hpp"
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
    void store_mesh_vis_vector(
        const char *name, const py_vector &vec) const
    {
      this->get_override("store_mesh_vis_vector")(name, vec);
    }
    void store_particle_vis_vector(
        const char *name, const py_vector &vec) const
    {
      this->get_override("store_particle_vis_vector")(name, vec);
    }
  };



  /*

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




  struct boundary_hit_listener_wrap : 
    boundary_hit_listener,
    python::wrapper<boundary_hit_listener>
  {
      void note_boundary_hit(particle_number pn) const 
      {
        this->get_override("note_boundary_hit")(pn);
      }
  };

  */



  struct warning_listener_wrap : 
    warning_listener,
    python::wrapper<warning_listener>
  {
      void note_warning(
          std::string const &message,
          std::string const &filename,
          unsigned lineno
          ) const 
      {
        this->get_override("note_warning")(message, filename, lineno);
      }
  };




  // lu -----------------------------------------------------------------------
  python::handle<> lu_wrapper(const py_matrix &a)
  {
    namespace lapack = boost::numeric::bindings::lapack;

    const unsigned piv_len = std::min(a.size1(), a.size2());

    py_int_vector piv(piv_len);
    py_fortran_matrix a_copy(1*a);

    int info = lapack::getrf(a_copy.as_ublas(), piv);
    if (info < 0)
      throw std::runtime_error("invalid argument to getrf");
    if (info > 0)
      throw std::runtime_error("singular matrix in getrf");
    
    return python::handle<>(Py_BuildValue("(NN)", 
        a_copy.to_python().release(), 
        piv.to_python().release()
        ));
  }



  template <class VecT>
  void expose_box(std::string const &tp_name)
  {
    typedef box<VecT> cl;
    python::class_<cl>(("Box"+tp_name).c_str(),
        python::init<VecT const &, VecT const &>())
      .DEF_BYVAL_RW_MEMBER(lower)
      .DEF_BYVAL_RW_MEMBER(upper)
      .def("contains", 
          (bool (cl::*)(py_vector const &) const) &cl::contains,
          python::args("point"))
      .def("contains", 
          (bool (cl::*)(py_vector const &, double) const) &cl::contains,
          python::args("point", "threshold"))
      .DEF_SIMPLE_METHOD(is_empty)
      .def("intersect", &cl::template intersect<VecT>)

      .def(python::self == python::self)
      .def(python::self != python::self)
      ;
  }
}




void expose_tools()
{
  python::def("asinh", (double (*)(double)) boost::math::asinh);
  python::def("acosh", (double (*)(double)) boost::math::acosh);
  python::def("gamma", (double (*)(double)) boost::math::tgamma);

  expose_box<bounded_vector>("Float");
  expose_box<bounded_int_vector>("Int");
  {
    typedef event_counter cl;
    python::class_<cl>("EventCounter")
      .DEF_SIMPLE_METHOD(get)
      .DEF_SIMPLE_METHOD(pop)
      .DEF_SIMPLE_METHOD(tick)
      ;
  }

  {
    typedef std::vector<unsigned int> cl;
    python::class_<cl>("UnsignedVector")
      .def(python::vector_indexing_suite<cl>())
      .DEF_SIMPLE_METHOD(clear)
      .DEF_SIMPLE_METHOD(reserve)
      ;
  }
  
  {
    typedef std::vector<unsigned long> cl;
    python::class_<cl>("Unsigned32Vector")
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
      .DEF_PURE_VIRTUAL_METHOD(store_mesh_vis_vector)
      .DEF_PURE_VIRTUAL_METHOD(store_particle_vis_vector)
      ;
  }

  /*
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
    typedef boundary_hit_listener cl;
    python::class_<boundary_hit_listener_wrap, 
      boost::shared_ptr<boundary_hit_listener_wrap>,
      boost::noncopyable>
      ("BoundaryHitListener")
      .DEF_PURE_VIRTUAL_METHOD(note_boundary_hit)
      ;
  }
  */

  {
    typedef warning_listener cl;
    python::class_<warning_listener_wrap, 
      boost::shared_ptr<warning_listener_wrap>,
      boost::noncopyable>
      ("WarningListener")
      .DEF_PURE_VIRTUAL_METHOD(note_warning)
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
      .DEF_SIMPLE_METHOD(variance_sample)
      .DEF_SIMPLE_METHOD(standard_deviation_sample)
      .DEF_SIMPLE_METHOD(variance)
      .DEF_SIMPLE_METHOD(standard_deviation)
      ;
  }

  python::def("lu", lu_wrapper);

  {
    typedef polynomial_shape_function cl;
    python::class_<cl>("PolynomialShapeFunction", 
        python::init<double, unsigned, python::optional<double> >(
          (python::args("radius", "dimensions", "alpha"))))
      .add_property("normalizer", &cl::normalizer)
      .add_property("radius", &cl::radius)
      .add_property("exponent", &cl::exponent)
      .def("__call__", 
          (const double (cl::*)(const py_vector &) const)
          &cl::operator())
      .DEF_SIMPLE_METHOD(name)
      ;
  }

  {
    typedef c_infinity_shape_function cl;
    python::class_<cl>("CInfinityShapeFunction", 
        python::init<double, unsigned, double>(
          (python::args("radius", "dimensions", "integral_for_rad1"))))
      .add_property("radius", &cl::radius)
      .def("__call__", 
          (const double (cl::*)(const py_vector &) const)
          &cl::operator())
      .DEF_SIMPLE_METHOD(name)
      ;
  }

  python::def("get_shape_function_name", &shape_function::name);

  python::register_tuple<boost::tuple<py_vector, py_vector> >();
}
