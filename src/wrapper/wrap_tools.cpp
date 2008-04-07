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
        const py_vector &vec,
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
  template <typename ValueType>
  python::handle<> lu_wrapper(const pyublas::numpy_matrix<ValueType> &a)
  {
    namespace lapack = boost::numeric::bindings::lapack;

    typedef pyublas::numpy_matrix<ValueType> matrix_t;
    typedef boost::numeric::ublas::matrix<
      ValueType, boost::numeric::ublas::column_major> col_matrix_t;

    const unsigned piv_len = std::min(a.size1(), a.size2());

    col_matrix_t temp(a);
    boost::numeric::ublas::vector<int> piv(piv_len);

    int info = lapack::getrf(temp, piv);
    if (info < 0)
      throw std::runtime_error("invalid argument to getrf");
    
    matrix_t l(a.size1(), a.size2());
    l.clear();
    matrix_t u(a.size1(), a.size2());
    u.clear();

    for (typename matrix_t::size_type i = 0; i < a.size1(); i++)
    {
      unsigned j = 0;
      for (; j < std::min(i, a.size2()); j++) l(i,j) = temp(i,j);
      l(i,i) = 1;
      for (; j < a.size2(); j++) u(i,j) = temp(i,j);
    }

    boost::numeric::ublas::vector<int> permut(piv_len);
    for (unsigned i = 0; i < piv_len; i++) 
      permut[i] = i;
    for (unsigned i = 0; i < piv_len; i++) 
      std::swap(permut[i], permut[piv[i]-1]);

    python::list py_permut;
    for (unsigned i = 0; i < piv_len; i++)
      py_permut.append(permut[i]);
    
    return python::handle<>(Py_BuildValue("(NNO)", 
        l.to_python().release(), 
        u.to_python().release(),
        py_permut.ptr()));
  }
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
      .DEF_SIMPLE_METHOD(variance)
      .DEF_SIMPLE_METHOD(standard_deviation)
      ;
  }

  python::def("lu", lu_wrapper<double>);
}
