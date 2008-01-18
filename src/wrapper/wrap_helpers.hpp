// Pyrticle - Particle in Cell in Python
// Little bits of helpfulness for wrapping things using Boost Python
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





#ifndef _BADFJAH_PYRTICLE_WRAP_HELPERS_HPP_INCLUDED
#define _BADFJAH_PYRTICLE_WRAP_HELPERS_HPP_INCLUDED




#include <string>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>




#define PYTHON_ERROR(TYPE, REASON) \
{ \
  PyErr_SetString(PyExc_##TYPE, REASON); \
  throw boost::python::error_already_set(); \
}

#define COPY_PY_LIST(CPP_LIST, PY_LIST, EXTRACT_TYPE) \
        for (unsigned i = 0; i < unsigned(len(PY_LIST)); i++) \
          (CPP_LIST).push_back( \
              python::extract<EXTRACT_TYPE>((PY_LIST)[i]));

#define DEF_SIMPLE_METHOD(NAME) \
  def(#NAME, &cl::NAME)
#define DEF_PURE_VIRTUAL_METHOD(NAME) \
  def(#NAME, boost::python::pure_virtual(&cl::NAME))
#define DEF_SIMPLE_FUNCTION(NAME) \
  def(#NAME, &NAME)
#define DEF_RO_MEMBER(NAME) \
  def_readonly(#NAME, &cl::m_##NAME)
#define DEF_RW_MEMBER(NAME) \
  def_readwrite(#NAME, &cl::m_##NAME)




namespace pyrticle {
  template <class T>
  inline PyObject *manage_new_object(T *obj)
  {
    typename boost::python::manage_new_object::apply<T *>::type 
      result_converter;
    return result_converter(obj);
  }




  template <class T>
  class no_compare_indexing_suite :
    public boost::python::vector_indexing_suite<T, false, no_compare_indexing_suite<T> >
  {
    public:
      static bool contains(T &container, typename T::value_type const &key)
      { PYTHON_ERROR(NotImplementedError, "containment checking not supported on this container"); }
  };




  template <class T>
  std::auto_ptr<std::vector<T> > construct_vector(boost::python::object iterable)
  {
    std::auto_ptr<std::vector<T> > result(new std::vector<T>());
    copy(
        boost::python::stl_input_iterator<T>(iterable),
        boost::python::stl_input_iterator<T>(),
        std::back_inserter(*result));
    return result;
  }




  template <class ValueType>
  void expose_std_vector(const char *name)
  {
    typedef std::vector<ValueType> cl;
    boost::python::class_<cl>((std::string(name)+"Vector").c_str())
      .def("__init__", boost::python::make_constructor(
            construct_vector<ValueType>))
      .def(no_compare_indexing_suite<cl>())
      .DEF_SIMPLE_METHOD(reserve)
      ;
  }
}





#endif
