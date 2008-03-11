// Pyrticle - Particle in Cell in Python
// Python wrapper for reconstruction bits
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





#include "rec_shape.hpp"
#include "wrap_helpers.hpp"




using namespace pyrticle;
namespace python = boost::python;




void expose_reconstructor()
{
  using python::arg;

  {
    typedef shape_function cl;
    python::class_<cl>("ShapeFunction", 
        python::init<double, unsigned, python::optional<double> >(
          (arg("radius"), arg("dimensions"), arg("alpha"))))
      .add_property("radius", &cl::radius)
      .add_property("exponent", &cl::exponent)
      .def("__call__", &cl::operator())
      ;
  }
}
