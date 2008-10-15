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




void expose_tools();
void expose_grid();
void expose_meshdata();
void expose_pic();
void expose_pusher();
void expose_reconstructor();
/*
void expose_shape_pic();
void expose_normshape_pic();
void expose_advective_pic();
void expose_grid_pic();
void expose_grid_find_pic();
*/




BOOST_PYTHON_MODULE(_internal)
{
  expose_tools();
  expose_grid();
  expose_meshdata();
  expose_pic();
  expose_pusher();
  expose_reconstructor();
  /*
  expose_shape_pic();
  expose_normshape_pic();
  expose_advective_pic();
  expose_grid_pic();
  expose_grid_find_pic();
  */
}
