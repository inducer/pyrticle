"""Python interface for reconstructors"""

__copyright__ = "Copyright (C) 2007, 2008 Andreas Kloeckner"

__license__ = """
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see U{http://www.gnu.org/licenses/}.
"""




class ShapeFunctionReconstructor:
    name = "Shape"

    def initialize(self, cloud):
        cloud.pic_algorithm.set_radius(0.5*cloud.mesh_data.min_vertex_distance())



class AdvectiveReconstructor:
    name = "Advective"

    def initialize(self, cloud):
        pass
