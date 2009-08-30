"""Python interface for depositors"""

from __future__ import division

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




from pyrticle.deposition import Depositor
import pyrticle._internal as _internal




class DiracDepositor(Depositor):
    name = "Shape"

    def initialize(self, method):
        Depositor.initialize(self, method)

        backend_class = getattr(_internal, "DiracDepositor" 
                + method.get_dimensionality_suffix())
        self.backend = backend_class(method.mesh_data)

    def make_state(self, state):
        return self.backend.DepositorState()

    def set_shape_function(self, state, sf):
        Depositor.set_shape_function(self, state, sf)
        self.backend.shape_function = sf

    def note_move(self, state, orig, dest, size):
        pass

    def note_change_size(self, state, count):
        pass
