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




import pyrticle._internal as _internal




class Depositor(object):
    def __init__(self):
        self.log_constants = {}
    
    def initialize(self, method):
        self.method = method
        self.shape_function = None

    def set_shape_function(self, sf):
        self.shape_function = sf

    def add_instrumentation(self, mgr):
        mgr.set_constant("depositor", self.name)

        for key, value in self.log_constants.iteritems():
            mgr.set_constant(key, value)

        from pytools.log import IntervalTimer, EventCounter,\
                time_and_count_function

        self.deposit_timer = IntervalTimer(
                "t_deposit",
                "Time spent depositing")
        self.deposit_counter = EventCounter(
                "n_deposit",
                "Number of depositions")

        self.deposit_densities = time_and_count_function(
                self.deposit_densites,
                self.deposit_timer,
                self.deposit_counter,
                1+self.method.dimensions_velocity)

        self.deposit_j = time_and_count_function(
                self.deposit_j,
                self.deposit_timer,
                self.deposit_counter,
                self.method.dimensions_velocity)

        self.deposit_rho = time_and_count_function(
                self.deposit_rho,
                self.deposit_timer,
                self.deposit_counter)

        mgr.add_quantity(self.deposit_timer)
        mgr.add_quantity(self.deposit_counter)

    def clear_particles(self):
        pass

    def deposit_hook(self):
        if self.shape_function is None:
            raise RuntimeError, "shape function never set"

    def make_state(self, state):
        return self.backend.DepositorState()

    def note_move(self, state, orig, dest, size):
        state.depositor_state.note_move(orig, dest, size)

    def note_change_size(self, state, count):
        state.depositor_state.note_change_size(count)

    def _deposit_densities(self, state, velocities, pslice):
        return _internal.deposit_densities(
                self.backend,
                state.depositor_state,
                state.particle_state,
                len(self.method.discretization),
                velocities, pslice)

    def _deposit_j(self, state, velocities, pslice):
        return _internal.deposit_j(
                self.backend,
                state.depositor_state,
                state.particle_state,
                len(self.method.discretization),
                velocities, pslice)

    def _deposit_rho(self, state, pslice):
        return _internal.deposit_rho(
            self.backend,
            state.depositor_state,
            state.particle_state,
            len(self.method.discretization),
            pslice)

    def deposit_densites(self, state, velocities):
        self.deposit_hook()
        rho, j =  self._deposit_densities(state, velocities, slice(None))

        return rho, j

    def deposit_j(self, state, velocities):
        self.deposit_hook()
        j = self._deposit_j(state, velocities, slice(None))

        return j

    def deposit_rho(self, state):
        self.deposit_hook()

        rho = self._deposit_rho(state, slice(None))

        return rho

    def upkeep(self):
        pass

    def rhs(self, state):
        return 0

    def advance_state(self, state, rhs):
        return state





