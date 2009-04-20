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




from pytools.log import LogQuantity
from pyrticle.deposition import Depositor
import pyrticle._internal as _internal
import numpy




class ActiveAdvectiveElements(LogQuantity):
    def __init__(self, observer, name="n_advec_elements"):
        LogQuantity.__init__(self, name, "1", "#active advective elements")
        self.observer = observer

    def __call__(self):
        return self.observer.state.depositor_state.active_elements




class AdvectiveDepositor(Depositor):
    name = "Advective"

    def __init__(self, activation_threshold=1e-5, kill_threshold=1e-3, 
            filter_amp=None, filter_order=None, 
            upwind_alpha=1,
            ):
        Depositor.__init__(self)

        self.activation_threshold = activation_threshold
        self.kill_threshold = kill_threshold
        self.upwind_alpha = upwind_alpha

        self.shape_function = None

        self.filter_amp = filter_amp
        self.filter_order = filter_order

        if filter_amp is not None:
            from hedge.discretization import ExponentialFilterResponseFunction
            self.filter_response = ExponentialFilterResponseFunction(
                    filter_amp, filter_order)
        else:
            self.filter_response = None

    def initialize(self, method):
        Depositor.initialize(self, method)

        discr = method.discretization
        
        eg, = discr.element_groups
        fg, = discr.face_groups
        ldis = eg.local_discretization

        from hedge.mesh import TAG_ALL
        bdry = discr.get_boundary(TAG_ALL)

        bdry_fg, = bdry.face_groups

        if self.filter_response:
            from hedge.discretization import Filter
            filter = Filter(discr, self.filter_response)
            filter_mat, = filter.filter_matrices
        else:
            filter_mat = numpy.zeros((0,0))

        backend_class = getattr(_internal, "AdvectiveDepositor" 
                + method.get_dimensionality_suffix())
        self.backend = backend_class(method.mesh_data,
                len(ldis.face_indices()),
                ldis.node_count(),
                ldis.mass_matrix(),
                ldis.inverse_mass_matrix(),
                filter_mat,
                ldis.face_mass_matrix(),
                fg,
                bdry_fg,
                self.activation_threshold,
                self.kill_threshold,
                self.upwind_alpha)

        for i, diffmat in enumerate(ldis.differentiation_matrices()):
            self.backend.add_local_diff_matrix(i, diffmat)

    def make_state(self, state):
        state = self.backend.DepositorState()

        from pyrticle.tools import NumberShiftMultiplexer
        state.rho_dof_shift_listener = NumberShiftMultiplexer()

        eg, = self.method.discretization.element_groups
        ldis = eg.local_discretization
        #state.resize_rho(eg.local_discretization.node_count() * 1024)

        return state

    def add_instrumentation(self, mgr, observer):
        Depositor.add_instrumentation(self, mgr, observer)

        # instrumentation 
        from pytools.log import IntervalTimer, EventCounter
        self.element_activation_counter = EventCounter(
                "n_el_activations",
                "#Advective rec. elements activated this timestep")
        self.element_kill_counter = EventCounter(
                "n_el_kills",
                "#Advective rec. elements retired this timestep")
        self.advective_rhs_timer = IntervalTimer(
                "t_advective_rhs",
                "Time spent evaluating advective RHS")
        self.active_elements_log = ActiveAdvectiveElements(observer)

        mgr.add_quantity(self.element_activation_counter)
        mgr.add_quantity(self.element_kill_counter)
        mgr.add_quantity(self.advective_rhs_timer)
        mgr.add_quantity(self.active_elements_log)

        mgr.set_constant("el_activation_threshold", self.activation_threshold)
        mgr.set_constant("el_kill_threshold", self.kill_threshold)
        mgr.set_constant("adv_upwind_alpha", self.upwind_alpha)

        mgr.set_constant("filter_amp", self.filter_amp)
        mgr.set_constant("filter_amp", self.filter_order)

    def set_shape_function(self, state, sf):
        Depositor.set_shape_function(self, state, sf)

        state.depositor_state.clear()
        for pn in xrange(len(state)):
            self.backend.add_advective_particle(
                    state.depositor_state, 
                    state.particle_state, 
                    sf, pn)

    def note_move(self, state, orig, dest, size):
        self.backend.note_move(state.depositor_state, orig, dest, size)

    def note_change_size(self, state, new_size):
        self.backend.note_change_size(state.depositor_state, new_size)

        if (self.shape_function is not None 
                and new_size > state.depositor_state.count_advective_particles()):
            for pn in range(
                    state.depositor_state.count_advective_particles(), 
                    new_size):
                self.backend.add_advective_particle(
                        state.depositor_state,
                        state.particle_state,
                        self.shape_function, pn)

    def clear_particles(self):
        Depositor.clear_particles(self)
        self.cloud.pic_algorithm.clear_advective_particles()

    def upkeep(self, state):
        self.backend.perform_depositor_upkeep(
                state.depositor_state,
                state.particle_state)

    def rhs(self, state):
        from pyrticle.tools import NumberShiftableVector
        self.advective_rhs_timer.start()
        result =  NumberShiftableVector(
                self.backend.get_advective_particle_rhs(
                    state.depositor_state,
                    state.particle_state,
                    self.method.velocities(state)),
                signaller=state.depositor_state.rho_dof_shift_listener
                )
        self.advective_rhs_timer.stop()
        self.element_activation_counter.transfer(
                state.depositor_state.element_activation_counter)
        self.element_kill_counter.transfer(
                state.depositor_state.element_kill_counter)

        return result

    def advance_state(self, state, rhs):
        from pyrticle.tools import NumberShiftMultiplexer, NumberShiftableVector
        return self.backend.apply_advective_particle_rhs(
                state.depositor_state,
                state.particle_state,
                NumberShiftableVector.unwrap(rhs),
                NumberShiftMultiplexer())
