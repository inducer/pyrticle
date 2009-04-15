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
    def __init__(self, depositor, name="n_advec_elements"):
        LogQuantity.__init__(self, name, "1", "#active advective elements")
        self.depositor = depositor

    def __call__(self):
        return self.depositor.cloud.pic_algorithm.active_elements




class AdvectiveDepositor(Depositor):
    name = "Advective"

    def __init__(self, activation_threshold=1e-5, kill_threshold=1e-3, 
            filter_amp=None, filter_order=None, 
            upwind_alpha=1,
            ):
        Depositor.__init__(self)
        _internal.NumberShiftListener.__init__(self)

        from pyrticle.tools import NumberShiftMultiplexer
        self.rho_shift_signaller = NumberShiftMultiplexer()

        self.activation_threshold = activation_threshold
        self.kill_threshold = kill_threshold
        self.upwind_alpha = upwind_alpha

        self.shape_function = None

        resize_state(m_dofs_per_element * 1024);

        self.filter_amp = filter_amp
        self.filter_order = filter_order

        if filter_amp is not None:
            from hedge.discretization import ExponentialFilterResponseFunction
            self.filter_response = ExponentialFilterResponseFunction(
                    filter_amp, filter_order)
        else:
            self.filter_response = None

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
        self.active_elements_log = ActiveAdvectiveElements(self)


    def initialize(self, cloud):
        Depositor.initialize(self, cloud)

        cloud.particle_number_shift_signaller.subscribe(self)

        discr = cloud.discretization
        
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

        cloud.pic_algorithm.setup_advective_depositor(
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
            cloud.pic_algorithm.add_local_diff_matrix(i, diffmat)

        cloud.pic_algorithm.rho_dof_shift_listener = self.rho_shift_signaller

    def add_instrumentation(self, mgr):
        Depositor.add_instrumentation(self, mgr)

        mgr.add_quantity(self.element_activation_counter)
        mgr.add_quantity(self.element_kill_counter)
        mgr.add_quantity(self.advective_rhs_timer)
        mgr.add_quantity(self.active_elements_log)

        mgr.set_constant("el_activation_threshold", self.activation_threshold)
        mgr.set_constant("el_kill_threshold", self.kill_threshold)
        mgr.set_constant("adv_upwind_alpha", self.upwind_alpha)

        mgr.set_constant("filter_amp", self.filter_amp)
        mgr.set_constant("filter_amp", self.filter_order)

    def set_shape_function(self, sf):
        Depositor.set_shape_function(self, sf)

        self.cloud.pic_algorithm.clear_advective_particles()
        for pn in xrange(len(self.cloud)):
            self.cloud.pic_algorithm.add_advective_particle(sf, pn)

    def note_change_size(self, new_size):
        pic = self.cloud.pic_algorithm

        if (self.shape_function is not None 
                and new_size > pic.count_advective_particles()):
            for pn in range(pic.count_advective_particles(), new_size):
                pic.add_advective_particle(self.shape_function, pn)

    def clear_particles(self):
        Depositor.clear_particles(self)
        self.cloud.pic_algorithm.clear_advective_particles()

    def upkeep(self):
        self.cloud.pic_algorithm.perform_depositor_upkeep()

    def rhs(self, state):
        from pyrticle.tools import NumberShiftableVector
        self.advective_rhs_timer.start()
        result =  NumberShiftableVector(
                self.cloud.pic_algorithm.get_advective_particle_rhs(self.cloud.velocities()),
                signaller=self.rho_shift_signaller
                )
        self.advective_rhs_timer.stop()
        self.element_activation_counter.transfer(
                self.cloud.pic_algorithm.element_activation_counter)
        self.element_kill_counter.transfer(
                self.cloud.pic_algorithm.element_kill_counter)

        return result

    def advance_state(self, rhs):
        from pyrticle.tools import NumberShiftableVector
        self.cloud.pic_algorithm.apply_advective_particle_rhs(
                NumberShiftableVector.unwrap(rhs))
