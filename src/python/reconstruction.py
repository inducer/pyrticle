"""Python interface for reconstructors"""

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




import pytools.log
import pylinear.array as num
import pylinear.computation as comp




class ShapeFunctionReconstructor(object):
    name = "Shape"

    def initialize(self, cloud):
        cloud.pic_algorithm.set_radius(0.5*cloud.mesh_data.min_vertex_distance())

    def add_instrumentation(self, mgr):
        pass

    def add_particle_hook(self, pn):
        pass

    def rhs(self):
        return 0

    def add_rhs(self, rhs):
        return 0









class ActiveAdvectiveElements (pytools.log.LogQuantity):
    def __init__(self, reconstructor, name="n_advec"):
        pytools.log.LogQuantity.__init__(self, name, "1", "#active advective_elements")
        self.reconstructor = reconstructor

    def __call__(self):
        return self.reconstructor.cloud.pic_algorithm.active_elements




class AdvectiveReconstructor(object):
    name = "Advective"

    def __init__(self):
        from pyrticle.tools import DOFShiftForwarder
        self.rho_shift_signaller = DOFShiftForwarder()

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

    def add_instrumentation(self, mgr):
        mgr.add_quantity(self.element_activation_counter)
        mgr.add_quantity(self.element_kill_counter)
        mgr.add_quantity(self.advective_rhs_timer)
        mgr.add_quantity(self.active_elements_log)

    def initialize(self, cloud):
        self.cloud = cloud
        discr = cloud.mesh_data.discr
        
        eg, = discr.element_groups
        (fg, fmm), = discr.face_groups
        ldis = eg.local_discretization

        from hedge.mesh import TAG_ALL
        bdry = discr._get_boundary(TAG_ALL)

        (bdry_fg, _), = bdry.face_groups_and_ldis

        cloud.pic_algorithm.setup_advective_reconstructor(
                len(ldis.face_indices()),
                ldis.node_count(),
                ldis.mass_matrix(),
                ldis.inverse_mass_matrix(),
                fmm,
                fg,
                bdry_fg)

        for i, diffmat in enumerate(ldis.differentiation_matrices()):
            cloud.pic_algorithm.add_local_diff_matrix(i, diffmat)

        self.radius = 0.5*cloud.mesh_data.min_vertex_distance()

        cloud.pic_algorithm.rho_dof_shift_listener = self.rho_shift_signaller

    def add_particle_hook(self, pn):
        self.cloud.pic_algorithm.add_advective_particle(self.radius, pn)

    def rhs(self):
        from pyrticle.tools import DOFShiftableVector
        self.advective_rhs_timer.start()
        result =  DOFShiftableVector(
                self.cloud.pic_algorithm.get_advective_particle_rhs(self.cloud.raw_velocities()),
                self.rho_shift_signaller
                )
        self.advective_rhs_timer.stop()
        self.element_activation_counter.transfer(
                self.cloud.pic_algorithm.element_activation_counter)
        self.element_kill_counter.transfer(
                self.cloud.pic_algorithm.element_kill_counter)

        return result

    def add_rhs(self, rhs):
        from pyrticle.tools import DOFShiftableVector
        self.cloud.pic_algorithm.apply_advective_particle_rhs(
                DOFShiftableVector.unwrap(rhs))
