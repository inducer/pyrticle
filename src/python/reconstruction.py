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




import pyrticle._internal as _internal
import pytools.log
import numpy
import numpy.linalg as la




class Reconstructor(object):
    def __init__(self):
        pass
    
    def initialize(self, cloud):
        self.cloud = cloud
        self.shape_function = None

    def set_shape_function(self, sf):
        self.shape_function = sf

    def add_instrumentation(self, mgr):
        pass

    def clear_particles(self):
        pass

    def reconstruct_hook(self):
        if self.shape_function is None:
            raise RuntimeError, "shape function never set"

    def rhs(self):
        return 0

    def add_rhs(self, rhs):
        return 0




class ShapeFunctionReconstructor(Reconstructor):
    name = "Shape"

    def set_shape_function(self, sf):
        Reconstructor.set_shape_function(self, sf)
        self.cloud.pic_algorithm.shape_function = sf





class NormalizedShapeFunctionReconstructor(Reconstructor):
    name = "NormShape"

    def initialize(self, cloud):
        Reconstructor.initialize(self, cloud)

        eg, = cloud.mesh_data.discr.element_groups
        ldis = eg.local_discretization

        cloud.pic_algorithm.setup_normalized_shape_reconstructor(
                ldis.mass_matrix())

    def add_instrumentation(self, mgr):
        Reconstructor.add_instrumentation(self, mgr)

        from pyrticle.log import StatsGathererLogQuantity
        mgr.add_quantity(StatsGathererLogQuantity(
            self.cloud.pic_algorithm.normalization_stats,
            "normshape_norm", "1", 
            "normalization constants applied during reconstruction"))

        mgr.add_quantity(StatsGathererLogQuantity(
            self.cloud.pic_algorithm.centroid_distance_stats,
            "normshape_centroid_dist", "m", 
            "distance of shape center from element centroid"))

        mgr.add_quantity(StatsGathererLogQuantity(
            self.cloud.pic_algorithm.el_per_particle_stats,
            "normshape_el_per_particle", "1", 
            "number of elements per particle"))

    def set_shape_function(self, sf):
        Reconstructor.set_shape_function(self, sf)
        self.cloud.pic_algorithm.shape_function = sf





# advective reconstruction ----------------------------------------------------
class ActiveAdvectiveElements(pytools.log.LogQuantity):
    def __init__(self, reconstructor, name="n_advec_elements"):
        pytools.log.LogQuantity.__init__(self, name, "1", "#active advective elements")
        self.reconstructor = reconstructor

    def __call__(self):
        return self.reconstructor.cloud.pic_algorithm.active_elements




class AdvectiveReconstructor(Reconstructor, _internal.NumberShiftListener):
    name = "Advective"

    def __init__(self, activation_threshold=1e-5, kill_threshold=1e-3, 
            filter_amp=None, filter_order=None, 
            upwind_alpha=1):
        Reconstructor.__init__(self)
        _internal.NumberShiftListener.__init__(self)

        from pyrticle.tools import NumberShiftForwarder
        self.rho_shift_signaller = NumberShiftForwarder()

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
        Reconstructor.initialize(self, cloud)

        cloud.particle_number_shift_signaller.subscribe(self)

        discr = cloud.mesh_data.discr
        
        eg, = discr.element_groups
        (fg, fmm), = discr.face_groups
        ldis = eg.local_discretization

        from hedge.mesh import TAG_ALL
        bdry = discr._get_boundary(TAG_ALL)

        (bdry_fg, _), = bdry.face_groups_and_ldis

        if self.filter_response:
            from hedge.discretization import Filter
            filter = Filter(discr, self.filter_response)
            filter_mat, = filter.filter_matrices
        else:
            filter_mat = numpy.zeros((0,0))

        cloud.pic_algorithm.setup_advective_reconstructor(
                len(ldis.face_indices()),
                ldis.node_count(),
                ldis.mass_matrix(),
                ldis.inverse_mass_matrix(),
                filter_mat,
                fmm,
                fg,
                bdry_fg,
                self.activation_threshold,
                self.kill_threshold,
                self.upwind_alpha)

        for i, diffmat in enumerate(ldis.differentiation_matrices()):
            cloud.pic_algorithm.add_local_diff_matrix(i, diffmat)

        cloud.pic_algorithm.rho_dof_shift_listener = self.rho_shift_signaller

    def add_instrumentation(self, mgr):
        Reconstructor.add_instrumentation(self, mgr)

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
        Reconstructor.set_shape_function(self, sf)

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
        Reconstructor.clear_particles(self)
        self.cloud.pic_algorithm.clear_advective_particles()

    def rhs(self):
        from pyrticle.tools import NumberShiftableVector
        self.advective_rhs_timer.start()
        result =  NumberShiftableVector(
                self.cloud.pic_algorithm.get_advective_particle_rhs(self.cloud.velocities()),
                multiplier=1,
                signaller=self.rho_shift_signaller
                )
        self.advective_rhs_timer.stop()
        self.element_activation_counter.transfer(
                self.cloud.pic_algorithm.element_activation_counter)
        self.element_kill_counter.transfer(
                self.cloud.pic_algorithm.element_kill_counter)

        return result

    def add_rhs(self, rhs):
        from pyrticle.tools import NumberShiftableVector
        self.cloud.pic_algorithm.apply_advective_particle_rhs(
                NumberShiftableVector.unwrap(rhs))




# grid reconstruction ---------------------------------------------------------
class SingleBrick:
    def __init__(self, overresolve=1.5):
        self.overresolve = overresolve

    def __call__(self, discr):
        from hedge.discretization import integral, ones_on_volume
        mesh_volume = integral(discr, ones_on_volume(discr))
        dx =  (mesh_volume/len(discr))**(1/discr.dimensions) \
                / self.overresolve

        mesh = discr.mesh
        bbox_min, bbox_max = mesh.bounding_box
        bbox_size = bbox_max-bbox_min
        dims = numpy.asarray(bbox_size/dx, dtype=numpy.uint32)
        stepwidths = bbox_size/(dims-1)
        yield stepwidths, bbox_min, dims




class GridReconstructor(Reconstructor):
    name = "Grid"

    def __init__(self, brick_generator=SingleBrick(), el_tolerance=0):
        self.brick_generator = brick_generator
        self.el_tolerance = el_tolerance

    def initialize(self, cloud):
        Reconstructor.initialize(self, cloud)

        discr = cloud.mesh_data.discr

        from pyublas import why_not
        for stepwidths, origin, dims in self.brick_generator(discr):
            self.cloud.pic_algorithm.add_brick(
                    why_not(stepwidths), why_not(origin), why_not(dims, dtype=numpy.uint32))

        eg, = discr.element_groups
        ldis = eg.local_discretization

        self.cloud.pic_algorithm.commit_bricks(
                ldis.vandermonde(),
                ldis.basis_functions(),
                self.el_tolerance)

    def set_shape_function(self, sf):
        Reconstructor.set_shape_function(self, sf)
        self.cloud.pic_algorithm.shape_function = sf

    def write_grid_rho(self, silo):
        dims = self.cloud.dimensions_mesh
        pic = self.cloud.pic_algorithm

        for i_brick, brick in enumerate(pic.bricks):
            coords = [
                numpy.arange(
                    brick.origin[axis], 
                    brick.origin[axis] + (brick.dimensions[axis]-0.5) * brick.stepwidths[axis],
                    brick.stepwidths[axis])
                for axis in xrange(dims)]
            for axis in xrange(dims):
                assert len(coords[axis]) == brick.dimensions[axis]

            mname = "structmesh%d" % i_brick
            vname = "rho_struct%d" % i_brick
            silo.put_quadmesh(mname, coords)

            from pylo import DB_NODECENT
            silo.put_quadvar1(vname, mname, pic.get_rec_debug_quantity("rho_grid"),
                    [int(x) for x in brick.dimensions], # get rid of array scalars
                    DB_NODECENT)
