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
import pyublas




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

        from pyrticle.tools import NumberShiftMultiplexer
        self.rho_shift_signaller = NumberShiftMultiplexer()

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
class SingleBrickGenerator:
    def __init__(self, overresolve=1.5, mesh_margin=0):
        self.overresolve = overresolve
        self.mesh_margin = mesh_margin

    def __call__(self, discr):
        from hedge.discretization import integral, ones_on_volume
        mesh_volume = integral(discr, ones_on_volume(discr))
        dx =  (mesh_volume / len(discr)/ self.overresolve)**(1/discr.dimensions) \

        mesh = discr.mesh
        bbox_min, bbox_max = mesh.bounding_box()

        bbox_min -= self.mesh_margin
        bbox_max += self.mesh_margin

        bbox_size = bbox_max-bbox_min
        dims = numpy.asarray(bbox_size/dx, dtype=numpy.int32)
        stepwidths = bbox_size/dims
        yield stepwidths, bbox_min, dims




class GridReconstructor(Reconstructor):
    name = "Grid"

    def __init__(self, brick_generator=SingleBrickGenerator(), 
            el_tolerance=0.2,
            max_extra_points=20,
            enforce_continuity=False,
            method="simplex_enlarge"):
        self.brick_generator = brick_generator
        self.el_tolerance = el_tolerance
        self.max_extra_points = max_extra_points
        self.enforce_continuity = enforce_continuity
        self.method = method

    def initialize(self, cloud):
        Reconstructor.initialize(self, cloud)

        discr = cloud.mesh_data.discr

        pic = self.cloud.pic_algorithm

        if self.enforce_continuity:
            self.prepare_average_groups()
        else:
            pic.average_group_starts.append(0)

        from pyrticle._internal import Brick
        for i, (stepwidths, origin, dims) in enumerate(
                self.brick_generator(discr)):
            pic.bricks.append(Brick(i, pic.grid_node_count(),
                        stepwidths, origin, dims))

        if self.method == "simplex_extra":
            self.prepare_with_pointwise_projection_and_extra_points()
        elif self.method == "simplex_enlarge":
            self.prepare_with_pointwise_projection_and_enlargement()
        elif self.method == "brick":
            self.prepare_with_brick_interpolation()
        else:
            raise RuntimeError, "invalid rec_grid submethod specified"

    def find_containing_brick(self, pt):
        for brk in self.cloud.pic_algorithm.bricks:
            if brk.bounding_box().contains(pt):
                return brk
        raise RuntimeError, "no containing brick found for point"

    def prepare_average_groups(self):
        discr = self.cloud.mesh_data.discr
        pic = self.cloud.pic_algorithm

        avg_group_finder = {}
        avg_groups = []

        for fg, fmm in discr.face_groups:
            for fp in fg.face_pairs:
                for el_idx, opp_el_idx in zip(
                        fg.index_lists[fp.face_index_list_number],
                        fg.index_lists[fp.opp_face_index_list_number]):
                    idx1 = fp.el_base_index + el_idx
                    idx2 = fp.opp_el_base_index + opp_el_idx

                    ag1 = avg_group_finder.get(idx1)
                    ag2 = avg_group_finder.get(idx2)

                    if ag1 is None and ag2 is None:
                        ag = set([idx1, idx2])
                        avg_groups.append(ag)
                    elif (ag1 is not None and ag2 is not None and
                            ag1 is not ag2):
                        # need to merge
                        ag1.update(ag2)
                        ag2.clear()
                        ag = ag1
                    else:
                        ag = ag1 or ag2
                        ag.add(idx1)
                        ag.add(idx2)

                    for idx in ag:
                        avg_group_finder[idx] = ag

        for ag in avg_groups:
            pic.average_groups.extend(ag)
            pic.average_group_starts.append(len(pic.average_groups))

        print len(avg_groups), "average groups"

    def find_points_in_element(self, el, el_tolerance):
        from pyrticle._internal import BrickIterator, ElementOnGrid, BoxFloat
        eog = ElementOnGrid()
        eog.element_number = el.id

        discr = self.cloud.mesh_data.discr
        pic = self.cloud.pic_algorithm

        el_bbox = BoxFloat(*el.bounding_box(discr.mesh.points))

        grid_node_count = pic.grid_node_count()

        # enlarge the element bounding box by the mapped tolerance
        scaled_tolerance = el_tolerance * la.norm(el.map.matrix, 2)
        el_bbox.lower -= scaled_tolerance
        el_bbox.upper += scaled_tolerance

        points = []

        # For each element, find all structured points inside the element.
        for brk in pic.bricks:
            brk_and_el = brk.bounding_box().intersect(el_bbox)

            if brk_and_el.is_empty():
                continue

            for coord in BrickIterator(brk, brk.index_range(brk_and_el)):
                point = brk.point(coord)

                in_el = True
                md_elinfo = pic.mesh_data.element_info[el.id]
                for f in md_elinfo.faces:
                    if (numpy.dot(f.normal, point) 
                            - f.face_plane_eqn_rhs > scaled_tolerance):
                        in_el = False
                        break

                if in_el:
                    points.append(point)
                    grid_node_index = brk.index(coord)
                    assert grid_node_index < grid_node_count
                    eog.grid_nodes.append(grid_node_index)

        return eog, points

    def make_pointwise_interpolation_matrix(self, el, ldis, svd, structured_vdm):
        u, s, vt = svd

        point_count = u.shape[0]
        node_count = vt.shape[1]

        thresh = (numpy.finfo(float).eps * max(s.shape) * s[0])

        nonzero_flags = numpy.abs(s) >= thresh
        inv_s = numpy.zeros((len(s),), dtype=float)
        inv_s[nonzero_flags] = 1/s[nonzero_flags]

        # compute the pseudoinverse of the structured
        # Vandermonde matrix
        inv_s_diag = numpy.zeros(
                (node_count, point_count), 
                dtype=float)
        inv_s_diag[:len(s),:len(s)] = numpy.diag(1/s)

        svdm_pinv = numpy.dot(numpy.dot(vt.T, inv_s_diag), u.T)

        # check that it's reasonable
        pinv_resid = la.norm(
            numpy.dot(svdm_pinv, structured_vdm)
            - numpy.eye(node_count))

        if pinv_resid > 1e-8:
            from warnings import warn
            warn("rec_grid: bad pseudoinv precision, element=%d, "
                    "#nodes=%d, #sgridpts=%d, resid=%.5g centroid=%s"
                % (el.id, node_count, point_count, pinv_resid,
                        el.centroid(self.cloud.mesh_data.discr.mesh.points)))
        return numpy.asarray(
                numpy.dot(ldis.vandermonde(), svdm_pinv),
                order="F")

    def prepare_with_pointwise_projection_and_extra_points(self):
        discr = self.cloud.mesh_data.discr
        pic = self.cloud.pic_algorithm

        pic.elements_on_grid.reserve(
                sum(len(eg.members) for eg in discr.element_groups))

        # map brick numbers to [ (point, el_id, el_structured_point_index),...]
        # This is used to write out the C++ extra_points structure further
        # down.
        ep_brick_map = {}

        total_points = 0

        # Iterate over all elements
        for eg in discr.element_groups:
            ldis = eg.local_discretization

            for el in eg.members:
                eog, points = self.find_points_in_element(el, self.el_tolerance)

                # If the structured Vandermonde matrix is singular,
                # add "extra points" to prevent that.
                ep_count = 0
                while True:
                    from hedge.polynomial import generic_vandermonde
                    structured_vdm = generic_vandermonde(
                            [el.inverse_map(x) for x in points], 
                            ldis.basis_functions())

                    u, s, vt = svd = la.svd(structured_vdm)

                    if len(points) >= ldis.node_count():
                        # case 1: theoretically enough points found
                        thresh = (numpy.finfo(float).eps
                                * max(structured_vdm.shape) * s[0])
                        zero_indices = [i for i, si in enumerate(s)
                            if abs(si) < thresh]
                    else:
                        zero_indices = range(vt.shape[0]-len(s), vt.shape[0])
                        assert zero_indices

                    if not zero_indices:
                        break

                    ep_count += 1
                    if ep_count > self.max_extra_points:
                        from warnings import warn
                        warn("rec_grid: could not regularize structured "
                                "vandermonde matrix for el#%d with #ep bound" % el.id)
                        break

                    # Getting here means that a mode
                    # maps to zero on the structured grid.
                    # Find it.

                    # mode linear combination available for zeroed 
                    # mode: use it
                    zeroed_mode = vt[zero_indices[0]]
                    zeroed_mode_nodal = numpy.dot(ldis.vandermonde(), 
                            zeroed_mode)

                    # Then, find the point in that mode with
                    # the highest absolute value.
                    from pytools import argmax
                    max_node_idx = argmax(abs(xi) for xi in zeroed_mode_nodal)
                    start, stop = discr.find_el_range(el.id)

                    new_point = discr.nodes[start+max_node_idx]

                    ep_brick_map.setdefault(
                            self.find_containing_brick(new_point).number,
                            []).append(
                                    (new_point, el.id, len(points)))

                    points.append(new_point)

                    # the final grid_node_number at which this point
                    # will end up is as yet unknown. insert a
                    # placeholder
                    eog.grid_nodes.append(0)

                if ep_count:
                    print "element %d #nodes=%d sgridpt=%d, extra=%d" % (
                            el.id, ldis.node_count(), len(points), ep_count)
                total_points += len(points)

                eog.interpolation_matrix = self.make_pointwise_interpolation_matrix(
                        el, ldis, svd, structured_vdm)

                pic.elements_on_grid.append(eog)

        # fill in the extra points
        ep_brick_starts = [0]
        extra_points = []
        
        grid_node_count = pic.grid_node_count()
        for brk in pic.bricks:
            for pt, el_id, struc_idx in ep_brick_map.get(brk.number, []):
                # replace zero placeholder from above
                pic.elements_on_grid[el_id].grid_nodes[struc_idx] = \
                        grid_node_count + len(extra_points)
                extra_points.append(pt)

            ep_brick_starts.append(len(extra_points))

        pic.first_extra_point = grid_node_count
        pic.extra_point_brick_starts.extend(ep_brick_starts)
        pic.extra_points = numpy.array(extra_points)

        # print some statistics
        print("rec_grid.simplex_extra stats: #nodes: %d, "
                "#points: %d, #points added: %d" % (
                len(discr), total_points, len(extra_points))
        )

    def prepare_with_pointwise_projection_and_enlargement(self):
        tolerance_bound = 1.5

        discr = self.cloud.mesh_data.discr
        pic = self.cloud.pic_algorithm

        pic.elements_on_grid.reserve(
                sum(len(eg.members) for eg in discr.element_groups))

        total_points = 0
        points_added = 0

        # Iterate over all elements
        for eg in discr.element_groups:
            ldis = eg.local_discretization
            basis = ldis.basis_functions()

            for el in eg.members:
                # If the structured Vandermonde matrix is singular,
                # enlarge the element tolerance

                my_tolerance = self.el_tolerance

                orig_point_count = None

                while True:
                    eog, points = self.find_points_in_element(el, my_tolerance)
                    if orig_point_count is None:
                        orig_point_count = len(points)

                    from hedge.polynomial import generic_vandermonde
                    structured_vdm = generic_vandermonde(
                            [el.inverse_map(x) for x in points], 
                            basis)

                    bad_vdm = len(points) < len(basis)
                    if not bad_vdm:
                        try:
                            u, s, vt = svd = la.svd(structured_vdm)
                            thresh = (numpy.finfo(float).eps
                                    * max(structured_vdm.shape) * s[0])
                            zero_indices = [i for i, si in enumerate(s)
                                if abs(si) < thresh]
                            bad_vdm = bool(zero_indices)
                        except la.LinAlgError:
                            bad_vdm = True

                    if not bad_vdm:
                        break

                    my_tolerance += 0.03
                    if my_tolerance >= tolerance_bound:
                        from warnings import warn
                        warn("rec_grid: could not regularize structured "
                                "vandermonde matrix for el#%d by enlargement" % el.id)
                        break

                if my_tolerance > self.el_tolerance:
                    print "element %d #nodes=%d sgridpt=%d, extra tol=%g, #extra points=%d" % (
                            el.id, ldis.node_count(), len(points), 
                            my_tolerance-self.el_tolerance,
                            len(points)-orig_point_count)
                    points_added += len(points)-orig_point_count
                total_points += len(points)

                eog.interpolation_matrix = self.make_pointwise_interpolation_matrix(
                        el, ldis, svd, structured_vdm)

                pic.elements_on_grid.append(eog)

        # we don't need no stinkin' extra points
        pic.extra_point_brick_starts.extend([0]*(len(pic.bricks)+1))

        # print some statistics
        from pytools import average
        print("rec_grid.simplex_enlarge stats: #nodes: %d, "
                "#points: %d, #points added: %d" % (
                len(discr), total_points, points_added)
        )

    def prepare_with_brick_interpolation(self):
        class TensorProductLegendreBasisFunc:
            def __init__(self, n):
                from hedge.polynomial import LegendreFunction
                self.funcs = [LegendreFunction(n_i) for n_i in n]

            def __call__(self, x):
                result = 1
                for f_i, x_i in zip(self.funcs, x):
                    result *= f_i(x_i)
                return result

        def make_legendre_basis(dimensions):
            from pytools import generate_nonnegative_integer_tuples_below
            return [TensorProductLegendreBasisFunc(n)
                for n in generate_nonnegative_integer_tuples_below(dimensions)]

        discr = self.cloud.mesh_data.discr
        pic = self.cloud.pic_algorithm
        bricks = pic.bricks

        pic.elements_on_grid.reserve(
                sum(len(eg.members) for eg in discr.element_groups))

        from pyrticle._internal import BrickIterator, ElementOnGrid, BoxFloat

        total_points = 0

        # Iterate over all elements
        for eg in discr.element_groups:
            ldis = eg.local_discretization

            for el in eg.members:
                el_bbox = BoxFloat(*el.bounding_box(discr.mesh.points))

                scaled_tolerance = self.el_tolerance * la.norm(el.map.matrix, 2)
                el_bbox.lower -= scaled_tolerance
                el_bbox.upper += scaled_tolerance

                # For each brick, find all element nodes that lie in it
                for brk in bricks:
                    eog = ElementOnGrid()
                    eog.element_number = el.id

                    brk_bbox = brk.bounding_box()
                    brk_and_el = brk_bbox.intersect(el_bbox)

                    if brk_and_el.is_empty():
                        continue

                    el_nodes = []
                    el_node_indices = []
                    el_start, el_end = discr.find_el_range(el.id)
                    el_length = el_end-el_start
                    if brk_and_el == el_bbox:
                        # whole element in brick? fantastic.
                        el_nodes = discr.nodes[el_start:el_end]
                        el_node_indices = range(el_length)
                    else:
                        # no? go through the nodes one by one.
                        for i, node in enumerate(discr.nodes[el_start:el_end]):
                            # this containment check has to be exact,
                            # we positively cannot have nodes belong to
                            # two bricks
                            if brk_bbox.contains(node, threshold=0):
                                el_nodes.append(node)
                                el_node_indices.append(i)

                    idx_range = brk.index_range(brk_and_el)
                    lb = make_legendre_basis(idx_range.upper-idx_range.lower)

                    from hedge.polynomial import generic_vandermonde
                    brk_and_el_points = [brk.point(c) 
                            for c in BrickIterator(brk, idx_range)]
                    svdm = generic_vandermonde(
                            points=brk_and_el_points,
                            functions=lb)
                    total_points += len(brk_and_el_points)

                    mixed_vdm_pre = generic_vandermonde(
                            points=el_nodes, functions=lb)

                    if len(el_nodes) < el_length:
                        mixed_vdm = numpy.zeros((el_length, len(lb)),
                                dtype=float)
                        for i, vdm_row in enumerate(mixed_vdm_pre):
                            mixed_vdm[el_node_indices[i]] = vdm_row
                    else:
                        mixed_vdm = mixed_vdm_pre

                    from hedge.tools import leftsolve
                    eog.interpolation_matrix = numpy.asarray(
                            leftsolve(svdm, mixed_vdm),
                            order="Fortran")
                    #print eog.interpolation_matrix.shape
                    #raw_input()

                    eog.grid_nodes.extend(brk.index(c) 
                            for c in BrickIterator(brk, idx_range))
                    pic.elements_on_grid.append(eog)

        # we don't need no stinkin' extra points
        pic.extra_point_brick_starts.extend([0]*(len(pic.bricks)+1))

        print("rec_grid.brick stats: #nodes: %d, #points: %d" % (
                len(discr), total_points))

    def set_shape_function(self, sf):
        Reconstructor.set_shape_function(self, sf)
        self.cloud.pic_algorithm.shape_function = sf

    def write_grid_quantities(self, silo, quantities):
        dims = self.cloud.dimensions_mesh
        vdims = self.cloud.dimensions_velocity
        pic = self.cloud.pic_algorithm

        extra_points = pic.extra_points
        if len(extra_points):
            extra_points = numpy.reshape(extra_points,
                    (len(extra_points)//dims,dims))

            silo.put_pointmesh("rec_grid_extra", dims, 
                    numpy.asarray(extra_points.T, order="C"))

        for i_brick, brk in enumerate(pic.bricks):
            coords = [
                numpy.arange(
                    brk.origin[axis] + brk.stepwidths[axis]/2, 
                    brk.origin[axis] 
                    + brk.dimensions[axis] * brk.stepwidths[axis], 
                    brk.stepwidths[axis])
                for axis in xrange(dims)]
            for axis in xrange(dims):
                assert len(coords[axis]) == brk.dimensions[axis]

            mname = "structmesh%d" % i_brick
            silo.put_quadmesh(mname, coords)

            from pylo import DB_NODECENT

            brk_start = brk.start_index
            brk_stop = brk.start_index + len(brk)

            for quant in quantities:
                if quant == "rho":
                    vname = "rho_struct%d" % i_brick

                    silo.put_quadvar1(vname, mname, 
                            pic.get_grid_rho()[brk_start:brk_stop],
                            brk.dimensions, DB_NODECENT)
                elif quant == "j":
                    vname = "j_struct%d" % i_brick

                    vnames = [
                        "%s_coord%d" % (vname, axis) 
                        for axis in range(vdims)]

                    from pytools import product
                    j_grid = numpy.reshape(
                            pic.get_grid_j(self.cloud.velocities()),
                            (pic.grid_node_count(), vdims))[brk_start:brk_stop]

                    j_grid_compwise = numpy.asarray(j_grid.T, order="C")
                    
                    silo.put_quadvar(vname, mname, vnames,
                            j_grid_compwise,
                            brk.dimensions, DB_NODECENT)
                elif quant == "usecount":
                    vname = "usecount_struct%d" % i_brick
                    usecount = numpy.zeros((len(brk),), dtype=float)
                    for eog in pic.elements_on_grid:
                        for idx in eog.grid_nodes:
                            if brk_start <= idx < brk_stop:
                                usecount[idx-brk_start] += 1
                    silo.put_quadvar1(vname, mname, 
                            usecount, brk.dimensions, DB_NODECENT)
                else:
                    raise ValueError, "invalid vis quantity: %s" % quant
