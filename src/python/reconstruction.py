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
        self.log_constants = {}
    
    def initialize(self, cloud):
        self.cloud = cloud
        self.shape_function = None

    def set_shape_function(self, sf):
        self.shape_function = sf

    def add_instrumentation(self, mgr):
        mgr.set_constant("reconstructor", self.name)

        for key, value in self.log_constants.iteritems():
            mgr.set_constant(key, value)

        from pytools.log import IntervalTimer, EventCounter,\
                time_and_count_function

        self.reconstruct_timer = IntervalTimer(
                "t_reconstruct",
                "Time spent reconstructing")
        self.reconstruct_counter = EventCounter(
                "n_reconstruct",
                "Number of reconstructions")

        self.reconstruct_densities = time_and_count_function(
                self.reconstruct_densites,
                self.reconstruct_timer,
                self.reconstruct_counter,
                1+self.cloud.dimensions_velocity)

        self.reconstruct_j = time_and_count_function(
                self.reconstruct_j,
                self.reconstruct_timer,
                self.reconstruct_counter,
                self.cloud.dimensions_velocity)

        self.reconstruct_rho = time_and_count_function(
                self.reconstruct_rho,
                self.reconstruct_timer,
                self.reconstruct_counter)

        mgr.add_quantity(self.reconstruct_timer)
        mgr.add_quantity(self.reconstruct_counter)

    def clear_particles(self):
        pass

    def reconstruct_hook(self):
        if self.shape_function is None:
            raise RuntimeError, "shape function never set"

    def reconstruct_densites(self, velocities):
        self.reconstruct_hook()
        return self.cloud.pic_algorithm.reconstruct_densities(velocities)

    def reconstruct_j(self, velocities):
        self.reconstruct_hook()
        return self.cloud.pic_algorithm.reconstruct_j(velocities)

    def reconstruct_rho(self):
        self.reconstruct_hook()
        return self.cloud.pic_algorithm.reconstruct_rho()

    def upkeep(self):
        pass

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

        cloud.pic_algorithm.setup_advective_reconstructor(
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

    def upkeep(self):
        self.cloud.pic_algorithm.perform_reconstructor_upkeep()

    def rhs(self):
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

    def add_rhs(self, rhs):
        from pyrticle.tools import NumberShiftableVector
        self.cloud.pic_algorithm.apply_advective_particle_rhs(
                NumberShiftableVector.unwrap(rhs))




# grid reconstruction ---------------------------------------------------------
class SingleBrickGenerator(object):
    def __init__(self, overresolve=1.5, mesh_margin=0):
        self.overresolve = overresolve
        self.mesh_margin = mesh_margin

    def log_data(self, mgr):
        mgr.set_constant("rec_grid_brick_gen", self.__class__.__name__)
        mgr.set_constant("rec_grid_overresolve", self.overresolve)
        mgr.set_constant("rec_grid_mesh_margin", self.mesh_margin)

    def __call__(self, discr):
        from hedge.discretization import ones_on_volume
        mesh_volume = discr.integral(ones_on_volume(discr))
        dx =  (mesh_volume / len(discr)/ self.overresolve)**(1/discr.dimensions)

        mesh = discr.mesh
        bbox_min, bbox_max = mesh.bounding_box()

        bbox_min -= self.mesh_margin
        bbox_max += self.mesh_margin

        bbox_size = bbox_max-bbox_min
        dims = numpy.asarray(bbox_size/dx, dtype=numpy.int32)
        stepwidths = bbox_size/dims
        yield stepwidths, bbox_min, dims




class FineCoreBrickGenerator(object):
    def __init__(self, overresolve=1.5, mesh_margin=0, 
            core_axis=None, core_fraction=0.1, core_factor=2):
        self.overresolve = overresolve
        self.mesh_margin = mesh_margin
        self.core_axis = core_axis
        self.core_fraction = core_fraction
        self.core_factor = core_factor

        assert isinstance(core_factor, int)

    def log_data(self, mgr):
        mgr.set_constant("rec_grid_brick_gen", self.__class__.__name__)
        mgr.set_constant("rec_grid_overresolve", self.overresolve)
        mgr.set_constant("rec_grid_mesh_margin", self.mesh_margin)
        mgr.set_constant("rec_grid_core_fraction", self.core_fraction)
        mgr.set_constant("rec_grid_core_factor", self.core_factor)

    def __call__(self, discr):
        mesh = discr.mesh
        bbox_min, bbox_max = mesh.bounding_box()
        d = len(bbox_min)

        core_axis = self.core_axis
        if core_axis is None:
            core_axis = d-1

        # calculate outer bbox, as above
        from hedge.discretization import ones_on_volume
        mesh_volume = discr.integral(ones_on_volume(discr))
        dx =  (mesh_volume / len(discr)/ self.overresolve)**(1/discr.dimensions)
                
        bbox_min -= self.mesh_margin
        bbox_max += self.mesh_margin

        bbox_size = bbox_max-bbox_min
        dims = numpy.asarray(bbox_size/dx, dtype=numpy.int32)
        stepwidths = bbox_size/dims

        # calculate inner bbox
        core_margin_dims = (dims*(1-self.core_fraction)/2).astype(numpy.int32)
        core_margin_dims[core_axis] = 0
        core_dx = dx/self.core_factor
        core_min = bbox_min + stepwidths*core_margin_dims
        core_max = bbox_max - stepwidths*core_margin_dims
        core_size = core_max-core_min
        core_dims = numpy.asarray(core_size/core_dx, dtype=numpy.int32)
        core_stepwidths = core_size/core_dims

        # yield the core
        yield core_stepwidths, core_min, core_dims

        # yield the surrounding bricks
        from hedge.tools import unit_vector
        axis_vec = unit_vector(d, core_axis, dtype=numpy.int32)
        if d == 2:
            margin_dims = core_margin_dims + dims*axis_vec
            yield stepwidths, bbox_min, margin_dims
            yield -stepwidths, bbox_max, margin_dims
        elif d == 3:
            other_axes = set(range(d)) - set([core_axis])
            x, y = other_axes

            # ^ y
            # |
            # +----+----------+-----+
            # |    |     2    |     |
            # |    +----------+     |
            # | 1  |   core   |  1  |
            # |    +----------+     |
            # |    |     2    |     |
            # +----+----------+-----+--> x
            # 

            x_vec = unit_vector(d, x, dtype=numpy.int32)

            dims1 = dims.copy()
            dims1[x] = core_margin_dims[x]
            dims2 = dims.copy()
            dims2[x] = dims[x]-2*core_margin_dims[x]
            dims2[y] = core_margin_dims[y]

            yield stepwidths, bbox_min, dims1
            yield -stepwidths, bbox_max, dims1

            x_offset_2 = x_vec*core_margin_dims[x]*stepwidths[x]
            yield stepwidths, bbox_min+x_offset_2, dims2
            yield -stepwidths, bbox_max-x_offset_2, dims2

            from pytools import product
            print "surround pts: %d, core pts: %d" % (
                    product(dims2)*2+product(dims1)*2,
                    product(core_dims))
        else:
            raise ValueError, "invalid dimensionality"




class GridReconstructor(Reconstructor):
    name = "Grid"

    def __init__(self, brick_generator=SingleBrickGenerator(), 
            el_tolerance=0.2,
            max_extra_points=20,
            enforce_continuity=False,
            method="simplex_enlarge",
            filter_min_amplification=None,
            filter_order=None):
        Reconstructor.__init__(self)
        self.brick_generator = brick_generator
        self.el_tolerance = el_tolerance
        self.max_extra_points = max_extra_points
        self.enforce_continuity = enforce_continuity
        self.method = method

        self.filter_min_amplification = filter_min_amplification
        self.filter_order = filter_order

    def initialize(self, cloud):
        Reconstructor.initialize(self, cloud)

        discr = cloud.mesh_data.discr

        pic = self.cloud.pic_algorithm

        if self.enforce_continuity:
            self.prepare_average_groups()
        else:
            pic.average_group_starts.append(0)

        if self.filter_min_amplification is not None:
            from hedge.discretization import Filter, ExponentialFilterResponseFunction
            self.filter = Filter(discr, ExponentialFilterResponseFunction(
                    self.filter_min_amplification, self.filter_order))
        else:
            self.filter = None

        grid_nodes = []
        from pyrticle._internal import RecBrick, RecBrickIterator, BoxInt
        for i, (stepwidths, origin, dims) in enumerate(
                self.brick_generator(discr)):
            brk = RecBrick(i, pic.grid_node_count(), stepwidths, origin, dims)
            pic.bricks.append(brk)

            grid_nodes.extend(brk.point(c) 
                    for c in RecBrickIterator(brk, 
                        BoxInt(numpy.zeros((len(dims),), dtype=numpy.int32),
                            brk.dimensions)))
        self.grid_nodes = numpy.array(grid_nodes)

        if self.method == "simplex_extra":
            self.prepare_with_pointwise_projection_and_extra_points()
        elif self.method == "simplex_enlarge":
            self.prepare_with_pointwise_projection_and_enlargement()
        elif self.method == "simplex_reduce":
            self.prepare_with_pointwise_projection_and_basis_reduction()
        elif self.method == "brick":
            self.prepare_with_brick_interpolation()
        else:
            raise RuntimeError, "invalid rec_grid submethod specified"

    def set_shape_function(self, sf):
        Reconstructor.set_shape_function(self, sf)
        self.cloud.pic_algorithm.shape_function = sf

    def add_instrumentation(self, mgr):
        Reconstructor.add_instrumentation(self, mgr)

        mgr.set_constant("rec_grid_el_tolerance", self.el_tolerance)
        mgr.set_constant("rec_grid_enforce_continuity", self.enforce_continuity)
        mgr.set_constant("rec_grid_method", self.method)

        self.brick_generator.log_data(mgr)




    # preparation helpers -----------------------------------------------------
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
        pic = self.cloud.pic_algorithm

        from pyrticle._internal import ElementOnGrid
        eog = ElementOnGrid()
        eog.element_number = el.id

        points = pic.find_points_in_element(eog, el_tolerance * la.norm(el.map.matrix, 2))
        
        return eog, list(points)

    def scaled_vandermonde(self, el, eog, points, basis):
        """The remapping procedure in rec_grid relies on a pseudoinverse
        minimizing the pointwise error on all found interpolation points.
        But if an element spans bricks with different resolution, the 
        high-res points get an unfair advantage, because there's so many
        more of them. Therefore, each point is weighted by its 
        corresponding cell volume, meaning that the corresponding row of
        the structured Vandermonde matrix must be scaled by this volume, 
        as computed by L{find_points_in_element}. This routine returns the
        scaled Vandermonde matrix.
        """

        from hedge.polynomial import generic_vandermonde
        vdm = generic_vandermonde(
                [el.inverse_map(x) for x in points], 
                basis)

        for i, weight in enumerate(eog.weight_factors):
            vdm[i] *= weight

        return vdm

    def make_pointwise_interpolation_matrix(self, eog, eg, el, ldis, svd, scaled_vdm,
            basis_subset=None):
        u, s, vt = svd

        point_count = u.shape[0]
        node_count = vt.shape[1]

        thresh = (numpy.finfo(float).eps * max(s.shape) * s[0])

        nonzero_flags = numpy.abs(s) >= thresh
        inv_s = numpy.zeros((len(s),), dtype=float)
        inv_s[nonzero_flags] = 1/s[nonzero_flags]

        # compute the pseudoinverse of the structured Vandermonde matrix
        inv_s_diag = numpy.zeros(
                (node_count, point_count), 
                dtype=float)
        inv_s_diag[:len(s),:len(s)] = numpy.diag(1/s)

        svdm_pinv = numpy.dot(numpy.dot(vt.T, inv_s_diag), u.T)

        # check that it's reasonable
        pinv_resid = la.norm(
            numpy.dot(svdm_pinv, scaled_vdm)
            - numpy.eye(node_count))

        if pinv_resid > 1e-8:
            from warnings import warn
            warn("rec_grid: bad pseudoinv precision, element=%d, "
                    "#nodes=%d, #sgridpts=%d, resid=%.5g centroid=%s"
                % (el.id, node_count, point_count, pinv_resid,
                        el.centroid(self.cloud.mesh_data.discr.mesh.points)))

        el_vdm = ldis.vandermonde()
        if basis_subset is not None:
            el_vdm = el_vdm[:,basis_subset]
        imat = numpy.dot(el_vdm, svdm_pinv)

        if self.filter is not None:
            imat = numpy.dot(self.filter.get_filter_matrix(eg), imat)

        eog.interpolation_matrix = numpy.asarray(imat, order="F")

        if basis_subset is None:
            from hedge.tools import leftsolve
            eog.inverse_interpolation_matrix = numpy.asarray(
                    leftsolve(el_vdm, scaled_vdm), order="F")

    def generate_point_statistics(self, cond_claims=0):
        pic = self.cloud.pic_algorithm
        discr = self.cloud.mesh_data.discr

        point_claimers = {}

        for eog in pic.elements_on_grid:
            for i in eog.grid_nodes:
                point_claimers.setdefault(i, []).append(eog.element_number)

        claims = 0
        multiple_claims = 0
        
        for i_pt, claimers in point_claimers.iteritems():
            claims += len(claimers)
            if len(claimers) > 1:
                multiple_claims += 1

        claimed_pts = set(point_claimers.iterkeys())
        all_pts = set(xrange(pic.grid_node_count()))
        unclaimed_pts = all_pts-claimed_pts

        print("rec_grid.%s stats: #nodes: %d, #points total: %d" 
                % (self.method, len(discr), len(all_pts)))
        print("  #points unclaimed: %d #points claimed: %d, "
                "#points multiply claimed: %d, #claims: %d"
                % (len(unclaimed_pts), len(claimed_pts), multiple_claims, claims))
        print("  #claims made for conditioning: %d, #extra points: %d"
                % (cond_claims, len(pic.extra_points)/discr.dimensions))

        self.log_constants["rec_grid_points"] = len(all_pts)
        self.log_constants["rec_grid_claimed"] = len(claimed_pts)
        self.log_constants["rec_grid_claims"] = claims
        self.log_constants["rec_grid_mulclaims"] = multiple_claims
        self.log_constants["rec_grid_4cond"] = cond_claims
        self.log_constants["rec_grid_extra"] = len(pic.extra_points)/discr.dimensions




    # preparation methods -----------------------------------------------------
    def prepare_with_pointwise_projection_and_extra_points(self):
        discr = self.cloud.mesh_data.discr
        pic = self.cloud.pic_algorithm

        pic.elements_on_grid.reserve(
                sum(len(eg.members) for eg in discr.element_groups))

        # map brick numbers to [ (point, el_id, el_structured_point_index),...]
        # This is used to write out the C++ extra_points structure further
        # down.
        ep_brick_map = {}

        min_s_values = []

        # Iterate over all elements
        for eg in discr.element_groups:
            ldis = eg.local_discretization

            for el in eg.members:
                eog, points = self.find_points_in_element(el, self.el_tolerance)

                # If the structured Vandermonde matrix is singular,
                # add "extra points" to prevent that.
                ep_count = 0
                while True:
                    scaled_vdm = self.scaled_vandermonde(el, eog, points, 
                            ldis.basis_functions())

                    u, s, vt = svd = la.svd(scaled_vdm)

                    if len(points) >= ldis.node_count():
                        # case 1: theoretically enough points found
                        thresh = (numpy.finfo(float).eps
                                * max(scaled_vdm.shape) * s[0])
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
                                "vandermonde matrix for el #%d with #ep bound" % el.id)
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
                    new_point = discr.nodes[
                            discr.find_el_range(el.id).start
                            +max_node_idx]

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
                min_s_values.append(min(s))

                self.make_pointwise_interpolation_matrix(eog, eg, el, ldis, svd, scaled_vdm)

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
        self.generate_point_statistics(len(extra_points))




    def prepare_with_pointwise_projection_and_enlargement(self):
        tolerance_bound = 1.5

        discr = self.cloud.mesh_data.discr
        pic = self.cloud.pic_algorithm

        pic.elements_on_grid.reserve(
                sum(len(eg.members) for eg in discr.element_groups))

        cond_claims = 0
        min_s_values = []
        max_s_values = []
        cond_s_values = []

        from hedge.discretization import Projector, Discretization
        fine_discr = Discretization(discr.mesh, order=8)
        proj = Projector(discr, fine_discr)

        from hedge.visualization import SiloVisualizer
        vis = SiloVisualizer(fine_discr)

        # Iterate over all elements
        for eg in discr.element_groups:
            ldis = eg.local_discretization

            for el in eg.members:
                # If the structured Vandermonde matrix is singular,
                # enlarge the element tolerance

                my_tolerance = self.el_tolerance

                orig_point_count = None

                while True:
                    eog, points = self.find_points_in_element(el, my_tolerance)
                    if orig_point_count is None:
                        orig_point_count = len(points)

                    scaled_vdm = self.scaled_vandermonde(el, eog, points, 
                            ldis.basis_functions())

                    bad_vdm = len(points) < ldis.node_count()
                    if not bad_vdm:
                        try:
                            u, s, vt = svd = la.svd(scaled_vdm)
                            thresh = (numpy.finfo(float).eps
                                    * max(scaled_vdm.shape) * s[0])
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
                                "vandermonde matrix for el #%d by enlargement" % el.id)
                        break

                #from pytools import average
                #print average(eog.weight_factors), min(s)
                #raw_input()

                if my_tolerance > self.el_tolerance:
                    print "element %d: #nodes=%d, orig #sgridpt=%d, extra tol=%g, #extra points=%d" % (
                            el.id, ldis.node_count(), orig_point_count, 
                            my_tolerance-self.el_tolerance,
                            len(points)-orig_point_count)
                    cond_claims += len(points)-orig_point_count
                min_s_values.append(min(s))
                max_s_values.append(max(s))
                cond_s_values.append(max(s)/min(s))

                if max(s)/min(s) > 1e2:
                    for i in range(len(s)):
                        if s[0]/s[i] > 1e2:
                            zeroed_mode = vt[i]
                            zeroed_mode_nodal = numpy.dot(ldis.vandermonde(), 
                                    zeroed_mode)
                            print el.id, i, s[0]/s[i]
                            fromvec = discr.volume_zeros()
                            fromvec[discr.find_el_range(el.id)] = zeroed_mode_nodal

                            gn = list(eog.grid_nodes)
                            assert len(gn) == len(u[i])

                            tovec = numpy.zeros((pic.grid_node_count(),), dtype=float)
                            tovec[gn] = s[i]*u[i]

                            usevec = numpy.zeros((pic.grid_node_count(),), dtype=float)
                            usevec[gn] = 1

                            visf = vis.make_file("nulled-%04d%02d" % (el.id, i))
                            vis.add_data(visf, [("meshmode", proj(fromvec))],
                                    expressions=[
                                    ("absmesh", "abs(meshmode)"),
                                    ("absgrid", "abs(gridmode)"),
                                    ])
                            self.visualize_grid_quantities(visf, [
                                    ("gridmode", tovec),
                                    ("usevec", usevec),
                                    ],
                                    )
                            visf.close()
                    #print s
                    #raw_input()

                self.make_pointwise_interpolation_matrix(eog, eg, el, ldis, svd, scaled_vdm)

                pic.elements_on_grid.append(eog)

        # we don't need no stinkin' extra points
        pic.extra_point_brick_starts.extend([0]*(len(pic.bricks)+1))

        # print some statistics
        self.generate_point_statistics(cond_claims)




    def prepare_with_pointwise_projection_and_basis_reduction(self):
        discr = self.cloud.mesh_data.discr
        pic = self.cloud.pic_algorithm

        pic.elements_on_grid.reserve(
                sum(len(eg.members) for eg in discr.element_groups))

        min_s_values = []
        max_s_values = []
        cond_s_values = []

        basis_len_vec = discr.volume_zeros()
        el_condition_vec = discr.volume_zeros()
        point_count_vec = discr.volume_zeros()

        # Iterate over all elements
        for eg in discr.element_groups:
            ldis = eg.local_discretization

            mode_id_to_index = dict(
                    (bid, i) for i, bid in enumerate(ldis.generate_mode_identifiers()))

            for el in eg.members:
                basis = list(zip(
                    ldis.generate_mode_identifiers(), 
                    ldis.basis_functions()))

                eog, points = self.find_points_in_element(el, self.el_tolerance)

                while True:
                    scaled_vdm = self.scaled_vandermonde(el, eog, points, 
                            [bf for bid, bf in basis])

                    max_bid_sum = max(sum(bid) for bid, bf in basis)
                    killable_basis_elements = [
                            (i, bid) for i, (bid, bf) in enumerate(basis)
                            if sum(bid) == max_bid_sum]

                    try:
                        u, s, vt = svd = la.svd(scaled_vdm)

                        thresh = (numpy.finfo(float).eps
                                * max(scaled_vdm.shape) * s[0])

                        assert s[-1] == numpy.min(s)
                        assert s[0] == numpy.max(s)
                        
                        if len(basis) > len(points) or s[0]/s[-1] > 10:
                            retry = True

                            # badly conditioned, kill a basis entry
                            vti = vt[-1]

                            from pytools import argmax2
                            kill_idx, kill_bid = argmax2(
                                    ((j, bid), abs(vti[j])) 
                                    for j, bid in killable_basis_elements)

                            assert kill_bid == basis[kill_idx][0]
                            basis.pop(kill_idx)
                        else:
                            retry = False

                    except la.LinAlgError:
                        # SVD likely didn't converge. Lacking an SVD, we don't have 
                        # much guidance on what modes to kill. Any of the killable
                        # ones will do.

                        # Bang, you're dead.
                        basis.pop(killable_basis_elements[0][0])

                        retry = True

                    if not retry:
                        break

                    if len(basis) == 1:
                        raise RuntimeError(
                                "basis reduction has killed almost the entire basis on element %d"
                                % el.id)

                print "element %d: #nodes=%d, leftover modes=%d" % (
                        el.id, ldis.node_count(), len(basis),)

                basis_len_vec[discr.find_el_range(el.id)] = len(basis)
                el_condition_vec[discr.find_el_range(el.id)] = s[0]/s[-1]
                point_count_vec[discr.find_el_range(el.id)] = len(points)

                min_s_values.append(min(s))
                max_s_values.append(max(s))
                cond_s_values.append(max(s)/min(s))

                self.make_pointwise_interpolation_matrix(eog, eg, el, ldis, svd, scaled_vdm,
                        basis_subset=[mode_id_to_index[bid] for bid, bf in basis])

                pic.elements_on_grid.append(eog)


        # visualize basis length for each element
        from hedge.visualization import SiloVisualizer
        vis = SiloVisualizer(discr)
        visf = vis.make_file("rec-debug")
        vis.add_data(visf, [
            ("basis_len", basis_len_vec),
            ("el_condition", el_condition_vec),
            ("point_count", point_count_vec),
            ])
        visf.close()

        # we don't need no stinkin' extra points
        pic.extra_point_brick_starts.extend([0]*(len(pic.bricks)+1))

        # print some statistics
        self.generate_point_statistics()




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

        from pyrticle._internal import RecBrickIterator, ElementOnGrid, BoxFloat

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
                            for c in RecBrickIterator(brk, idx_range)]
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
                            for c in RecBrickIterator(brk, idx_range))
                    pic.elements_on_grid.append(eog)

        # we don't need no stinkin' extra points
        pic.extra_point_brick_starts.extend([0]*(len(pic.bricks)+1))

        # stats
        self.generate_point_statistics(0)





    # reconstruction onto mesh ------------------------------------------------
    def remap_grid_to_mesh(self, q_grid):
        discr = self.cloud.mesh_data.discr

        eff_shape = q_grid.shape[1:]
        if len(eff_shape) == 0:
            result = discr.volume_zeros()
            self.cloud.pic_algorithm.remap_grid_to_mesh(
                    q_grid, result, 0, 1)
        elif len(eff_shape) == 1:
            result = numpy.zeros((len(discr),)+eff_shape, dtype=float)
            for i in range(eff_shape[0]):
                self.cloud.pic_algorithm.remap_grid_to_mesh(
                        q_grid, result, i, eff_shape[0])
        else:
            raise ValueError, "invalid effective shape for remap"
        return result

    def reconstruct_densites(self, velocities):
        return tuple(
                self.remap_grid_to_mesh(q_grid) 
                for q_grid in self.reconstruct_grid_densities(velocities))

    def reconstruct_j(self, velocities):
        return self.remap_grid_to_mesh(self.reconstruct_grid_j(velocities))

    def reconstruct_rho(self):
        return self.remap_grid_to_mesh(self.reconstruct_grid_rho())

    # reconstruction onto grid ------------------------------------------------
    def reconstruct_grid_densities(self, velocities):
        self.reconstruct_hook()
        return self.cloud._get_derived_quantities_from_cache(
                ["rho_grid", "j_grid"],
                [self.reconstruct_grid_rho, self.reconstruct_grid_j],
                lambda: self.cloud.pic_algorithm.reconstruct_grid_densities(velocities))

    def reconstruct_grid_j(self, velocities):
        self.reconstruct_hook()
        return self.cloud._get_derived_quantity_from_cache("j_grid", 
                lambda: self.cloud.pic_algorithm.reconstruct_grid_j(velocities))

    def reconstruct_grid_rho(self):
        self.reconstruct_hook()
        return self.cloud._get_derived_quantity_from_cache("rho_grid", 
                self.cloud.pic_algorithm.reconstruct_grid_rho)

    # grid debug quantities ---------------------------------------------------
    def ones_on_grid(self):
        return numpy.ones((self.cloud.pic_algorithm.grid_node_count(),), 
                dtype=float)

    def grid_usecount(self):
        usecount = numpy.zeros((self.cloud.pic_algorithm.grid_node_count(),), 
                dtype=float)
        for eog in self.cloud.pic_algorithm.elements_on_grid:
            for idx in eog.grid_nodes:
                usecount[idx] += 1
        return usecount

    def remap_residual(self, q_grid):
        discr = self.cloud.mesh_data.discr

        gnc = self.cloud.pic_algorithm.grid_node_count()

        eff_shape = q_grid.shape[1:]
        if len(eff_shape) == 0:
            result = numpy.zeros((gnc,), dtype=float)
            self.cloud.pic_algorithm.remap_residual(
                    q_grid, result, 0, 1)
        elif len(eff_shape) == 1:
            result = numpy.zeros((gnc,)+eff_shape, dtype=float)
            for i in range(eff_shape[0]):
                self.cloud.pic_algorithm.remap_residual(
                        q_grid, result, i, eff_shape[0])
        else:
            raise ValueError, "invalid effective shape for remap"
        return result




    # grid visualization ------------------------------------------------------
    def visualize_grid_quantities(self, silo, names_and_quantities):
        dims = self.cloud.dimensions_mesh
        vdims = self.cloud.dimensions_velocity
        pic = self.cloud.pic_algorithm

        extra_points = pic.extra_points
        if len(extra_points):
            extra_points = numpy.reshape(extra_points,
                    (len(extra_points)//dims,dims))

            silo.put_pointmesh("rec_grid_extra", 
                    numpy.asarray(extra_points.T, order="C"))

        silo.put_pointmesh("rec_grid_nodes", 
                numpy.asarray(self.grid_nodes.T, order="C"))

        from pylo import DB_ZONECENT, DB_QUAD_RECT, DBObjectType

        if len(pic.bricks) > 1:
            def name_mesh(brick): return "_rec_grid_b%d" % brick
            def name_var(name, brick): return "_%s_b%d" % (name, brick)
        else:
            def name_mesh(brick): return "rec_grid"
            def name_var(name, brick): return name

        for brk in pic.bricks:
            coords = [
                numpy.arange(
                    brk.origin[axis], 
                    brk.origin[axis] 
                    + brk.dimensions[axis] * brk.stepwidths[axis] 
                    + brk.stepwidths[axis]/2, 
                    brk.stepwidths[axis])
                for axis in xrange(dims)]
            for axis in xrange(dims):
                assert len(coords[axis]) == brk.dimensions[axis]+1

            mname = name_mesh(brk.number)
            silo.put_quadmesh(mname, coords)

            brk_start = brk.start_index
            brk_stop = brk.start_index + len(brk)

            for name, quant in names_and_quantities:
                eff_shape = quant.shape[1:]
                vname = name_var(name, brk.number)
                if len(eff_shape) == 0:
                    silo.put_quadvar1(vname, mname, 
                            quant[brk_start:brk_stop], brk.dimensions, DB_ZONECENT)
                elif len(eff_shape) == 1:
                    d = eff_shape[0]
                    silo.put_quadvar(vname, mname, 
                            ["%s_c%d" % (vname, axis) for axis in range(d)],
                            numpy.asarray(quant[brk_start:brk_stop].T, order="C"),
                            brk.dimensions, DB_ZONECENT)
                else:
                    raise ValueError, "invalid effective shape for vis"

        if len(pic.bricks) > 1:
            silo.put_multimesh("rec_grid", 
                    [(name_mesh(brk.number), DB_QUAD_RECT) for brk in pic.bricks])

            for name, quant in names_and_quantities:
                silo.put_multivar(name, 
                        [(name_var(name, brk.number), DBObjectType.DB_QUADVAR) 
                            for brk in pic.bricks])
