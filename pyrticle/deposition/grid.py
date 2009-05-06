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
from pyrticle.deposition.grid_base import GridVisualizer, SingleBrickGenerator
import pyrticle._internal as _internal
import numpy
import numpy.linalg as la




# pure grid deposition ----------------------------------------------------
class GridDepositor(Depositor, GridVisualizer):
    def __init__(self, brick_generator=SingleBrickGenerator(), 
            el_tolerance=0.12,
            max_extra_points=20,
            enforce_continuity=False,
            submethod="simplex_reduce",
            filter_min_amplification=None,
            filter_order=None,
            jiggle_radius=0.0,
            ):
        Depositor.__init__(self)
        self.brick_generator = brick_generator
        self.el_tolerance = el_tolerance
        self.max_extra_points = max_extra_points
        self.enforce_continuity = enforce_continuity
        self.submethod = submethod

        self.filter_min_amplification = filter_min_amplification
        self.filter_order = filter_order

        self.jiggle_radius = jiggle_radius

    @property
    def name(self):
        if self.jiggle_radius:
            return "GridJiggly"
        else:
            return "GridRegular"

    @property
    def iterator_type(self):
        if self.jiggle_radius:
            return _internal.JigglyBrickIterator
        else:
            return _internal.BrickIterator

    def initialize(self, method):
        Depositor.initialize(self, method)

        if self.jiggle_radius:
            dep_type = "Jiggly"
        else:
            dep_type = "Regular"

        backend_class = getattr(_internal, 
                dep_type + "GridDepositor" + method.get_dimensionality_suffix())
        backend = self.backend = backend_class(method.mesh_data)

        discr = method.discretization

        if self.enforce_continuity:
            self.prepare_average_groups()
        else:
            backend.average_group_starts.append(0)

        if self.filter_min_amplification is not None:
            from hedge.discretization import Filter, ExponentialFilterResponseFunction
            self.filter = Filter(discr, ExponentialFilterResponseFunction(
                    self.filter_min_amplification, self.filter_order))
        else:
            self.filter = None

        for i, (stepwidths, origin, dims) in enumerate(
                self.brick_generator(discr)):
            if self.jiggle_radius:
                brk = _internal.JigglyBrick(i, backend.grid_node_count_with_extra(), 
                        stepwidths, origin, dims,
                        jiggle_radius=self.jiggle_radius)
            else:
                brk = _internal.Brick(i, backend.grid_node_count_with_extra(), 
                        stepwidths, origin, dims)
            backend.bricks.append(brk)

        if self.submethod == "simplex_extra":
            self.prepare_with_pointwise_projection_and_extra_points()
        elif self.submethod == "simplex_enlarge":
            self.prepare_with_pointwise_projection_and_enlargement()
        elif self.submethod == "simplex_reduce":
            self.prepare_with_pointwise_projection_and_basis_reduction()
        elif self.submethod == "brick":
            self.prepare_with_brick_interpolation()
        else:
            raise RuntimeError, "invalid rec_grid submethod specified"

    def set_shape_function(self, state, sf):
        Depositor.set_shape_function(self, state, sf)
        self.backend.shape_function = sf

    def add_instrumentation(self, mgr):
        Depositor.add_instrumentation(self, mgr)

        mgr.set_constant("rec_grid_el_tolerance", self.el_tolerance)
        mgr.set_constant("rec_grid_enforce_continuity", self.enforce_continuity)
        mgr.set_constant("rec_grid_method", self.method)
        mgr.set_constant("rec_grid_jiggle_radius", self.jiggle_radius)

        self.brick_generator.log_data(mgr)




    # preparation helpers -----------------------------------------------------
    def find_containing_brick(self, pt):
        for brk in self.backend.bricks:
            if brk.bounding_box().contains(pt):
                return brk
        raise RuntimeError, "no containing brick found for point"

    def prepare_average_groups(self):
        discr = self.method.discretization

        avg_group_finder = {}
        avg_groups = []

        for fg in discr.face_groups:
            for fp in fg.face_pairs:
                for el_idx, opp_el_idx in zip(
                        fg.index_lists[fp.loc.face_index_list_number],
                        fg.index_lists[fp.opp.face_index_list_number]):
                    idx1 = fp.loc.el_base_index + el_idx
                    idx2 = fp.opp.el_base_index + opp_el_idx

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
            self.backend.average_groups.extend(ag)
            self.backend.average_group_starts.append(
                    len(self.backend.average_groups))

        print len(avg_groups), "average groups"

    def find_points_in_element(self, el, el_tolerance):
        from pyrticle._internal import ElementOnGrid
        eog = ElementOnGrid()
        eog.element_number = el.id

        points = self.backend.find_points_in_element(
                eog, el_tolerance * la.norm(el.map.matrix, 2))
        
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
                        el.centroid(self.method.discretization.mesh.points)))

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
        discr = self.method.discretization

        point_claimers = {}

        for eog in self.backend.elements_on_grid:
            for i in eog.grid_nodes:
                point_claimers.setdefault(i, []).append(eog.element_number)

        claims = 0
        multiple_claims = 0
        
        for i_pt, claimers in point_claimers.iteritems():
            claims += len(claimers)
            if len(claimers) > 1:
                multiple_claims += 1

        claimed_pts = set(point_claimers.iterkeys())
        all_pts = set(xrange(self.backend.grid_node_count_with_extra()))
        unclaimed_pts = all_pts-claimed_pts

        print("rec_grid.%s stats: #nodes: %d, #points total: %d" 
                % (self.submethod, len(discr), len(all_pts)))
        print("  #points unclaimed: %d #points claimed: %d, "
                "#points multiply claimed: %d, #claims: %d"
                % (len(unclaimed_pts), len(claimed_pts), multiple_claims, claims))
        print("  #claims made for conditioning: %d, #extra points: %d"
                % (cond_claims, len(self.backend.extra_points)/discr.dimensions))

        self.log_constants["rec_grid_points"] = len(all_pts)
        self.log_constants["rec_grid_claimed"] = len(claimed_pts)
        self.log_constants["rec_grid_claims"] = claims
        self.log_constants["rec_grid_mulclaims"] = multiple_claims
        self.log_constants["rec_grid_4cond"] = cond_claims
        self.log_constants["rec_grid_extra"] = len(self.backend.extra_points)/discr.dimensions




    # preparation methods -----------------------------------------------------
    def prepare_with_pointwise_projection_and_extra_points(self):
        discr = self.method.discretization
        backend = self.backend

        backend.elements_on_grid.reserve(
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

                backend.elements_on_grid.append(eog)

        # fill in the extra points
        ep_brick_starts = [0]
        extra_points = []
        
        gnc = backend.grid_node_count_with_extra()
        for brk in backend.bricks:
            for pt, el_id, struc_idx in ep_brick_map.get(brk.number, []):
                # replace zero placeholder from above
                backend.elements_on_grid[el_id].grid_nodes[struc_idx] = \
                        gnc + len(extra_points)
                extra_points.append(pt)

            ep_brick_starts.append(len(extra_points))

        backend.first_extra_point = gnc
        backend.extra_point_brick_starts.extend(ep_brick_starts)
        backend.extra_points = numpy.array(extra_points)

        # print some statistics
        self.generate_point_statistics(len(extra_points))




    def prepare_with_pointwise_projection_and_enlargement(self):
        tolerance_bound = 1.5

        discr = self.method.discretization
        backend = self.backend

        backend.elements_on_grid.reserve(
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

                            tovec = numpy.zeros((backend.grid_node_count_with_extra(),), dtype=float)
                            tovec[gn] = s[i]*u[i]

                            usevec = numpy.zeros((backend.grid_node_count_with_extra(),), dtype=float)
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

                backend.elements_on_grid.append(eog)

        # we don't need no stinkin' extra points
        backend.extra_point_brick_starts.extend([0]*(len(backend.bricks)+1))

        # print some statistics
        self.generate_point_statistics(cond_claims)




    def prepare_with_pointwise_projection_and_basis_reduction(self):
        discr = self.method.discretization
        backend = self.backend

        backend.elements_on_grid.reserve(
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

                if ldis.node_count() > len(basis):
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

                backend.elements_on_grid.append(eog)

        # visualize basis length for each element
        if set(["depositor", "vis_files"]) < self.method.debug:
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
        backend.extra_point_brick_starts.extend([0]*(len(backend.bricks)+1))

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

        discr = self.method.discretization
        backend = self.backend

        backend.elements_on_grid.reserve(
                sum(len(eg.members) for eg in discr.element_groups))

        from pyrticle._internal import ElementOnGrid, BoxFloat

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
                for brk in backend.bricks:
                    eog = ElementOnGrid()
                    eog.element_number = el.id

                    brk_bbox = brk.bounding_box()
                    brk_and_el = brk_bbox.intersect(el_bbox)

                    if brk_and_el.is_empty():
                        continue

                    el_nodes = []
                    el_node_indices = []
                    el_slice = discr.find_el_range(el.id)
                    el_length = el_slice.stop-el_slice.start
                    if brk_and_el == el_bbox:
                        # whole element in brick? fantastic.
                        el_nodes = discr.nodes[el_slice]
                        el_node_indices = range(el_length)
                    else:
                        # no? go through the nodes one by one.
                        for i, node in enumerate(discr.nodes[el_slice]):
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
                            for c in self.iterator_type(brk, idx_range)]
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
                            for c in self.iterator_type(brk, idx_range))
                    eog.weight_factors = numpy.ones(
                            (len(eog.grid_nodes),), dtype=numpy.float64)
                    backend.elements_on_grid.append(eog)

        # we don't need no stinkin' extra points
        backend.extra_point_brick_starts.extend([0]*(len(backend.bricks)+1))

        # stats
        self.generate_point_statistics(0)





    # deposition onto mesh ------------------------------------------------
    def remap_grid_to_mesh(self, q_grid):
        discr = self.method.discretization

        eff_shape = q_grid.shape[1:]
        if len(eff_shape) == 0:
            result = discr.volume_zeros()
            self.backend.remap_grid_to_mesh(
                    q_grid, result, 0, 1)
        elif len(eff_shape) == 1:
            result = numpy.zeros((len(discr),)+eff_shape, dtype=float)
            for i in range(eff_shape[0]):
                self.backend.remap_grid_to_mesh(
                        q_grid, result, i, eff_shape[0])
        else:
            raise ValueError, "invalid effective shape for remap"
        return result

    def _deposit_densites(self, state, velocities, pslice):
        return tuple(
                self.remap_grid_to_mesh(q_grid) 
                for q_grid in self.deposit_grid_densities(
                    state, velocities, pslice))

    def _deposit_j(self, state, velocities, pslice):
        return self.remap_grid_to_mesh(self.deposit_grid_j(
            state, velocities, pslice))

    def _deposit_rho(self, state, pslice):
        return self.remap_grid_to_mesh(
                self.deposit_grid_rho(state, pslice))

    # deposition onto grid ------------------------------------------------
    def deposit_grid_densities(self, state, velocities, pslice=slice(None)):
        self.deposit_hook()
        return state.get_derived_quantities_from_cache(
                [("rho_grid", pslice.start, pslice.stop, pslice.step),
                    ("j_grid", pslice.start, pslice.stop, pslice.step)],
                [lambda: self.deposit_grid_rho(pslice), 
                    lambda: self.deposit_grid_j(state, velocities, pslice)],
                lambda: self.backend.deposit_grid_densities(
                    state.depositor_state, state.particle_state, 
                    velocities, pslice))

    def deposit_grid_j(self, state, velocities, pslice=slice(None)):
        self.deposit_hook()
        return state.get_derived_quantity_from_cache(
                ("j_grid", pslice.start, pslice.stop, pslice.step), 
                lambda: self.backend.deposit_grid_j(
                    state.depositor_state, state.particle_state, 
                    velocities, pslice))

    def deposit_grid_rho(self, state, pslice=slice(None)):
        self.deposit_hook()
        return state.get_derived_quantity_from_cache(
                ("rho_grid", pslice.start, pslice.stop, pslice.step), 
                lambda: self.backend.deposit_grid_rho(
                    state.depositor_state, state.particle_state, 
                    pslice))

    # grid debug quantities ---------------------------------------------------
    def ones_on_grid(self):
        return numpy.ones((self.backend.grid_node_count_with_extra(),), 
                dtype=float)

    def grid_usecount(self):
        usecount = numpy.zeros((self.backend.grid_node_count_with_extra(),), 
                dtype=float)
        for eog in self.backend.elements_on_grid:
            for idx in eog.grid_nodes:
                usecount[idx] += 1
        return usecount

    def remap_residual(self, q_grid):
        discr = self.method.discretization

        gnc = self.backend.grid_node_count_with_extra()

        eff_shape = q_grid.shape[1:]
        if len(eff_shape) == 0:
            result = numpy.zeros((gnc,), dtype=float)
            self.backend.remap_residual(
                    q_grid, result, 0, 1)
        elif len(eff_shape) == 1:
            result = numpy.zeros((gnc,)+eff_shape, dtype=float)
            for i in range(eff_shape[0]):
                self.backend.remap_residual(
                        q_grid, result, i, eff_shape[0])
        else:
            raise ValueError, "invalid effective shape for remap"
        return result






