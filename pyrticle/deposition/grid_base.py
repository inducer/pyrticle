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



import numpy




# grid deposition : brick generation --------------------------------------
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
        else:
            raise ValueError, "invalid dimensionality"




# grid visualization ----------------------------------------------------------
class GridVisualizer(object):
    def visualize_grid_quantities(self, silo, names_and_quantities):
        dims = self.method.dimensions_mesh
        backend = self.backend

        try:
            extra_points = backend.extra_points
        except AttributeError:
            extra_points = []

        if len(extra_points):
            extra_points = numpy.reshape(extra_points,
                    (len(extra_points)//dims,dims))

            silo.put_pointmesh("rec_grid_extra", 
                    numpy.asarray(extra_points.T, order="C"))

        try:
            grid_nodes_t = self.grid_nodes_t
        except AttributeError:
            grid_nodes = []
            from pyrticle._internal import BoxInt
            for brk in backend.bricks:
                grid_nodes.extend(brk.point(c) 
                        for c in self.iterator_type(brk, 
                            BoxInt(0*brk.dimensions, brk.dimensions)))
            grid_nodes_t = self.grid_nodes_t = \
                    numpy.asarray(numpy.array(grid_nodes).T, order="C")

        silo.put_pointmesh("rec_grid_nodes", grid_nodes_t)

        from pylo import DB_ZONECENT, DB_QUAD_RECT, DBObjectType

        if len(backend.bricks) > 1:
            def name_mesh(brick): return "_rec_grid_b%d" % brick
            def name_var(name, brick): return "_%s_b%d" % (name, brick)
        else:
            def name_mesh(brick): return "rec_grid"
            def name_var(name, brick): return name

        for brk in backend.bricks:
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

        if len(backend.bricks) > 1:
            silo.put_multimesh("rec_grid", 
                    [(name_mesh(brk.number), DB_QUAD_RECT) for brk in backend.bricks])

            for name, quant in names_and_quantities:
                silo.put_multivar(name, 
                        [(name_var(name, brk.number), DBObjectType.DB_QUADVAR) 
                            for brk in backend.bricks])





