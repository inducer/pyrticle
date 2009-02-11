"""Various mesh setups."""

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
import numpy.linalg as la




# 2D --------------------------------------------------------------------------
def make_glued_rect_mesh(a=(0,0), b=(1,1), max_area=None, 
        boundary_tagger=(lambda fvi, el, fn: []),
        periodicity=None, subdivisions=None,
        refine_func=None):
    """Create two unstructured rectangular meshes, side-by-side.

    @arg a: the lower left hand point of the base rectangle
    @arg b: the upper right hand point of the base rectangle
    @arg max_area: maximum area of each triangle.
    @arg periodicity: either None, or a tuple of bools specifying whether
      the mesh is to be periodic in x and y.
    @arg subdivisions: If not C{None}, this is a 2-tuple specifying
      the number of facet subdivisions in X and Y.
    @arg refine_func: A refinement function as taken by C{meshpy.triangle.build}.
    """
    import meshpy.triangle as triangle

    def round_trip_connect(start, end):
        for i in range(start, end):
            yield i, i+1
        yield end, start

    if max_area is not None:
        if refine_func is not None:
            raise ValueError, "cannot specify both refine_func and max_area"
        def refine_func(vertices, area):
            return area > max_area

    marker2tag = {
            1: "minus_x", 
            2: "minus_y", 
            3: "plus_x", 
            4: "plus_y", 
            }

    w = b[0]-a[0]
    h = b[1]-a[1]

    assert w > 0
    assert h > 0

    # 3----4----5
    # |    |    |
    # 0----1----2

    points = [
            a, 
            (a[0]+w,a[1]),
            (a[0]+2*w,a[1]),
            (a[0],a[1]+h),
            (a[0]+w,a[1]+h),
            (a[0]+2*w,a[1]+h),
            ]
    facets = [(0,1),(1,2),(2,5),(5,4),(4,3),(3,0), (1,4)]
    facet_markers = [2, 2, 3, 4, 4, 1, 0]

    if subdivisions is not None:
        points, facets, facet_markers = triangle.subdivide_facets(
                [subdivisions[0], subdivisions[0], 
                    subdivisions[1],
                    subdivisions[0], subdivisions[0], 
                    subdivisions[1],
                    subdivisions[1]],
                points, facets, facet_markers)
            
    from hedge.mesh import finish_2d_rect_mesh
    return finish_2d_rect_mesh(points, facets, facet_markers, marker2tag, 
            refine_func, periodicity, boundary_tagger)




def make_fine_center_rect_mesh(a=(0,0), b=(1,1), 
        boundary_tagger=(lambda fvi, el, fn: []),
        periodicity=None, 
        inner_width=0.1,
        max_area=None, inner_max_area=0.02,
        subdivisions=None, inner_subdivisions=None,
        refine_func=None):
    """Create an unstructured rectangular mesh.

    @arg a: the lower left hand point of the rectangle
    @arg b: the upper right hand point of the rectangle
    @arg periodicity: either None, or a tuple of bools specifying whether
      the mesh is to be periodic in x and y.
    @arg inner_width: total width of the refined inner "channel".
    @arg max_area: maximum area of each triangle outside the "channel".
    @arg inner_max_area: maximum area of each triangle inside the "channel".
    @arg subdivisions: If not C{None}, this is a 2-tuple specifying
      the number of facet subdivisions in X and Y.
    @arg inner_subdivisions: If not C{None}, this is a 2-tuple specifying
      the number of facet subdivisions in X and Y used on the channel boundary.
    @arg refine_func: A refinement function as taken by C{meshpy.triangle.build}.
    """
    import meshpy.triangle as triangle

    def round_trip_connect(start, end):
        for i in range(start, end):
            yield i, i+1
        yield end, start

    if max_area is not None:
        if refine_func is not None:
            raise ValueError, "cannot specify both refine_func and max_area"
        def refine_func(vertices, area):
            from pytools import average
            centroid_y = average(vertex[1] for vertex in vertices)
            if a[1]+fine_start_h <= centroid_y <= a[1]+fine_end_h:
                return area > inner_max_area
            else:
                return area > max_area

    marker2tag = {
            1: "minus_x", 
            2: "minus_y", 
            3: "plus_x", 
            4: "plus_y", 
            }

    w = b[0]-a[0]
    h = b[1]-a[1]

    assert w > 0
    assert h > 0
    assert inner_width > 0

    fine_start_h = h*0.5-inner_width*0.5
    fine_end_h = h*0.5+inner_width*0.5

    # 6-----------7 a[1]+h = b[1]
    # | coarse    |
    # 4-----------5 a[1]+fine_end_h
    # | fine      |
    # 2-----------3 a[1]+fine_start_h
    # | coarse    |
    # 0-----------1 a[1]

    a = numpy.asarray(a)
    b = numpy.asarray(b)
    x = numpy.array([1,0])
    y = numpy.array([0,1])

    points = [
            a, 
            a+w*x,
            a+fine_start_h*y,
            a+fine_start_h*y+w*x,
            a+fine_end_h*y,
            a+fine_end_h*y+w*x,
            a+h*y,
            a+h*y+w*x,
            ]

    facets = [
            (0,1),(0,2),(1,3),
            (2,3),(2,4),(3,5),
            (4,5),(4,6),(5,7),
            (6,7)
            ]
    facet_markers = [
            2,1,3,
            0,1,3,
            0,1,3,
            4
            ]

    if subdivisions is not None:
        assert inner_subdivisions is not None
        points, facets, facet_markers = triangle.subdivide_facets([
            subdivisions[0], subdivisions[1], subdivisions[1],
            inner_subdivisions[0], inner_subdivisions[1], inner_subdivisions[1],
            inner_subdivisions[0], subdivisions[1], subdivisions[1],
            subdivisions[0]],
            points, facets, facet_markers)

    from hedge.mesh import finish_2d_rect_mesh
    return finish_2d_rect_mesh(points, facets, facet_markers, marker2tag, refine_func,
            periodicity, boundary_tagger)





# 3D --------------------------------------------------------------------------
def make_extrusion_with_fine_core(rz, inner_r, 
        max_volume_inner=1e-4, max_volume_outer=5e-2,
        radial_subdiv=20):

    min_z = min(rz_pt[1] for rz_pt in rz)
    max_z = max(rz_pt[1] for rz_pt in rz)

    from meshpy.tet import MeshInfo, build
    from meshpy.geometry import generate_surface_of_revolution

    MINUS_Z_MARKER = 1
    PLUS_Z_MARKER = 2

    inner_points, inner_facets, inner_holes, inner_markers = \
            generate_surface_of_revolution(
                    [
                        (0, min_z),
                        (inner_r, min_z), 
                        (inner_r, max_z),
                        (0, max_z),
                        ],
                    ring_markers=[
                        MINUS_Z_MARKER,
                        0,
                        PLUS_Z_MARKER
                        ],
                    radial_subdiv=radial_subdiv,
                    )

    inner_point_indices = tuple(range(len(inner_points)))

    outer_points, outer_facets, outer_holes, outer_markers = \
            generate_surface_of_revolution(
                    [(inner_r,min_z)] + rz + [(inner_r, max_z)], 
                    ring_markers=[MINUS_Z_MARKER] + [0]*(len(rz)-1) + [PLUS_Z_MARKER],
                    point_idx_offset=len(inner_points),
                    radial_subdiv=radial_subdiv,
                    ring_point_indices= 
                    [ inner_point_indices[:radial_subdiv] ] 
                    + [None]*len(rz)
                    + [inner_point_indices[radial_subdiv:]]
                    )


    mesh_info = MeshInfo()
    mesh_info.set_points(inner_points + outer_points)
    mesh_info.set_facets_ex(
            inner_facets + outer_facets,
            inner_holes + outer_holes,
            inner_markers + outer_markers,
            )

    # set regional max. volume
    mesh_info.regions.resize(2)
    mesh_info.regions[0] = [0, 0, (max_z+min_z)/2, 0, 
            max_volume_inner]
    mesh_info.regions[1] = [inner_r+(rz[0][0]-inner_r)/2, 0, (max_z+min_z)/2, 0, 
            max_volume_outer]

    # add periodicity
    mesh_info.pbc_groups.resize(1)
    pbcg = mesh_info.pbc_groups[0]

    pbcg.facet_marker_1 = MINUS_Z_MARKER
    pbcg.facet_marker_2 = PLUS_Z_MARKER

    pbcg.set_transform(translation=[0,0,max_z-min_z])

    mesh = build(mesh_info, verbose=True, volume_constraints=True)
    #print "%d elements" % len(mesh.elements)
    #mesh.write_vtk("gun.vtk")

    fvi2fm = mesh.face_vertex_indices_to_face_marker
        
    def zper_boundary_tagger(fvi, el, fn):
        face_marker = fvi2fm[frozenset(fvi)]
        if face_marker == MINUS_Z_MARKER:
            return ["minus_z"]
        elif face_marker == PLUS_Z_MARKER:
            return ["plus_z"]
        else:
            return ["shell"]

    from hedge.mesh import make_conformal_mesh
    return make_conformal_mesh(mesh.points, mesh.elements,
            zper_boundary_tagger,
            periodicity=[None, None, ("minus_z", "plus_z")])




def make_cylinder_with_fine_core(r, inner_r, min_z, max_z, 
        max_volume_inner=1e-4, max_volume_outer=5e-2,
        radial_subdiv=20):

    return make_extrusion_with_fine_core(
            [ (r, min_z), (r, max_z), ], inner_r,
            max_volume_inner=max_volume_inner,
            max_volume_outer=max_volume_outer,
            radial_subdiv=radial_subdiv)
