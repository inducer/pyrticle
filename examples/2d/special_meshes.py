import pylinear.array as num




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

    a = num.asarray(a)
    b = num.asarray(b)
    x = num.array([1,0])
    y = num.array([0,1])

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
