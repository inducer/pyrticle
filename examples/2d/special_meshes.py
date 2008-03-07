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
        def refine_func(vert_origin, vert_destination, vert_apex, area):
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
            
    mesh_info = triangle.MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets(facets, facet_markers)
    
    if periodicity is None:
        periodicity = (False, False)

    axes = ["x", "y"]
    mesh_periodicity = []
    periodic_tags = set()
    for i, axis in enumerate(axes):
        if periodicity[i]:
            minus_tag = "minus_"+axis
            plus_tag = "plus_"+axis
            mesh_periodicity.append((minus_tag, plus_tag))
            periodic_tags.add(minus_tag)
            periodic_tags.add(plus_tag)
        else:
            mesh_periodicity.append(None)

    generated_mesh = triangle.build(mesh_info, 
            refinement_func=refine_func,
            allow_boundary_steiner=not (periodicity[0] or periodicity[1]))

    fvi2fm = dict((frozenset(fvi), marker) for fvi, marker in
        zip(generated_mesh.facets, generated_mesh.facet_markers))

    def wrapped_boundary_tagger(fvi, el, fn):
        btag = marker2tag[fvi2fm[frozenset(fvi)]]
        if btag in periodic_tags:
            return [btag]
        else:
            return [btag] + boundary_tagger(fvi, el, fn)

    from hedge.mesh import make_conformal_mesh
    return make_conformal_mesh(
            generated_mesh.points,
            generated_mesh.elements,
            wrapped_boundary_tagger,
            periodicity=mesh_periodicity)

