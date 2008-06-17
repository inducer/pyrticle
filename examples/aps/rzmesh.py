def make_mesh_info(rz, radial_subdiv):
    from meshpy.tet import MeshInfo, EXT_OPEN, generate_surface_of_revolution

    points, facets, facet_holestarts, facet_markers = \
            generate_surface_of_revolution(rz, closure=EXT_OPEN,
                    radial_subdiv=radial_subdiv)

    mesh_info = MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets_ex(facets, facet_holestarts, facet_markers)

    return mesh_info




def bounding_box(tuples):
    assert tuples

    l = len(tuples[0])
    return [(min(t[i] for t in tuples), max(t[i] for t in tuples))
            for i in range(l)]



    
def make_inverse_mesh_info(rz, radial_subdiv):
    # chop off points with zero radius
    while rz[0][0] == 0:
        rz.pop(0)
    while rz[-1][0] == 0:
        rz.pop(-1)

    # construct outer cylinder
    ((min_r, max_r), (min_z, max_z)) = bounding_box(rz)

    if rz[0][1] < rz[-1][1]:
        # built in positive z direction
        rz.extend([
                (max_r+2, max_z),
                (max_r+2, min_z),
                ])
    else:
        rz.extend([
                (max_r+2, min_z),
                (max_r+2, max_z),
                ])

    from meshpy.tet import MeshInfo, EXT_CLOSED_IN_RZ, \
            generate_surface_of_revolution
            
    points, facets, facet_holestarts, facet_markers = \
            generate_surface_of_revolution(rz, closure=EXT_CLOSED_IN_RZ,
                    radial_subdiv=radial_subdiv)

    mesh_info = MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets_ex(facets, facet_holestarts, facet_markers)

    return mesh_info




def make_mesh_info_with_inner_tube(rz, tube_r, radial_subdiv):
    # chop off points with zero radius
    while rz[0][0] == 0:
        rz.pop(0)
    while rz[-1][0] == 0:
        rz.pop(-1)

    # construct outer cylinder
    first_z = rz[0][1]
    last_z = rz[-1][1]

    rz.insert(0, (tube_r, first_z))
    rz.append((tube_r, last_z))

    from meshpy.tet import MeshInfo, EXT_OPEN, generate_surface_of_revolution

    outer_points, outer_facets, outer_facet_holestarts, outer_facet_markers = \
            generate_surface_of_revolution(rz, 
                    closure=EXT_OPEN, radial_subdiv=radial_subdiv)

    outer_point_indices = tuple(range(len(outer_points)))

    inner_points, inner_facets, inner_facet_holestarts, inner_facet_markers = \
            generate_surface_of_revolution(
                    [(0,first_z), 
                        (tube_r, first_z), 
                        (tube_r, last_z),
                        (0, last_z)], 
                    point_idx_offset=len(outer_points),
                    radial_subdiv=radial_subdiv,
                    ring_point_indices=[
                        None,
                        outer_point_indices[:radial_subdiv],
                        outer_point_indices[-radial_subdiv:],
                        None,
                        ]
                    )

    points = outer_points + inner_points
    facets = outer_facets + inner_facets
    facet_holestarts = outer_facet_holestarts + inner_facet_holestarts
    facet_markers = outer_facet_markers + inner_facet_markers

    mesh_info = MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets_ex(facets, facet_holestarts, facet_markers)

    # set regional max. volume
    mesh_info.regions.resize(1)
    mesh_info.regions[0] = [0, 0,(first_z+last_z)/2, 0, 1e-4]

    return mesh_info
