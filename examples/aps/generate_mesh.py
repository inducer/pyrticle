from meshpy.tet import MeshInfo, build, \
        generate_surface_of_revolution, \
        EXT_CLOSED_IN_RZ, EXT_OPEN




def bounding_box(tuples):
    assert tuples

    l = len(tuples[0])
    return [(min(t[i] for t in tuples), max(t[i] for t in tuples))
            for i in range(l)]



    
def main():
    from superfish import parse_superfish_format
    from optparse import OptionParser

    description = "generate a vtk mesh from a superfish input file"
    parser = OptionParser(description=description)
    parser.add_option(
	    "--point-dist", dest="point_dist", default="0.3",
	    help="Spacing of intermediate points in circular arcs")
    parser.add_option(
	    "--inverse", dest="generate_inverse", action="store_true",
	    help="Whether to mesh the inverse of the gun (useless, but fun)")
    parser.add_option(
	    "--radial-subdiv", dest="radial_subdiv", action="store", type="int",
            default=16, help="How many radial subdivisions")
    parser.add_option(
	    "--tube-radius", dest="tube_radius", action="store", type="float",
            default=0.1, help="Radius of the interior tube")

    options, args = parser.parse_args()

    mesh_info = MeshInfo()
    rz = parse_superfish_format(args[0], float(options.point_dist))

    closure = EXT_OPEN

    build_kwargs = {}

    if options.generate_inverse:
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
        closure=EXT_CLOSED_IN_RZ

    if options.tube_radius:
        assert closure is EXT_OPEN

        # chop off points with zero radius
        while rz[0][0] == 0:
            rz.pop(0)
        while rz[-1][0] == 0:
            rz.pop(-1)

        # construct outer cylinder
        first_z = rz[0][1]
        last_z = rz[-1][1]
        tube_r = options.tube_radius
        assert tube_r is not None

        rz.insert(0, (tube_r, first_z))
        rz.append((tube_r, last_z))

        outer_points, outer_facets, outer_facet_holestarts, outer_facet_markers = \
                generate_surface_of_revolution(rz, 
                        closure=closure,
                        radial_subdiv=options.radial_subdiv
                        )

        outer_point_indices = tuple(range(len(outer_points)))

        inner_points, inner_facets, inner_facet_holestarts, inner_facet_markers = \
                generate_surface_of_revolution(
                        [(0,first_z), 
                            (tube_r, first_z), 
                            (tube_r, last_z),
                            (0, last_z)], 
                        point_idx_offset=len(outer_points),
                        radial_subdiv=options.radial_subdiv,
                        ring_point_indices=[
                            None,
                            outer_point_indices[:options.radial_subdiv],
                            outer_point_indices[-options.radial_subdiv:],
                            None,
                            ]
                        )

        points = outer_points + inner_points
        facets = outer_facets + inner_facets
        facet_holestarts = outer_facet_holestarts + inner_facet_holestarts
        facet_markers = outer_facet_markers + inner_facet_markers

        # set regional max. volume
        mesh_info.regions.resize(1)
        mesh_info.regions[0] = [0, 0,(first_z+last_z)/2, 0, 1e-4]

        build_kwargs["volume_constraints"] = True
    else:
        points, facets, facet_holestarts, facet_markers = \
                generate_surface_of_revolution(rz, closure=closure,
                        radial_subdiv=options.radial_subdiv)

    #print rz
    #for i, p in enumerate(points):
        #print i, p
    #for i, f in enumerate(facets):
        #print i, f

    mesh_info.set_points(points)
    mesh_info.set_facets_ex(facets, facet_holestarts, facet_markers)
    #mesh_info.save_poly("gun")
    #mesh_info.save_nodes("gun")

    mesh = build(mesh_info, verbose=True, **build_kwargs)
    print "%d elements" % len(mesh.elements)
    mesh.write_vtk("gun.vtk")




if __name__ == "__main__":
    main()
