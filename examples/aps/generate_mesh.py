def main():
    from optparse import OptionParser

    description = "generate a vtk mesh from a superfish input file"
    parser = OptionParser(description=description)
    parser.add_option(
	    "--point-dist", dest="point_dist", default=0.3, type="float",
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

    from superfish import parse_superfish_format
    rz = parse_superfish_format(args[0], options.point_dist)

    build_kwargs = {}

    if options.generate_inverse:
        from rzmesh import inverse_mesh_info_from_rz
        mesh_info = inverse_mesh_info_from_rz(rz, options.radial_subdiv)
    elif options.tube_radius:
        from rzmesh import make_mesh_info_with_inner_tube
        mesh_info = make_mesh_info_with_inner_tube(rz, 
                options.tube_radius, options.radial_subdiv)

        build_kwargs["volume_constraints"] = True
    else:
        from rzmesh import make_mesh_info
        mesh_info = make_mesh_info_with_inner_tube(rz, options.radial_subdiv)

    #mesh_info.save_poly("gun")
    #mesh_info.save_nodes("gun")

    from meshpy.tet import build
    mesh = build(mesh_info, verbose=True, **build_kwargs)
    print "%d elements" % len(mesh.elements)
    mesh.write_vtk("gun.vtk")




if __name__ == "__main__":
    main()
