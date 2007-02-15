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

    mesh_info = MeshInfo()
    rz = parse_superfish_format("gun.am", 0.3)

    closure = EXT_OPEN

    build_kwargs = {}

    generate_inverse = False
    if generate_inverse:
        # chop off points with zero radius
        while rz[0][0] == 0:
            rz.pop(0)
        while rz[-1][0] == 0:
            rz.pop(-1)

        # construct outer cylinder
        ((min_r, max_r), (min_z, max_z)) = bounding_box(rz)
        first_z = rz[0][1]
        last_z = rz[-1][1]
        rz.extend([
                (max_r+2, last_z),
                (max_r+2, first_z),
                ])
        closure=EXT_CLOSED_IN_RZ

    generate_fine_tube = True
    if generate_fine_tube:
        # chop off points with zero radius
        while rz[0][0] == 0:
            rz.pop(0)
        while rz[-1][0] == 0:
            rz.pop(-1)

        # construct outer cylinder
        ((min_r, max_r), (min_z, max_z)) = bounding_box(rz)
        first_z = rz[0][1]
        last_z = rz[-1][1]
        first_r = rz[0][0]
        tube_r = first_r*0.6

        rz.insert(0, (tube_r, min_z))
        rz.append((tube_r, max_z))

        points, facets = generate_surface_of_revolution(rz,
                closure=closure)

        tube_points, tube_facets = generate_surface_of_revolution(
                [(0,min_z), (tube_r, min_z), (tube_r, max_z),
                    (0, max_z)], point_idx_offset=len(points))

        points.extend(tube_points)
        facets.extend(tube_facets)

        # set regional max. volume
        mesh_info.regions.resize(1)
        mesh_info.regions[0] = [0, 0,(max_z+min_z)/2, 0, 1e-4]

        build_kwargs["varvolume"] = True
    else:
        points, facets = generate_surface_of_revolution(rz,
                closure=closure)


    mesh_info.set_points(points)
    mesh_info.set_facets(facets, [0 for i in range(len(facets))])
    mesh_info.save_poly("gun")
    mesh_info.save_nodes("gun")
    mesh = build(mesh_info, **build_kwargs)
    print "%d elements" % len(mesh.elements)
    mesh.write_vtk("gun.vtk")




if __name__ == "__main__":
    main()
