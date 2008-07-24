from __future__ import division




def generate_arc_points(radius, start_angle, end_angle, 
        subdiv_degrees, include_final=True):
    from math import ceil, sin, cos, pi
    point_count = int(ceil((end_angle-start_angle)/subdiv_degrees))
    dphi = (end_angle-start_angle)/point_count

    if include_final:
        end_step = point_count + 1
    else:
        end_step = point_count

    for i in range(end_step):
        phi = pi/180*(start_angle + i*dphi)
        yield cos(phi)*radius, sin(phi)*radius

def make_magnetron_outline(
        cavity_count,
        cavity_angle,
        radius_anode,
        radius_cathode,
        radius_outer,
        anode_marker,
        cathode_marker,
        subdiv_degrees,
        cathode_cavities={},
        ):

    from math import sin, cos

    def round_trip_connect(start, end):
        for i in range(start, end):
            yield i, i+1
        yield end, start

    cathode_points = []
    anode_points = []
    cathode_facet_markers = []

    angle_step = 360/cavity_count

    for cav_idx in range(cavity_count):
        start_angle = angle_step * cav_idx
        if cav_idx in cathode_cavities:
            cathode_cavities[cav_idx](
                    start_angle,
                    angle_step,
                    cathode_points,
                    cathode_facet_markers)
        else:
            cathode_point_count_before = len(cathode_points)
            cathode_points.extend(generate_arc_points(
                radius_outer, start_angle, start_angle+cavity_angle,
                subdiv_degrees))
            cathode_points.extend(generate_arc_points(
                radius_cathode, 
                start_angle+cavity_angle, 
                start_angle+angle_step,
                subdiv_degrees))
            cathode_facet_markers.extend([cathode_marker] * 
                    (len(cathode_points) - cathode_point_count_before))


    anode_points = list(
            generate_arc_points(radius_anode, 0, 360, subdiv_degrees,
                include_final=False))

    return (anode_points+cathode_points, 
            list(round_trip_connect(0, len(anode_points)-1))
            + list(round_trip_connect(
                len(anode_points), 
                len(anode_points)+len(cathode_points)-1)),
            len(anode_points)*[anode_marker] + cathode_facet_markers)




def make_poly_mesh(points, facets, facet_markers, refinement_func=None):
    import meshpy.triangle as triangle
    mesh_info = triangle.MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets(facets, facet_markers)
    mesh_info.holes.resize(1)
    mesh_info.holes[0] = [0,0]
    return triangle.build(mesh_info, 
            refinement_func=refinement_func,
            generate_edges=True)




def main():

    def refinement_func(tri, area):
        return area > 2e-6

    cavity_angle = 20
    radius_anode = 0.0158
    radius_cathode = 0.0211
    radius_outer = 0.0411

    anode_marker = 1
    cathode_marker = 2
    open_marker = 3

    horn_radius = 0.08

    subdiv_degrees = 5

    def make_horn(start_angle, angle_step,
            cathode_points, cathode_facet_markers):
        from math import sin, pi
        cathode_points.extend([
            (horn_radius,0),
            (horn_radius,sin(pi/180*cavity_angle)*horn_radius),
            ])
        cathode_facet_markers.extend([open_marker, cathode_marker])
        cathode_point_count_before = len(cathode_points)
        cathode_points.extend(generate_arc_points(
            radius_cathode, 
            start_angle+cavity_angle, 
            start_angle+angle_step,
            subdiv_degrees))
        cathode_facet_markers.extend([cathode_marker] * 
                (len(cathode_points) - cathode_point_count_before))

    points, facets, facet_markers = make_magnetron_outline(
        cavity_count=6,
        cavity_angle=cavity_angle,
        radius_anode=radius_anode,
        radius_cathode=radius_cathode,
        radius_outer=radius_outer,
        cathode_marker=cathode_marker,
        anode_marker=anode_marker,
        cathode_cavities={0:make_horn},
        subdiv_degrees=subdiv_degrees,
        )

    mesh = make_poly_mesh(
            points,
            facets,
            facet_markers,
            refinement_func)

    import meshpy.triangle as triangle
    triangle.write_gnuplot_mesh("magnetron.dat", mesh)
    mesh.write_neu(open("magnetron.neu", "w"),
            bc={
                anode_marker:("anode", anode_marker), 
                cathode_marker:("cathode", cathode_marker),
                open_marker:("open", open_marker)})




if __name__ == "__main__":
    main()
