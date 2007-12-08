def _angle_distance(start, end):
    from math import pi
    if end < start:
        end += 2*pi
    return end-start


def parse_superfish_format(filename, max_point_dist=0.1):
    import re
    
    def chop_comment(line):
        com_start = line.find(";")
        if com_start != -1:
            return line[:com_start]
        else:
            return line

    lines = [chop_comment(line.strip().lower()) for line in file(filename, "r").readlines()]
    lines = [line for line in lines if not line.startswith("!")]
    lines = [line for line in lines if line]
    lines = [line for line in lines if line.startswith("&po")]

    # parse the file
    line_re = re.compile(r"^\&po\s+(.*)\s*\&$")
    key_val_re = re.compile(r"^([0-9a-zA-Z]+)\s*\=\s*([-+0-9.a-zA-Z]+)")

    proplines = []

    for line in lines:
        line_match = line_re.match(line)
        assert line_match
        line = line_match.group(1)

        properties = {}

        pos = 0 
        while pos < len(line):
            if line[pos] in [" ", ","]:
                pos += 1
                continue
            datum_match = key_val_re.match(line[pos:])
            assert datum_match
            properties[datum_match.group(1)] = float(datum_match.group(2))
            pos += datum_match.end()

        proplines.append(properties)

    #for p in proplines:
        #print p

    # concoct x-y-points
    from math import atan2, sin, cos, pi, ceil
    import pylinear.array as num
    import pylinear.computation as comp
    from pytools import argmin

    def add_point(pt):
        points.append((pt[1], pt[0]))

    points = []
    for p in proplines:
        shape_type = 1
        if "nt" in p:
            shape_type = p["nt"]

        if shape_type == 1:
            points.append(num.array((p["x"], p["y"])))
        elif shape_type == 2:
            # draw arc
            last_point = points[-1]
            if "x0" in p:
                center = num.array((p["x0"], p["y0"]))
            else:
                center = num.zeros((2,))

            if "r" in p or "radius" in p:
                if "r" in p:
                    r = p["r"]
                    assert "radius" not in p
                else:
                    r = p["radius"]
                    assert "r" not in p

                theta = pi/180.*p["theta"]

                tangent = num.array((cos(theta), sin(theta)))
                upward_normal = num.array((-sin(theta), cos(theta)))

                #print 180*theta/pi, last_point, center, tangent, upward_normal, last_point-center
                if (last_point - center)*upward_normal < 0:
                    # draw ccw (positive) arc / upward
                    phi_start = 3*pi/2 + theta
                    phi_end = phi_start + pi/2
                else:
                    # draw cw (negative) arc / downward
                    phi_start = pi/2 + theta
                    phi_end = phi_start - pi/2
            elif "x" in p:
                start_pt = last_point - center
                end_pt = num.array((p["x"], p["y"]))

                r = comp.norm_2(end_pt)
                r2 = comp.norm_2(start_pt)

                assert abs(r-r2)/r < 1e-2

                phi_start = atan2(start_pt[1], start_pt[0])
                raw_phi_end = atan2(end_pt[1], end_pt[0])

                phi_end_choices = [
                        raw_phi_end, 
                        raw_phi_end-2*pi, 
                        raw_phi_end+2*pi]
                phi_end = phi_end_choices[argmin(
                    abs(phi_end-phi_start) for phi_end in phi_end_choices)]
            else:
                raise ValueError, "arcs (nt=2) must have either r/radius or x/y specified"

            steps = int(ceil((abs(phi_end-phi_start))*r/max_point_dist))
            dphi = (phi_end-phi_start) / steps
            phi = phi_start+dphi

            for i in range(steps):
                points.append(center + r*num.array((cos(phi), sin(phi))))
                phi += dphi
        else:
            raise ValueError, "unhandled shape type: nt=%s" % shape_type

    # assert that last point coincides with first
    assert ((points[0][0]-points[-1][0])**2 + (points[0][1]-points[-1][1])**2) < 1e-10
    # and kill it
    points.pop()

    # now find the on-axis line, if any, and shift it so straddles the array boundary
    assert len(points) > 1

    for i, (z, r) in enumerate(points):
        next_i = (i+1)%len(points)
        next_r = points[next_i][1]

        if r == 0 and next_r == 0:
            if i < len(points):
                # if the on-axis line is already the last line in the polygon,
                # then that's fine--no need to change. Otherwise, shift so that
                # this becomes the case.

                from pytools import shift
                points = shift(points, len(points)-1-i)
                break

    # assert start and end r are 0
    assert points[0][1] == 0
    assert points[-1][1] == 0

    outf = file("gun-lines.dat", "w")
    for p in points:
        outf.write("%f\t%f\n" % (p[0], p[1]))
    outf.close()

    # turn (x,y) into (r,z)
    return [(p[1], p[0]) for p in points]


            


