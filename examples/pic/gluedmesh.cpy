mesh = make_glued_rect_mesh(
        a=(-0.5, -tube_width/2),
        b=(-0.5+tube_length, tube_width/2),
        periodicity=(tube_periodic, False),
        subdivisions=(10,5),
        max_area=tube_max_tri_area)
