import special_meshes as user_meshes
mesh = user_meshes.make_fine_center_rect_mesh(
        a=(-0.5, -tube_width/2),
        b=(-0.5+tube_length, tube_width/2),
        periodicity=(tube_periodic, False),
        subdivisions=(10,5), inner_subdivisions=(20,3),
        inner_width=0.15,
        max_area=tube_max_tri_area, inner_max_area=0.15*tube_max_tri_area)

pusher = PushAverage()
reconstructor = RecShape()
shape_exponent = 2
nparticles = 10000
shape_bandwidth = "optimize"
vis_interval = 10

element_order = 5

user_c0 = units.VACUUM_LIGHT_SPEED()

sigma_x = 0.035 * num.ones((2,))
sigma_v = num.array([user_c0*0.9*1e-2, user_c0*0.9*1e-2])
