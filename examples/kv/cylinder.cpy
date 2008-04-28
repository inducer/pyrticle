import hedge.mesh as user_mesh
mesh = user_mesh.make_cylinder_mesh(radius=tube_radius, height=tube_length, 
    periodic=True,
    max_volume=1000*units.MM**3, 
    radial_subdivisions=radial_subdiv)

