from tubemesh import make_extrusion_with_fine_core
full_mesh = make_extrusion_with_fine_core(
        rz=[
            (1*setup.tube_radius,0),
            (1*setup.tube_radius,setup.tube_length*0.333),
            (2*setup.tube_radius,setup.tube_length*0.333),
            (2*setup.tube_radius,setup.tube_length*0.666),
            (1*setup.tube_radius,setup.tube_length*0.666),
            (1*setup.tube_radius,setup.tube_length),
            ],
        inner_r=setup.tube_radius_inner, 
        max_volume_inner=setup.max_volume_inner,
        max_volume_outer=setup.max_volume_outer,
        radial_subdiv=setup.radial_subdiv,
        )

