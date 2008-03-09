from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pylinear.operator as op
import cProfile as profile
import pytools




def main():
    from hedge.element import TriangularElement
    from hedge.timestep import RK4TimeStepper
    from hedge.discretization import \
            Discretization, \
            pair_with_boundary
    from hedge.visualization import VtkVisualizer, SiloVisualizer
    from hedge.tools import dot
    from math import sqrt, pi
    from pytools.arithmetic_container import join_fields

    from random import seed
    seed(0)

    from pyrticle.units import SI
    units = SI()

    # user interface ----------------------------------------------------------
    def make_setup():
        c0 = units.VACUUM_LIGHT_SPEED

        tube_width = 1
        variables = {
                "mesh": None,
                "tube_length": 2,
                "tube_width": tube_width,
                "tube_periodic": True,
                "tube_max_tri_area": 0.02,

                "element_order": 7,
                "shape_exponent": 2,

                "pusher": None,
                "reconstructor": None,

                "mean_x": num.array([0,0]),
                "sigma_x": num.array([0.05,0.5*tube_width]),
                "nparticles": 1000,
                "cloud_charge": -1e-9 * units.C,
                }
        
        from pyrticle.reconstruction import \
                ShapeFunctionReconstructor, \
                NormalizedShapeFunctionReconstructor, \
                AdvectiveReconstructor
        from pyrticle.pusher import \
                MonomialParticlePusher, \
                AverageParticlePusher
        from hedge.mesh import make_rect_mesh

        constants = {
                "num": num,
                "comp": comp,
                "units": units,

                "RecShape": ShapeFunctionReconstructor,
                "RecNormShape": NormalizedShapeFunctionReconstructor,
                "RecAdv": AdvectiveReconstructor,

                "PushMonomial": MonomialParticlePusher,
                "PushAverage": AverageParticlePusher,

                "make_rect_mesh": make_rect_mesh,
                }

        doc = {
                "chi": "relative speed of hyp. cleaning (None for no cleaning)",
                "shape_bandwidth": "either 'optimize', 'guess' or a positive real number",
                }

        from pytools import gather_parameters_from_user
        return gather_parameters_from_user(variables, constants, doc)

    setup = make_setup()

    # discretization setup ----------------------------------------------------
    if setup.mesh is None:
        from hedge.mesh import make_rect_mesh
        setup.mesh = make_rect_mesh(
                a=(-0.5, -setup.tube_width/2),
                b=(-0.5+setup.tube_length, setup.tube_width/2),
                periodicity=(setup.tube_periodic, False),
                subdivisions=(10,5),
                max_area=setup.tube_max_tri_area)
    elif setup.mesh == "glued":
        from special_meshes import make_glued_rect_mesh
        setup.mesh = make_glued_rect_mesh(
                a=(-0.5, -setup.tube_width/2),
                b=(-0.5+setup.tube_length, setup.tube_width/2),
                periodicity=(setup.tube_periodic, False),
                subdivisions=(10,5),
                max_area=setup.tube_max_tri_area)
    else:
        from hedge.mesh import Mesh
        assert isinstance(setup.mesh, Mesh)

    from hedge.parallel import guess_parallelization_context

    pcon = guess_parallelization_context()

    if pcon.is_head_rank:
        mesh = pcon.distribute_mesh(setup.mesh)
    else:
        mesh = pcon.receive_mesh()

    discr = pcon.make_discretization(mesh, TriangularElement(setup.element_order))
    vis = SiloVisualizer(discr)
    #vis = VtkVisualizer(discr, "pic")

    from hedge.operators import TEMaxwellOperator, DivergenceOperator
    from hedge.mesh import TAG_ALL, TAG_NONE

    max_op = TEMaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)

    # particles setup ---------------------------------------------------------
    from pyrticle._internal import StatsGatherer
    sg = StatsGatherer()

    from pyrticle.reconstruction import Reconstructor
    from pyrticle.pusher import Pusher

    assert isinstance(setup.reconstructor, Reconstructor), \
            "must specify valid reconstructor"
    assert isinstance(setup.pusher, Pusher), \
            "must specify valid reconstructor"

    from pyrticle.cloud import ParticleCloud
    cloud = ParticleCloud(discr, units, setup.reconstructor, setup.pusher,
                dimensions_pos=2, dimensions_velocity=2,
                verbose_vis=True)

    from pyrticle.cloud import guess_shape_bandwidth
    guess_shape_bandwidth(cloud, setup.shape_exponent)

    integrals = []

    bbox = mesh.bounding_box
    real_tube_length = bbox[1][0] - bbox[0][0]

    for particle in xrange(setup.nparticles):
        cloud.clear_particles()

        from random import gauss, uniform
        pos = num.array([
            gauss(setup.mean_x[0], setup.sigma_x[0]) % real_tube_length + bbox[0][0],
            gauss(setup.mean_x[1],0.1)
            ])
        cloud.add_particles(
                positions=[pos],
                velocities=[num.array([0,0])],
                charges=setup.cloud_charge, 
                masses=units.EL_MASS)

        assert len(cloud) == 1
        rho = cloud.reconstruct_rho()

        from hedge.discretization import integral
        int_rho = integral(discr, rho)
        sg.add(int_rho)
        integrals.append(int_rho)

        visf = vis.make_file("particle-%04d" % particle)
        cloud.add_to_vis(vis, visf)
        vis.add_data(visf, [
                    ("rho", rho),
                    ],)
        visf.close()

    descr = "mean_x=%s, mean=%g, stddev=%g, min=%g, max=%g" % (setup.mean_x,
            sg.mean(), sg.standard_deviation(),
            sg.minimum(), sg.maximum())

    print descr 

    from pylab import hist, show,title
    title(descr)
    hist(integrals, 100)
    show()




if __name__ == "__main__":
    main()

