from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pylinear.operator as op
import cProfile as profile




def main():
    from hedge.element import TetrahedralElement
    from hedge.timestep import RK4TimeStepper
    from hedge.mesh import \
            make_box_mesh, \
            make_cylinder_mesh
    from hedge.discretization import \
            Discretization, \
            pair_with_boundary
    from hedge.visualization import VtkVisualizer, SiloVisualizer
    from hedge.tools import dot
    from math import sqrt, pi
    from pytools.arithmetic_container import \
            ArithmeticList, join_fields
    from hedge.operators import MaxwellOperator, DivergenceOperator
    from pyrticle.cloud import ParticleCloud
    from kv import \
            KVZIntervalBeam, \
            KVRadiusPredictor, \
            MaxBeamRadiusLogger, \
            RMSBeamRadiusLogger
    from tubemesh import make_cylinder_with_fine_core
    from random import seed
    from pytools.stopwatch import Job

    from pyrticle.units import SI
    units = SI()

    seed(0)

    nparticles = 10000
    beam_radius = 2.5 * units.MM
    emittance = 5 * units.MM * units.MRAD
    final_time = 0.1*units.M/units.VACUUM_LIGHT_SPEED
    field_dump_interval = 1
    tube_length = 20*units.MM

    # discretization setup ----------------------------------------------------
    job = Job("mesh")

    full_mesh = make_cylinder_with_fine_core(
            r=10*beam_radius, inner_r=1*beam_radius, 
            min_z=0, max_z=tube_length,
            max_volume_inner=10*units.MM**3,
            max_volume_outer=100*units.MM**3,
            radial_subdiv=10,
            )
    job.done()

    from hedge.parallel import guess_parallelization_context

    pcon = guess_parallelization_context()

    if pcon.is_head_rank:
        mesh = pcon.distribute_mesh(full_mesh)
    else:
        mesh = pcon.receive_mesh()

    job = Job("discretization")
    discr = pcon.make_discretization(mesh, TetrahedralElement(3))
    job.done()

    vis = SiloVisualizer(discr)
    #vis = VtkVisualizer(discr, "pic")

    max_op = MaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)
    div_op = DivergenceOperator(discr)

    dt = discr.dt_factor(max_op.c) / 2
    nsteps = int(final_time/dt)+1
    dt = final_time/nsteps

    print "#elements=%d, dt=%s, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    def l2_norm(field):
        return sqrt(dot(field, discr.mass_operator*field))
    def l2_error(field, true):
        return l2_norm(field-true)/l2_norm(true)

    # particles setup ---------------------------------------------------------
    cloud = ParticleCloud(discr, units, 3, 3, verbose_vis=False)

    cloud_charge = -1e-9 * units.C
    electrons_per_particle = abs(cloud_charge/nparticles/units.EL_CHARGE)
    print "e-/particle = ", electrons_per_particle 

    #el_energy = 5.2e6 * units.EV
    el_energy = units.EL_REST_ENERGY
    el_lorentz_gamma = el_energy/units.EL_REST_ENERGY
    #el_lorentz_gamma = 100000
    beta = (1-1/el_lorentz_gamma**2)**0.5
    gamma = 1/sqrt(1-beta**2)
    print "beta = %g, gamma = %g" % (beta, gamma)

    beam = KVZIntervalBeam(units, nparticles, 
            p_charge=cloud_charge/nparticles, 
            p_mass=electrons_per_particle*units.EL_MASS,
            radii=2*[beam_radius],
            emittances=2*[5 * units.MM * units.MRAD], 
            z_length=tube_length,
            z_pos=tube_length/2,
            beta=beta)
    beam.add_to(cloud, discr)

    # intial condition --------------------------------------------------------
    from pyrticle.cloud import compute_initial_condition
    fields = compute_initial_condition(pcon, discr, cloud, 
            mean_beta=num.array([0, 0, beta]), max_op=max_op,
            debug=True)

    # check against theory ----------------------------------------------------
    from pytools.arithmetic_container import work_with_arithmetic_containers
    ac_multiply = work_with_arithmetic_containers(num.multiply)

    q_per_unit_z = cloud_charge/beam.z_length
    class TheoreticalEField():
        shape = (3,)

        def __call__(self, x):
            r = comp.norm_2(x[:2])
            if r >= max(beam.radii):
                xy_unit = x/r
                xy_unit[2] = 0
                return xy_unit*((gamma*q_per_unit_z)
                        /
                        (2*pi*r*max_op.epsilon))
            else:
                return num.zeros((3,))

    def theory_indicator(x):
        r = comp.norm_2(x[:2])
        if r >= max(beam.radii):
            return 1
        else:
            return 0

    e_theory = discr.interpolate_volume_function(TheoreticalEField())
    theory_ind = discr.interpolate_volume_function(theory_indicator)

    restricted_e = ac_multiply(fields.e, theory_ind)

    print "L2 error in outer E field: %g" % l2_error(restricted_e, e_theory)

    visf = vis.make_file("e_comparison")
    mesh_scalars, mesh_vectors = \
            cloud.add_to_vis(vis, visf)
    vis.add_data(visf, [
        ("e", restricted_e), 
        ("e_theory", e_theory), 
        ]
        + mesh_vectors
        + mesh_scalars
        )
    visf.close()





if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
