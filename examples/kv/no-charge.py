from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import cProfile as profile




def main():
    from hedge.element import TetrahedralElement
    from hedge.timestep import RK4TimeStepper
    from hedge.mesh import \
            make_box_mesh, \
            make_cylinder_mesh
    from hedge.discretization import Discretization
    from hedge.visualization import VtkVisualizer, SiloVisualizer
    from hedge.tools import dot, cross
    from pytools.arithmetic_container import ArithmeticList
    from pyrticle.cloud import ParticleCloud
    from kv import \
            KVZIntervalBeam, \
            ChargelessKVRadiusPredictor, \
            KVRadiusPredictor, \
            MaxBeamRadiusLogger, \
            RMSBeamRadiusLogger
    from random import seed
    seed(0)

    from pyrticle.units import SI
    units = SI()

    # discretization setup ----------------------------------------------------
    mesh = make_cylinder_mesh(radius=25*units.MM, height=100*units.MM, periodic=True)
    #mesh = make_box_mesh([1,1,2], max_volume=0.01)

    discr = Discretization(mesh, TetrahedralElement(3))
    vis = SiloVisualizer(discr)
    #vis = VtkVisualizer(discr, "pic")

    dt = discr.dt_factor(units.VACUUM_LIGHT_SPEED) / 2
    final_time = 1*units.M/units.VACUUM_LIGHT_SPEED
    nsteps = int(final_time/dt)+1
    dt = final_time/nsteps

    print "#elements=%d, dt=%s, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    # particles setup ---------------------------------------------------------
    nparticles = 1000

    from hedge.operators import MaxwellOperator
    max_op = MaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)
    cloud = ParticleCloud(max_op, 3, 3, verbose_vis=False)

    cloud_charge = 1e-9 * units.C
    particle_charge = cloud_charge/nparticles
    electrons_per_particle = cloud_charge/nparticles/units.EL_CHARGE
    print "e-/particle = ", electrons_per_particle 

    el_energy = 5.2e6 * units.EV
    #el_energy = units.EL_REST_ENERGY*1.00001
    gamma = el_energy/units.EL_REST_ENERGY
    #gamma = 100000
    beta = (1-1/gamma**2)**0.5
    print "v = %g%% c" % (beta*100)

    beam = KVZIntervalBeam(units, nparticles, 
            p_charge=0, 
            p_mass=electrons_per_particle*units.EL_MASS,
            radii=2*[2.5*units.MM],
            emittances=2*[5 * units.MM * units.MRAD], 
            z_length=5*units.MM,
            z_pos=10*units.MM,
            beta=beta)
    beam.add_to(cloud, discr)

    # timestepping ------------------------------------------------------------
    vel = cloud.velocities()
    def rhs(t, y):
        return ArithmeticList([
            vel, 
            0*vel, 
            ])

    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    r_logger = MaxBeamRadiusLogger(cloud.mesh_info.dimensions)

    for step in xrange(nsteps):
        if step % 100 == 0:
            print "timestep %d" % step
        r_logger.update(t, cloud.positions, cloud.velocities())
        cloud = stepper(cloud, t, dt, rhs)
        cloud.upkeep()
        t += dt

    vis.close()

    theory_no_charge = ChargelessKVRadiusPredictor(
            beam.radii[0], beam.emittances[0]
            )
    theory_with_charge = KVRadiusPredictor(
            beam.radii[0], beam.emittances[0],
            xi=6e-5)

    r_logger.generate_plot("Kapchinskij-Vladimirskij Beam Evolution, "
            "no space charge", 
            theories=[
                ("theoretical, no space charge", theory_no_charge), 
                ("theoretical, with space charge", theory_with_charge)
                ])
    
    print "Relative error: %g" % r_logger.relative_error(theory_no_charge)
             




if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()

