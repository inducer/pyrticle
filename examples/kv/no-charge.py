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
    cloud = ParticleCloud(discr, units, 3, 3, verbose_vis=False)

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

    max_logger = MaxBeamRadiusLogger(cloud.dimensions_pos)
    rms_logger = RMSBeamRadiusLogger(cloud.dimensions_pos, 0)

    for step in xrange(nsteps):
        if step % 100 == 0:
            print "timestep %d" % step

        vel = cloud.velocities()
        max_logger.update(t, cloud.positions, vel)
        rms_logger.update(t, cloud.positions, vel)

        cloud = stepper(cloud, t, dt, rhs)
        cloud.upkeep()
        t += dt

    vis.close()

    theory_no_charge_max = ChargelessKVRadiusPredictor(
            beam.radii[0], beam.emittances[0])
    theory_with_charge_max = KVRadiusPredictor(
            beam.radii[0], beam.emittances[0],
            xi=6e-5)

    theory_no_charge_rms = ChargelessKVRadiusPredictor(
            beam.rms_radii[0], beam.rms_emittances[0])
    theory_with_charge_max = KVRadiusPredictor(
            beam.rms_radii[0], beam.rms_emittances[0],
            xi=6e-5)

    rms_logger.generate_plot(title="Kapchinskij-Vladimirskij Beam Evolution",
            sim_label="RMS, simulated, no space charge", 
            outfile="beam-rad-rms.eps",
            theories=[
                ("RMS, theoretical, no space charge", theory_no_charge_rms), 
                #("RMS, theoretical, with space charge", theory_with_charge_rms)
                ])
    max_logger.generate_plot(title="Kapchinskij-Vladimirskij Beam Evolution",
            sim_label="100%, simulated, no space charge", 
            outfile="beam-rad-max.eps",
            theories=[
                ("100%, theoretical, no space charge", theory_no_charge_max), 
                ("100%, theoretical, with space charge", theory_with_charge_max)
                ])
    
    print "Relative error (max): %g" % max_logger.relative_error(theory_no_charge_max)
    print "Relative error (rms): %g" % rms_logger.relative_error(theory_no_charge_rms)
    max_logger.write_data("beam-radius-max.dat")
    rms_logger.write_data("beam-radius-rms.dat")
             



if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()

