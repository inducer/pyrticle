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
    from kv import KVZIntervalBeam
    from random import seed
    seed(0)

    from pyrticle.units import SI
    units = SI()

    # discretization setup ----------------------------------------------------
    tube_length = 100*units.MM
    mesh = make_cylinder_mesh(radius=25*units.MM, height=tube_length, periodic=True)
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
    nparticles = 100000

    from pyrticle.cloud import ParticleCloud
    from pyrticle.reconstruction import ShapeFunctionReconstructor
    from pyrticle.pusher import MonomialParticlePusher
    cloud = ParticleCloud(discr, units, 
            ShapeFunctionReconstructor(),
            MonomialParticlePusher(),
            3, 3, verbose_vis=False)

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

    # timestep setup ----------------------------------------------------------
    vel = cloud.raw_velocities()
    def rhs(t, y):
        return ArithmeticList([
            vel, 
            0*vel, 
            ])

    stepper = RK4TimeStepper()
    t = 0

    # diagnostics setup -------------------------------------------------------
    from pytools.log import LogManager, \
            add_simulation_quantities, \
            add_general_quantities, \
            add_run_info, ETA
    from pyrticle.log import add_particle_quantities, add_beam_quantities, \
            ParticleCurrent
    logmgr = LogManager("no-charge.dat")
    add_run_info(logmgr)
    add_general_quantities(logmgr)
    add_simulation_quantities(logmgr, dt)
    add_particle_quantities(logmgr, cloud)
    add_beam_quantities(logmgr, cloud, axis=0, beam_axis=2)
    logmgr.add_quantity(ParticleCurrent(cloud, (1,0), tube_length))

    stepper.add_instrumentation(logmgr)
    cloud.add_instrumentation(logmgr)

    logmgr.set_constant("beta", beta)
    logmgr.set_constant("gamma", gamma)
    logmgr.set_constant("Q0", cloud_charge)
    logmgr.set_constant("n_part_0", nparticles)
    logmgr.set_constant("pmass", electrons_per_particle*units.EL_MASS)

    from kv import KVPredictedRadius
    logmgr.add_quantity(KVPredictedRadius(dt, 
        beam_v=beta*units.VACUUM_LIGHT_SPEED,
        predictor=beam.get_rms_predictor(axis=0),
        suffix="x_rms"))
    logmgr.add_quantity(KVPredictedRadius(dt, 
        beam_v=beta*units.VACUUM_LIGHT_SPEED,
        predictor=beam.get_total_predictor(axis=0),
        suffix="x_total"))

    logmgr.add_quantity(ETA(nsteps))

    logmgr.add_watches(["step", "t_eta", "(rx_rms-rx_rms_theory)/rx_rms_theory"])

    # timestep loop -----------------------------------------------------------
    for step in xrange(nsteps):
        logmgr.tick()

        cloud = stepper(cloud, t, dt, rhs)
        cloud.upkeep()
        t += dt

    vis.close()

    logmgr.tick()
    logmgr.save()

    _, _, err_table = logmgr.get_expr_dataset("(rx_rms-rx_rms_theory)/rx_rms_theory")
    print "Relative error (rms): %g" % max(err for step, err in err_table)
             



if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()

