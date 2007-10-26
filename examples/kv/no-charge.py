from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pyrticle.units as units
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
            add_kv_xy_particles, \
            KVRadiusPredictor, \
            BeamRadiusLogger
    from random import seed
    seed(0)

    # discretization setup ----------------------------------------------------
    mesh = make_cylinder_mesh(radius=25*units.MM, height=100*units.MM, periodic=True)
    #mesh = make_box_mesh([1,1,2], max_volume=0.01)

    discr = Discretization(mesh, TetrahedralElement(3))
    vis = SiloVisualizer(discr)
    #vis = VtkVisualizer(discr, "pic")

    dt = discr.dt_factor(units.C0) / 2
    final_time = 1*units.M/units.C0
    nsteps = int(final_time/dt)+1
    dt = final_time/nsteps

    print "#elements=%d, dt=%s, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    # particles setup ---------------------------------------------------------
    nparticles = 10000

    cloud = ParticleCloud(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            verbose_vis=False)

    cloud_charge = 1e-9 * units.C
    particle_charge = cloud_charge/nparticles
    electrons_per_particle = cloud_charge/nparticles/units.EL_CHARGE
    print "e-/particle = ", electrons_per_particle 

    emittance = 5 * units.MM * units.MRAD
    initial_radius = 2.5*units.MM

    el_energy = 5.2e6 * units.EV
    #el_energy = units.EL_REST_ENERGY*1.00001
    el_lorentz_gamma = el_energy/units.EL_REST_ENERGY
    #el_lorentz_gamma = 100000
    beta = (1-1/el_lorentz_gamma**2)**0.5
    print "v = %g%% c" % (beta*100)

    add_kv_xy_particles(nparticles, cloud, discr, 
            charge=0, 
            mass=electrons_per_particle*units.EL_MASS,
            radii=[2.5*units.MM, 2.5*units.MM],
            z_length=5*units.MM,
            z_pos=10*units.MM,
            emittances=[emittance, emittance], 
            beta=beta)

    # timestepping ------------------------------------------------------------
    def rhs(t, y):
        return ArithmeticList([
            cloud.velocities, 
            0*cloud.velocities, 
            ])

    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    r_logger = BeamRadiusLogger(cloud.mesh_info.dimensions,
            initial_radius, emittance)

    for step in xrange(nsteps):
        r_logger.update(t, cloud.positions, cloud.velocities)
        cloud = stepper(cloud, t, dt, rhs)
        cloud.upkeep()
        t += dt

    vis.close()

    theory = KVRadiusPredictor(initial_radius, emittance)
    r_logger.generate_plot("Kapchinskij-Vladimirskij Beam Evolution, "
            "no space charge", theory=theory)
    
    print "Relative error: %g" % r_logger.relative_error(theory)
             




if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()

