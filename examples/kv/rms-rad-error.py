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
    
    pcounts = range(1000, 40000, 100)
    #pcounts = [1000, 5000, 10000, 20000, 30000, 40000, 50000, 100000] 
    rms_radii = []
    for nparticles in pcounts:
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

        from pyrticle.beam import rms_beam_size
        rms_radii.append(rms_beam_size(cloud, 0))

    from pylab import plot, show
    plot(pcounts, rms_radii)
    show()



             



if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()


