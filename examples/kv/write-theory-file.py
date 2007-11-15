
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

    import sys
    if len(sys.argv) != 3:
        print "usage: %s outfile infile" % sys.argv[0]
        return

    outfile = sys.argv[1]
    infile = sys.argv[2]

    seed(0)

    beam_radius = 2.5*units.MM

    nparticles = 3000
    cloud_charge = 1e-9 * units.C
    electrons_per_particle = cloud_charge/nparticles/units.EL_CHARGE
    print "e-/particle = ", electrons_per_particle 

    emittance = 5 * units.MM * units.MRAD
    initial_radius = 2.5*units.MM

    el_energy = 5.2e6 * units.EV
    #el_energy = units.EL_REST_ENERGY*1.00001
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
            z_length=5*units.MM,
            z_pos=10*units.MM,
            beta=beta)

    rms_r_logger = RMSBeamRadiusLogger(3, 0)
    rms_r_logger.read_data(infile)

    rms_theory_no_charge = KVRadiusPredictor(
            beam.rms_radii[rms_r_logger.axis], 
            beam.rms_emittances[rms_r_logger.axis])
    rms_theory_with_charge = KVRadiusPredictor(
            beam.rms_radii[rms_r_logger.axis], 
            beam.rms_emittances[rms_r_logger.axis],
            xi=beam.get_space_charge_parameter())

    rms_r_logger.generate_plot(
            title="Kapchinskij-Vladimirskij Beam Evolution",
            sim_label="RMS, simulated, with space charge", 
            outfile=outfile,
            theories=[
                ("RMS, theoretical, with space charge", rms_theory_with_charge),
                ("RMS, theoretical, no space charge", rms_theory_no_charge)
                ],
            no_overwrite=True)





if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
