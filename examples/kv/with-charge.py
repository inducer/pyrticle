from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pyrticle.units as units
import cProfile as profile




class RTLogger:
    def __init__(self, dimensions, a0, eps):
        self.outf = open("particle-r-t.dat", "w")
        self.outf_theory = open("particle-r-t-theory.dat", "w")
        self.dimensions = dimensions

        self.a0 = a0
        self.eps = eps

    def __call__(self, t, positions, velocities):
        from math import sqrt,pi
        from pytools import argmax

        dim = self.dimensions
        nparticles = len(positions) // dim
        pn = argmax(positions[i*dim+0]**2 +positions[i*dim+1]**2
                for i in xrange(nparticles))
        r = comp.norm_2(positions[pn*dim+0:pn*dim+2])
        vz = velocities[pn*dim+2]
        s = t*vz
        self.outf.write("%g\t%g\n" % (s,r))
        self.outf.flush()

        a = sqrt(self.a0**2+(self.eps/self.a0)**2 * s**2)
        self.outf_theory.write("%g\t%g\n" % (s,a))
        self.outf_theory.flush()

            


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
    from hedge.tools import dot, cross
    from math import sqrt, pi
    from pytools.arithmetic_container import \
            ArithmeticList, concatenate_fields
    from hedge.operators import MaxwellOperator
    from pyrticle.cloud import ParticleCloud
    from pyrticle.distribution import add_kv_xy_particles
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

    def l2_norm(field):
        return sqrt(dot(field, discr.mass_operator*field))

    op = MaxwellOperator(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            upwind_alpha=1)

    # particles setup ---------------------------------------------------------
    nparticles = 1000

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

    full_charge_at_time = final_time*0.2

    # timestepping ------------------------------------------------------------
    fields = concatenate_fields(
            [discr.volume_zeros() for i in range(2*discr.dimensions)],
            [cloud])

    def rhs(t, y):
        e = y[:3]
        h = y[3:6]

        if True:
            maxwell_rhs = op.rhs(t, y[0:6])
            rho, j = cloud.reconstruct_densities()
            cloud_rhs = [cloud.rhs(t, e, h)]
        else:
            maxwell_rhs = ArithmeticList(6*[discr.volume_zeros()])
            rho = discr.volume_zeros()
            j = ArithmeticList(3*[discr.volume_zeros()])
            cloud_rhs = [ArithmeticList([
                cloud.velocities, 
                0*cloud.velocities, 
                ])]

        rhs_e = maxwell_rhs[:3]
        rhs_h = maxwell_rhs[3:6]
        return concatenate_fields(
                rhs_e + 1/units.EPSILON0*j,
                rhs_h,
                cloud_rhs,
                )

    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    rt_logger = RTLogger(cloud.mesh_info.dimensions,
            initial_radius, emittance)

    for step in xrange(nsteps):
        if True:
            visf = vis.make_file("pic-%04d" % step)

            mesh_scalars, mesh_vectors = \
                    cloud.add_to_vis(vis, visf, time=t, step=step)
            vis.add_data(visf,
                    scalars=mesh_scalars,
                    vectors=[("e", fields[0:3]), 
                        ("h", fields[3:6]), ]
                    + mesh_vectors
                    ,
                    write_coarse_mesh=True,
                    time=t, step=step)
            visf.close()

        rt_logger(t, cloud.positions, cloud.velocities)

        if False:
            myfields = [fields]
            fields = profile.runctx("myfields[0] = stepper(fields, t, dt, rhs)", 
                    globals(), locals(), "pic-%04d.prof" % step)
            fields = myfields[0]
        else:
            fields = stepper(fields, t, dt, rhs)

        cloud.upkeep()

        print "timestep %d, t=%g l2[e]=%g l2[h]=%g secs=%f particles=%d" % (
                step, t, l2_norm(fields[0:3]), l2_norm(fields[3:6]),
                time()-last_tstep, len(cloud))
        if False:
            print "searches: same=%d, normal=%d, vertex=%d, global=%d, periodic=%d" % (
                    cloud.same_searches.pop(),
                    cloud.normal_searches.pop(),
                    cloud.vertex_searches.pop(),
                    cloud.global_searches.pop(),
                    cloud.periodic_hits.pop(),
                    )
            print "shape-adds: neighbor=%d vertex=%d" % (
                    cloud.neighbor_shape_adds.pop(),
                    cloud.vertex_shape_adds.pop(),
                    )

        last_tstep = time()

        t += dt

        if t < full_charge_at_time:
            charge_now = t/full_charge_at_time*particle_charge
            cloud.charges = charge_now * \
                    num.ones((len(cloud.containing_elements),))

    vis.close()
             




if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()