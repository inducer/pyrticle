from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pylinear.operator as op
import pyrticle.units as units
import cProfile as profile




class LarmorScrew:
    """Implements Section 12.2 in the third edition of Jackson."""
    def __init__(self, vpar, vperp, bz, mass, charge, nparticles):
        self.vpar = vpar
        self.vperp = vperp
        self.bz = bz
        self.nparticles = nparticles
        self.mass = mass
        self.charge = charge

        c = units.VACUUM_LIGHT_SPEED
        gamma = (1-(vpar**2+vperp**2)/c**2)**(-1/2)
        self.omega_b = charge*bz/(gamma*mass) # (12.39) -> SI
        self.gyration_rad = vperp/abs(self.omega_b)

        from math import pi 
        print self.gyration_rad, self.omega_b, 2*pi/self.omega_b

        # self-check
        for v in self.velocities(0):
            mygamma = (1-comp.norm_2_squared(v)/c**2)**(-1/2)
            assert abs(mygamma-gamma)/gamma < 1e-13

    def e(self):
        return num.array([0,0,0])

    def h(self):
        return 1/units.MU0*num.array([0, 0, self.bz])

    def fields(self, discr):
        class ZField:
            shape = (3,)
            
            def __init__(self, zval):
                self.zval = zval

            def __call__(self, x):
                return num.array([0, 0, self.zval])

        from hedge.data import GivenFunction
        from pytools.arithmetic_container import ArithmeticList
        e = ArithmeticList([discr.volume_zeros() for i in range(3)])
        h = 1/units.MU0*GivenFunction(ZField(self.bz)).volume_interpolant(discr)

        return e, h

    def positions(self, t):
        from math import sin, cos, pi
        return [
            num.array([
                self.gyration_rad*( sin(angle + self.omega_b*t)),
                self.gyration_rad*( cos(angle + self.omega_b*t)),
                self.vpar*t])
            for angle in num.arange(0, 2*pi, 2*pi/self.nparticles)
            ]

    def velocities(self, t):
        from math import sin, cos, pi
        return [
            num.array([
                self.omega_b*self.gyration_rad*( cos(angle + self.omega_b*t)),
                self.omega_b*self.gyration_rad*(-sin(angle + self.omega_b*t)),
                self.vpar])
            for angle in num.arange(0, 2*pi, 2*pi/self.nparticles)
            ]




class EBParallel:
    """Implements Problem 12.6b) in the third edition of Jackson."""
    def __init__(self, ez, bz, mass, charge):
        self.ez = ez
        self.bz = bz
        #self.R = 


    def fields(self):
        class ZField:
            shape = (3,)
            
            def __init__(self, zval):
                self.zval = zval

            def __call__(self, x):
                return num.array([0,0,zval])

        from hedge.data import GivenFunction
        e = GivenFunction(ZField(ez)).volume_interpolant(discr)
        h = 1/units.MU0*GivenFunction(ZField(bz)).volume_interpolant(discr)

        return e, h

    def phi(self, x):
        pass

    def positions(self, t):
        pass

    def velocities(self, t):
        pass


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
    from pytools.arithmetic_container import \
            ArithmeticList, join_fields
    from pyrticle.cloud import ParticleCloud
    from math import sqrt, sin, cos, pi


    # discretization setup ----------------------------------------------------
    radius = 1*units.M
    full_mesh = make_cylinder_mesh(radius=radius, height=2*radius, periodic=True,
            radial_subdivisions=30)

    from hedge.parallel import guess_parallelization_context

    pcon = guess_parallelization_context()

    if pcon.is_head_rank:
        mesh = pcon.distribute_mesh(full_mesh)
    else:
        mesh = pcon.receive_mesh()

    discr = pcon.make_discretization(mesh, TetrahedralElement(1))
    vis = SiloVisualizer(discr)
    #vis = VtkVisualizer(discr, "pic")


    # particles setup ---------------------------------------------------------
    cloud = ParticleCloud(discr, 
            epsilon=units.EPSILON0, 
            mu=units.MU0, 
            verbose_vis=True)

    c = units.VACUUM_LIGHT_SPEED

    setup = LarmorScrew(c*0.8, c*0.1, bz=1e-3, 
            charge=-units.EL_CHARGE, mass=units.EL_MASS, nparticles=4)

    e, h = setup.fields(discr)

    init_velocities=setup.velocities(0)

    cloud.add_particles(
            positions=setup.positions(0),
            velocities=init_velocities,
            charges=-units.EL_CHARGE, masses=units.EL_MASS)

    vz = init_velocities[0][2]
    final_time = 2*radius/vz
    nsteps = 10**4
    dt = final_time/nsteps

    print "#elements=%d, dt=%s, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    deriv_dt = 1e-12

    def check_setup():
        """Check whether computed derivatives and velocities fit together."""
        dim = discr.dimensions
        check_steps = 17
        for step in range(check_steps):
            t = final_time/check_steps*step
            x1s = setup.positions(t)
            x2s = setup.positions(t+deriv_dt)
            vs = setup.velocities(t)

            for x1, x2, v in zip(x1s, x2s, vs):
                assert comp.norm_2((x2-x1)/deriv_dt - v)/comp.norm_2(v) < 1e-5

    check_setup()

    # timestepping ------------------------------------------------------------
    def rhs(t, y):
        return cloud.rhs(t, e, h)

    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    def check_result():
        from hedge.tools import cross

        dim = discr.dimensions
        all_x = setup.positions(t)
        all_v = setup.velocities(t)
        all_sim_v = cloud.velocities
        all_a = [(v2-v1)/deriv_dt 
                for v1, v2 in zip(setup.velocities(t), setup.velocities(t+deriv_dt))]
        all_sim_a = cloud.vis_info["lorentz_acc"]

        h = setup.h()

        x_err = 0
        v_err = 0
        a_err = 0

        for i in range(len(cloud)):
            x = all_x[i]
            sim_x = cloud.positions[i*dim:(i+1)*dim]
            v = all_v[i]
            sim_v = cloud.velocities[i*dim:(i+1)*dim]
            a = all_a[i]
            sim_a = all_sim_a[i*dim:(i+1)*dim]

            charge_over_mass = (setup.charge/setup.mass
                *sqrt(1-(comp.norm_2(v)/units.VACUUM_LIGHT_SPEED)**2))

            real_a = num.array(cross(v, charge_over_mass*units.MU0*h))

            #print "pos%d:" % i, comp.norm_2(x-sim_x)

            #print "vel%d:" % i, comp.norm_2(v-sim_v)
            #print "vel%d:..." % i, v, sim_v

            #print "acc%d:" % i, comp.norm_2(a-sim_a)
            #u = num.vstack((v, a, sim_a, real_a))
            #print "acc%d:\n%s" % (i, u)

            x_err = max(x_err, comp.norm_2(v-sim_v)/comp.norm_2(v))
            v_err = max(v_err, comp.norm_2(v-sim_v)/comp.norm_2(v))
            a_err = max(a_err, comp.norm_2(a-sim_a)/comp.norm_2(a))

        return x_err, v_err, a_err

    # make sure verbose-vis fields are filled
    rhs(t, cloud)

    errors = (0, 0, 0)

    for step in xrange(nsteps):
        if step % 100 == 0:
            print "timestep %d, t=%g secs=%f particles=%d" % (
                    step, t, time()-last_tstep, len(cloud))

            errors = tuple(
                    max(old_err, new_err) 
                    for old_err, new_err in zip(errors, check_result()))

            last_tstep = time()
            if True:
                visf = vis.make_file("pic-%04d" % step)

                mesh_scalars, mesh_vectors = \
                        cloud.add_to_vis(vis, visf, time=t, step=step)

                if False:
                    vis.add_data(visf,
                            scalars=mesh_scalars,
                            vectors=[
                                ("e", e), 
                                ("h", h), 
                                ]
                            + mesh_vectors,
                            write_coarse_mesh=True,

                            time=t, step=step)
                else:
                    vis.add_data(visf, write_coarse_mesh=True, time=t, step=step)
                visf.close()

        cloud = stepper(cloud, t, dt, rhs)
        cloud.upkeep()

        t += dt

    print "l^inf errors (pos,vel,acc):", errors

    vis.close()




if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
