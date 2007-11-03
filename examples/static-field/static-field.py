from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pylinear.operator as op
import pyrticle.units as units
import cProfile as profile




class StaticFieldSetup:
    def _check_velocities(self, deriv_dt=1e-12, final_time=1, check_steps=17,
            threshold=1e-4):
        """Check whether computed derivatives and velocities fit together."""
        from math import pi

        for step in range(check_steps):
            t = final_time/check_steps*step
            x1s = self.positions(t-deriv_dt)
            x2s = self.positions(t+deriv_dt)
            vs = self.velocities(t)

            for x1, x2, v in zip(x1s, x2s, vs):
                vpos = (x2-x1)/(2*deriv_dt)
                #u = num.vstack((vpos, v, (vpos-v)/comp.norm_2(v)))
                #print (vpos - v)/comp.norm_2(v)
                assert comp.norm_2(vpos - v)/comp.norm_2(v) < threshold





class LarmorScrew(StaticFieldSetup):
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
            
        self._check_velocities(final_time = 2*pi/self.omega_b*1000)

    def nsteps(self):
        return 10**4

    def final_time(self):
        return 1/self.vperp

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




class EBParallel(StaticFieldSetup):
    """Implements Problem 12.6b) in the third edition of Jackson."""
    def __init__(self, ez, bz, mass, charge, radius, nparticles):
        self.ez = ez
        self.bz = bz
        self.mass = mass
        self.charge = charge

        if nparticles == 1:
            self.particle_offsets = [num.zeros((3,))]
        elif nparticles == 4:
            self.particle_offsets = [
                    num.array([-0.1, -0.1, 0]),
                    num.array([ 0.1, -0.1, 0]),
                    num.array([-0.1,  0.1, 0]),
                    num.array([ 0.1,  0.1, 0])
                    ]
        else:
            raise ValueError, "invalid number of particles"

        self.shift = num.zeros((3,))

        c = self.c = units.VACUUM_LIGHT_SPEED
        self.R = mass*c/(charge*bz) # converted to SI
        self.rho = ez/(c*bz) # converted to SI

        self.scaling = radius/self.R

        from pytools import average
        self.shift[2] = average(self.positions(0))[2]

        self._check_t_phi()
        self._check_velocities(deriv_dt=1e-14, threshold=1e-8, final_time=self.R/self.c)

    def _check_t_phi(self):
        dt = 1/10

        def approx_deriv(f, x, deriv_dt=1e-8):
            if abs(x) < deriv_dt:
                return (f(x+deriv_dt)-f(x-deriv_dt))/(2*deriv_dt)
            else:
                step = abs(x)*deriv_dt
                return (f(x+step)-f(x-step))/(2*step)

        for i in range(0, 10):
            # check that t and phi are inverses of each other
            t = dt*i
            phi = self.phi(t)
            t2 = self.t(phi)

            if abs(t):
                assert abs(t-t2)/t < 4e-15
            else:
                assert abs(t-t2) < 4e-15

            # check that tprime works
            tprime_approx = approx_deriv(self.t, phi)
            tprime_real = self.tprime(phi)
            tprime_relerr = abs(tprime_approx-tprime_real)/tprime_real

            assert tprime_relerr < 1e-7

            # check that we can compute phiprime
            phiprime_approx = approx_deriv(self.phi, t)
            phiprime_real = 1/self.tprime(phi)
            phiprime_relerr = abs(phiprime_approx-phiprime_real)/phiprime_real

            assert tprime_relerr < 1e-7

    def nsteps(self):
        return 10**4

    def final_time(self):
        return self.R/self.c

    def e(self):
        return num.array([0,0,self.ez])

    def h(self):
        return 1/units.MU0*num.array([0, 0, self.bz])

    def fields(self, discr):
        class ZField:
            shape = (3,)
            
            def __init__(self, zval):
                self.zval = zval

            def __call__(self, x):
                return num.array([0,0,self.zval])

        from hedge.data import GivenFunction
        e = GivenFunction(ZField(self.ez)).volume_interpolant(discr)
        h = 1/units.MU0*GivenFunction(ZField(self.bz)).volume_interpolant(discr)

        return e, h

    def t(self, phi):
        from math import sqrt, sinh
        return self.R/self.rho*sqrt(1+self.scaling**2)*sinh(self.rho*phi)/self.c

    def tprime(self, phi):
        from math import sqrt, cosh
        return self.R*sqrt(1+self.scaling**2)*cosh(self.rho*phi)/self.c

    def phi(self, t):
        from pyrticle._internal import asinh
        from math import sqrt

        return asinh(self.rho*t*self.c/(self.R*sqrt(1+self.scaling**2)))/self.rho

    def positions(self, t):
        from math import sin, cos, pi, sqrt, cosh

        phi = self.phi(t)

        zconst = self.R/self.rho*sqrt(1+self.scaling**2)
        pos = (num.array([
                    self.scaling*self.R*sin(phi),
                    self.scaling*self.R*cos(phi),
                    zconst*cosh(self.rho*phi)
                    ])
                -self.shift)
        return [x+pos for x in self.particle_offsets]

    def velocities(self, t):
        from math import sin, cos, pi, sqrt, sinh

        phi = self.phi(t)
        phi_prime = 1/self.tprime(phi)
        zconst = self.R/self.rho*sqrt(1+self.scaling**2)
        return len(self.particle_offsets)*[
                num.array([
                    self.scaling*self.R*cos(phi)*phi_prime,
                    self.scaling*self.R*(-sin(phi))*phi_prime,
                    zconst*sinh(self.rho*phi)*self.rho*phi_prime
                    ])
                ]




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

    #setup = LarmorScrew(c*0.8, c*0.1, bz=1e-3, 
            #charge=-units.EL_CHARGE, mass=units.EL_MASS, nparticles=4)
    setup = EBParallel(ez=1e+5, bz=1e-3,
            charge=units.EL_CHARGE, mass=units.EL_MASS,
            radius=0.5*radius,
            nparticles=1)

    e, h = setup.fields(discr)

    init_positions = setup.positions(0)
    init_velocities = setup.velocities(0)

    #print "x", init_positions
    print "v", init_velocities

    cloud.add_particles(
            positions=init_positions,
            velocities=init_velocities,
            charges=setup.charge, masses=units.EL_MASS)

    final_time = setup.final_time()
    nsteps = setup.nsteps()
    dt = final_time/nsteps

    print "#elements=%d, dt=%s, #steps=%d" % (
            len(discr.mesh.elements), dt, nsteps)

    # timestepping ------------------------------------------------------------
    def rhs(t, y):
        return cloud.rhs(t, e, h)

    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    def check_result():
        from hedge.tools import cross

        deriv_dt = 1e-12

        dim = discr.dimensions
        all_x = setup.positions(t)
        all_v = setup.velocities(t)
        all_sim_v = cloud.velocities
        all_a = [(v2-v1)/(2*deriv_dt)
                for v1, v2 in zip(setup.velocities(t-deriv_dt), setup.velocities(t+deriv_dt))]
        all_sim_a = cloud.vis_info["lorentz_acc"] + cloud.vis_info["el_acc"]

        e = setup.e()
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

            real_a = num.array(cross(v, charge_over_mass*units.MU0*h)) + \
                    charge_over_mass*e

            #print "pos%d:" % i, comp.norm_2(x-sim_x)

            #print "vel%d:" % i, comp.norm_2(v-sim_v)
            #print "vel%d:..." % i, v, sim_v

            #print "acc%d:" % i, comp.norm_2(a-sim_a)
            #u = num.vstack((v, a, sim_a, real_a))
            #print "acc%d:\n%s" % (i, u)
            #raw_input()

            x_err = max(x_err, comp.norm_2(v-sim_v)/comp.norm_2(v))
            v_err = max(v_err, comp.norm_2(v-sim_v)/comp.norm_2(v))
            a_err = max(a_err, comp.norm_2(a-sim_a)/comp.norm_2(a))

        return x_err, v_err, a_err

    # make sure verbose-vis fields are filled
    rhs(t, cloud)

    errors = (0, 0, 0)

    for step in xrange(nsteps):
        if step % 1000 == 0:
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

    print
    print "l_inf errors (pos,vel,acc):", errors

    vis.close()




if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
