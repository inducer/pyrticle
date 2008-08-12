from __future__ import division
import numpy
import numpy.linalg as la




class StaticFieldSetup:
    def __init__(self, mass, charge, c):

        self.mass = mass
        self.charge = charge
        self.c = c

    def gamma(self, v):
        value = (1-numpy.dot(v, v)/self.c**2)**(-0.5)
        if value < 0:
            raise RuntimeError, "particle velocity > speed of light"
        return value

    def momenta(self, t):
        return [self.mass*self.gamma(v)*v for v in self.velocities(t)]

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
                #u = numpy.vstack((vpos, v, (vpos-v)/la.norm(v)))
                #print (vpos - v)/la.norm(v)
                assert la.norm(vpos - v)/la.norm(v) < threshold





class LarmorScrew(StaticFieldSetup):
    """Implements Section 12.2 in the third edition of Jackson."""
    def __init__(self, units, mass, charge, c, vpar, vperp, bz, nparticles):
        from math import pi 

        self.units = units

        StaticFieldSetup.__init__(self, mass, charge, c)
        self.vpar = vpar
        self.vperp = vperp
        self.bz = bz
        self.nparticles = nparticles

        gamma = (1-(vpar**2+vperp**2)/c**2)**(-1/2)
        self.omega_b = charge*bz/(gamma*mass) # (12.39) -> SI
        self.gyration_rad = vperp/abs(self.omega_b)

        #print self.gyration_rad, self.omega_b, 2*pi/self.omega_b

        # self-check
        for v in self.velocities(0):
            mygamma = (1-numpy.dot(v, v)/c**2)**(-1/2)
            assert abs(mygamma-gamma)/gamma < 1e-13
            
        self._check_velocities(final_time = 2*pi/self.omega_b*1000)

    def nsteps(self):
        return 10**4

    def final_time(self):
        return 1/self.vperp

    def e(self):
        return numpy.array([0,0,0])

    def h(self):
        return 1/self.units.MU0*numpy.array([0, 0, self.bz])

    def fields(self, discr):
        class ZField:
            shape = (3,)
            
            def __init__(self, zval):
                self.zval = zval

            def __call__(self, x, el):
                return numpy.array([0, 0, self.zval])

        from hedge.data import GivenFunction
        e = discr.volume_zeros(shape=(3,))
        h = 1/self.units.MU0*GivenFunction(ZField(self.bz)).volume_interpolant(discr)

        return e, h

    def positions(self, t):
        from math import sin, cos, pi
        return [
            numpy.array([
                self.gyration_rad*( sin(angle + self.omega_b*t)),
                self.gyration_rad*( cos(angle + self.omega_b*t)),
                self.vpar*t])
            for angle in numpy.arange(0, 2*pi, 2*pi/self.nparticles)
            ]

    def velocities(self, t):
        from math import sin, cos, pi
        return [
            numpy.array([
                self.omega_b*self.gyration_rad*( cos(angle + self.omega_b*t)),
                self.omega_b*self.gyration_rad*(-sin(angle + self.omega_b*t)),
                self.vpar])
            for angle in numpy.arange(0, 2*pi, 2*pi/self.nparticles)
            ]




class EBParallel(StaticFieldSetup):
    """Implements Problem 12.6b) in the third edition of Jackson."""
    def __init__(self, units, mass, charge, c, ez, bz, radius, nparticles):
        self.units = units
        self.ez = ez
        self.bz = bz
        StaticFieldSetup.__init__(self, mass, charge, c)

        if nparticles == 1:
            self.particle_offsets = [numpy.zeros((3,))]
        elif nparticles == 4:
            self.particle_offsets = [
                    numpy.array([-0.1, -0.1, 0]),
                    numpy.array([ 0.1, -0.1, 0]),
                    numpy.array([-0.1,  0.1, 0]),
                    numpy.array([ 0.1,  0.1, 0])
                    ]
        else:
            raise ValueError, "invalid number of particles"

        self.shift = numpy.zeros((3,))

        c = self.c = self.units.VACUUM_LIGHT_SPEED
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
        return numpy.array([0,0,self.ez])

    def h(self):
        return 1/self.units.MU0*numpy.array([0, 0, self.bz])

    def fields(self, discr):
        class ZField:
            shape = (3,)
            
            def __init__(self, zval):
                self.zval = zval

            def __call__(self, x, el):
                return numpy.array([0,0,self.zval])

        from hedge.data import GivenFunction
        e = GivenFunction(ZField(self.ez)).volume_interpolant(discr)
        h = 1/self.units.MU0*GivenFunction(ZField(self.bz)).volume_interpolant(discr)

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
        pos = (numpy.array([
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
                numpy.array([
                    self.scaling*self.R*cos(phi)*phi_prime,
                    self.scaling*self.R*(-sin(phi))*phi_prime,
                    zconst*sinh(self.rho*phi)*self.rho*phi_prime
                    ])
                ]




def run_setup(units, casename, setup, discr, pusher, visualize=False):
    from hedge.timestep import RK4TimeStepper
    from hedge.visualization import SiloVisualizer
    from hedge.pde import MaxwellOperator

    vis = SiloVisualizer(discr)

    from pyrticle.cloud import ParticleCloud, FaceBasedElementFinder
    from pyrticle.reconstruction import \
            ShapeFunctionReconstructor, \
            NormalizedShapeFunctionReconstructor
    cloud = ParticleCloud(discr, units, 
            ShapeFunctionReconstructor(),
            pusher(),
            FaceBasedElementFinder(),
            3, 3, debug=set(["verbose_vis"]))

    e, h = setup.fields(discr)
    b = units.MU0 * h

    init_positions = setup.positions(0)
    init_velocities = setup.velocities(0)

    nparticles = len(init_positions)
    cloud.add_particles( 
            zip(init_positions, init_velocities, 
                nparticles * [setup.charge],
                nparticles  * [units.EL_MASS],
                ),
            nparticles)

    final_time = setup.final_time()
    nsteps = setup.nsteps()
    dt = final_time/nsteps

    # timestepping ------------------------------------------------------------
    def rhs(t, y):
        return cloud.rhs(t, e, b)

    stepper = RK4TimeStepper()
    from time import time
    t = 0

    bbox = discr.mesh.bounding_box()
    z_period = bbox[1][2] - bbox[0][2]

    def check_result():
        from hedge.tools import cross

        deriv_dt = 1e-12

        dim = discr.dimensions
        true_x = setup.positions(t)
        true_v = setup.velocities(t)
        true_f = [(p2-p1)/(2*deriv_dt)
                for p1, p2 in zip(setup.momenta(t-deriv_dt), setup.momenta(t+deriv_dt))]

        from pyrticle.tools import NumberShiftableVector
        vis_info = cloud.vis_listener.particle_vis_map
        sim_x = cloud.positions
        sim_v = cloud.velocities()
        sim_f = NumberShiftableVector.unwrap(
                vis_info["mag_force"] + vis_info["el_force"])
        sim_el_f = NumberShiftableVector.unwrap(vis_info["el_force"])
        sim_mag_f = NumberShiftableVector.unwrap(vis_info["mag_force"])

        local_e = setup.e()
        local_b = units.MU0 * setup.h()

        x_err = 0
        v_err = 0
        f_err = 0

        for i in range(len(cloud)):
            #real_f = numpy.array(cross(sim_v, setup.charge*local_b)) + setup.charge*local_e

            my_true_x = true_x[i]
            my_true_x[2] = my_true_x[2] % z_period
            if False and i == 0:
                #print "particle %d" % i
                print "pos:", la.norm(true_x[i]-sim_x[i])/la.norm(true_x[i])
                #print "vel:", la.norm(true_v[i]-sim_v[i])/la.norm(true_v[i])
                #print "force:", la.norm(true_f[i]-sim_f[i])/la.norm(true_f[i])
                print "pos:", true_x[i], sim_x[i]
                #print "vel:", true_v[i], sim_v[i]
                #print "force:", true_f[i], sim_f[i]
                #print "forces%d:..." % i, sim_el_f[i], sim_mag_f[i]
                #print "acc%d:" % i, la.norm(a-sim_a)
                #u = numpy.vstack((v, sim_v, f, sim_f, real_f))
                #print "acc%d:\n%s" % (i, u)
                #raw_input()

            def rel_err(sim, true):
                return la.norm(true-sim)/la.norm(true)

            x_err = max(x_err, rel_err(sim_x[i], my_true_x))
            v_err = max(v_err, rel_err(sim_v[i], true_v[i]))
            f_err = max(f_err, rel_err(sim_f[i], true_f[i]))

        return x_err, v_err, f_err

    # make sure verbose-vis fields are filled
    rhs(t, cloud)

    errors = (0, 0, 0)

    for step in xrange(nsteps):
        if step % int(setup.nsteps()/300) == 0:
            errors = tuple(
                    max(old_err, new_err) 
                    for old_err, new_err in zip(errors, check_result()))

            if visualize:
                visf = vis.make_file("%s-%04d" % (casename, step))

                cloud.add_to_vis(vis, visf, time=t, step=step)

                if True:
                    vis.add_data(visf, [ ("e", e), ("h", h), ],
                            time=t, step=step)
                else:
                    vis.add_data(visf, [], time=t, step=step)
                visf.close()

        cloud.upkeep()
        cloud = stepper(cloud, t, dt, rhs)

        t += dt

    assert errors[0] < 2e-12, casename+"-pos"
    assert errors[1] < 2e-13, casename+"-vel"
    assert errors[2] < 2e-4, casename+"-acc"

    vis.close()





