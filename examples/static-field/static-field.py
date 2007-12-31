from __future__ import division
import pylinear.array as num
import pylinear.computation as comp
import pylinear.operator as op
import cProfile as profile
from pyrticle.units import SI

units = SI()




class StaticFieldSetup:
    def __init__(self, mass, charge, c):
        self.mass = mass
        self.charge = charge
        self.c = c

    def gamma(self, v):
        value = (1-comp.norm_2_squared(v)/self.c**2)**(-0.5)
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
                #u = num.vstack((vpos, v, (vpos-v)/comp.norm_2(v)))
                #print (vpos - v)/comp.norm_2(v)
                assert comp.norm_2(vpos - v)/comp.norm_2(v) < threshold





class LarmorScrew(StaticFieldSetup):
    """Implements Section 12.2 in the third edition of Jackson."""
    def __init__(self, mass, charge, c, vpar, vperp, bz, nparticles):
        from math import pi 

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
    def __init__(self, mass, charge, c, ez, bz, radius, nparticles):
        self.ez = ez
        self.bz = bz
        StaticFieldSetup.__init__(self, mass, charge, c)

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




def run_setup(casename, setup, discr):
    from hedge.timestep import RK4TimeStepper
    from hedge.visualization import VtkVisualizer, SiloVisualizer
    from hedge.operators import MaxwellOperator

    vis = SiloVisualizer(discr)
    #vis = VtkVisualizer(discr, "pic")

    from pyrticle.cloud import ParticleCloud
    from pyrticle.reconstruction import ShapeFunctionReconstructor
    from pyrticle.pusher import MonomialParticlePusher
    cloud = ParticleCloud(discr, units, 
            ShapeFunctionReconstructor(),
            MonomialParticlePusher(),
            3, 3, verbose_vis=True)

    e, h = setup.fields(discr)
    b = units.MU0 * h

    init_positions = setup.positions(0)
    init_velocities = setup.velocities(0)

    #print "x", init_positions
    #print "v", init_velocities

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
        return cloud.rhs(t, e, b)

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
        all_sim_v = cloud.velocities()
        all_f = [(p2-p1)/(2*deriv_dt)
                for p1, p2 in zip(setup.momenta(t-deriv_dt), setup.momenta(t+deriv_dt))]

        all_sim_f = cloud.vis_info["lorentz_force"] + cloud.vis_info["el_force"]

        local_e = setup.e()
        local_b = units.MU0 * setup.h()

        x_err = 0
        v_err = 0
        f_err = 0

        for i in range(len(cloud)):
            x = all_x[i]
            sim_x = cloud.positions[i*dim:(i+1)*dim]
            v = all_v[i]
            sim_v = cloud.velocities()[i*dim:(i+1)*dim]
            f = all_f[i]
            sim_f = all_sim_f[i*dim:(i+1)*dim]

            real_f = num.array(cross(sim_v, setup.charge*local_b)) + setup.charge*local_e

            #print "pos%d:" % i, comp.norm_2(x-sim_x)

            #print "vel%d:" % i, comp.norm_2(v-sim_v)
            #print "vel%d:..." % i, v, sim_v

            #print "acc%d:" % i, comp.norm_2(a-sim_a)
            #u = num.vstack((v, sim_v, f, sim_f, real_f))
            #print "acc%d:\n%s" % (i, u)
            #raw_input()

            x_err = max(x_err, comp.norm_2(v-sim_v)/comp.norm_2(v))
            v_err = max(v_err, comp.norm_2(v-sim_v)/comp.norm_2(v))
            f_err = max(f_err, comp.norm_2(f-sim_f)/comp.norm_2(f))

        return x_err, v_err, f_err

    # make sure verbose-vis fields are filled
    rhs(t, cloud)

    errors = (0, 0, 0)

    for step in xrange(nsteps):
        if step % (setup.nsteps()/100) == 0:
            print "timestep %d, t=%g secs=%f particles=%d" % (
                    step, t, time()-last_tstep, len(cloud))

            errors = tuple(
                    max(old_err, new_err) 
                    for old_err, new_err in zip(errors, check_result()))

            last_tstep = time()
            if True:
                visf = vis.make_file("%s-%04d" % (casename, step))

                mesh_scalars, mesh_vectors = \
                        cloud.add_to_vis(vis, visf, time=t, step=step)

                if False:
                    vis.add_data(visf, [ ("e", e), ("h", h), ]
                            + mesh_scalars + mesh_vectors,
                            write_coarse_mesh=True,

                            time=t, step=step)
                else:
                    vis.add_data(visf, [], time=t, step=step)
                visf.close()

        cloud.upkeep()
        cloud = stepper(cloud, t, dt, rhs)

        t += dt

    print
    print "l_inf errors (pos,vel,acc):", errors

    vis.close()




def main():
    from hedge.element import TetrahedralElement
    from hedge.mesh import \
            make_box_mesh, \
            make_cylinder_mesh
    from hedge.discretization import Discretization

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

    # particles setup ---------------------------------------------------------
    def get_setup(case):
        c = units.VACUUM_LIGHT_SPEED
        if case == "screw":
            return LarmorScrew(mass=units.EL_MASS, charge=units.EL_CHARGE, c=c,
                    vpar=c*0.8, vperp=c*0.1, bz=1e-3, 
                    nparticles=4)
        elif case == "epb":
            return EBParallel(mass=units.EL_MASS, charge=units.EL_CHARGE, c=c,
                    ez=1e+5, bz=1e-3, radius=0.5*radius, nparticles=1)
        else:
            raise ValueError, "invalid test case"

    for case in ["screw", "epb"]:
        print "----------------------------------------------"
        print case.upper()
        print "----------------------------------------------"
        run_setup(case, get_setup(case), discr)





if __name__ == "__main__":
    #profile.run("main()", "pic.prof")
    main()
