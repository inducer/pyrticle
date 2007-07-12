from __future__ import division
import pylinear.array as num
import pylinear.computation as comp




def gaussian_vectors(count, dim, center, sigma):
    from random import normalvariate
    return [num.array([ci + normalvariate(0, sigma) for ci in center])
            for i in range(count)]




def main():
    from hedge.element import TetrahedralElement
    from hedge.timestep import RK4TimeStepper
    from hedge.mesh import  make_box_mesh
    from hedge.discretization import \
            Discretization, \
            bind_flux, \
            bind_nabla, \
            bind_mass_matrix, \
            bind_inverse_mass_matrix, \
            pair_with_boundary
    from hedge.visualization import SiloVisualizer
    from hedge.silo import DB_VARTYPE_VECTOR
    from hedge.tools import dot, cross
    from math import sqrt, pi
    from pytools.arithmetic_container import \
            ArithmeticList, concatenate_fields
    from hedge.operators import MaxwellOperator
    from pyrticle.pointcloud import PointCloud

    epsilon0 = 8.8541878176e-12 # C**2 / (N m**2)
    mu0 = 4*pi*1e-7 # N/A**2.
    epsilon = 1*epsilon0
    mu = 1*mu0
    c = 1/sqrt(epsilon*mu)

    el_mass = 9.10938215e-31 # kg
    el_charge = -1.602176487e-19 # C

    #mesh = make_cylinder_mesh(radius=1, height=2, max_volume=0.01)
    mesh = make_box_mesh(max_volume=0.01)

    discr = Discretization(mesh, TetrahedralElement(order))
    vis = SiloVisualizer(discr)

    print "%d elements" % len(discr.mesh.elements)

    dt = discr.dt_factor(1/sqrt(mu*epsilon))
    final_time = dt*60
    nsteps = int(final_time/dt)+1
    dt = final_time/nsteps

    print "dt", dt
    print "nsteps", nsteps

    mass = bind_mass_matrix(discr)

    def l2_norm(field):
        return sqrt(dot(field, mass*field))

    maxwell = MaxwellOperator(discr, epsilon, mu, upwind_alpha=0)

    nparticles = 1
    cloud = PointCloud(discr)
    zero = num.zeros((discr.dimensions,))
    cloud.add_points(
            gaussian_vectors(nparticles, discr.dimensions, zero, 0.1),
            gaussian_vectors(nparticles, discr.dimensions, zero, 0.5*c),
            el_charge, el_mass)


    fields = concatenate_fields(
            [discr.volume_zeros() for i in range(2*discr.dimensions)],
            [cloud])

    def rhs(t, y):
        e = y[:3]
        h = y[3:6]
        
        rho, j = cloud.reconstruct_densities()

        maxwell_rhs = maxwell.rhs(t, y[0:6])
        rhs_e = maxwell_rhs[:3]
        rhs_h = maxwell_rhs[3:6]
        return concatenate_fields(
                rhs_e + 1/epsilon*j,
                rhs_h,
                [cloud.rhs(t, e, h)]
                )

    stepper = RK4TimeStepper()
    from time import time
    last_tstep = time()
    t = 0

    for step in range(nsteps):
        print "timestep %d, t=%f l2[e]=%g l2[h]=%g secs=%f" % (
                step, t, l2_norm(fields[0:3]), l2_norm(fields[3:6]),
                time()-last_tstep)
        last_tstep = time()

        db = vis("pic-%04d.silo" % step,
                vectors=[("e", fields[0:3]), 
                    ("h", fields[3:6]), ],
                expressions=[
                    ],
                write_coarse_mesh=True,
                time=t, step=step
                )
        cloud.add_to_silo_db(db)
        del db

        fields = stepper(fields, t, dt, rhs)
        cloud.upkeep()

        t += dt

if __name__ == "__main__":
    #import cProfile as profile
    #profile.run("main()", "wave2d.prof")
    main()
