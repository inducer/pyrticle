"""Main state container Python interface"""

__copyright__ = "Copyright (C) 2007, 2008 Andreas Kloeckner"

__license__ = """
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see U{http://www.gnu.org/licenses/}.
"""




import pylinear.array as num
import pylinear.computation as comp
from pytools.arithmetic_container import work_with_arithmetic_containers
from pytools.arithmetic_container import ArithmeticList
from math import sqrt
import pyrticle._internal as _internal

from pyrticle.meshdata import MeshData




class MapStorageVisualizationListener(_internal.VisualizationListener):
    def __init__(self, storage_map, particle_number_shift_signaller):
        _internal.VisualizationListener.__init__(self)
        self.storage_map = storage_map
        self.particle_number_shift_signaller = particle_number_shift_signaller

    def store_particle_vis_vector(self, name, vec, entries_per_particle):
        from pyrticle.tools import NumberShiftableVector
        self.storage_map[name] = NumberShiftableVector(vec, 
                multiplier=entries_per_particle,
                signaller=self.particle_number_shift_signaller)




class ParticleCloud:
    """State container for a cloud of particles. Supports particle
    problems of any dimension, examples below are given for three
    dimensions for simplicity.
    """
    def __init__(self, discr, units, 
            reconstructor, pusher,
            dimensions_pos, dimensions_velocity,
            verbose_vis=False):

        self.units = units
        self.discretization = discr
        self.verbose_vis = verbose_vis

        self.reconstructor = reconstructor
        self.pusher = pusher

        self.dimensions_mesh = discr.dimensions
        self.dimensions_pos = dimensions_pos
        self.dimensions_velocity = dimensions_velocity

        dims = (dimensions_pos, dimensions_velocity)

        self.pic_algorithm = getattr(_internal, 
                "PIC%s%s%d%d" % (
                    self.reconstructor.name,
                    self.pusher.name,
                    self.dimensions_pos,
                    self.dimensions_velocity
                    ),
                )(discr.dimensions, units.VACUUM_LIGHT_SPEED)

        # We need to retain this particular Python wrapper
        # of our mesh_data object because we write new stuff
        # to its __dict__. (And every time we access it via
        # pic_algorithm, a new wrapper--with a fresh __dict__--gets
        # created, which is not what we want.)
        self.mesh_data = self.pic_algorithm.mesh_data
        self.mesh_data.fill_from_hedge(discr)

        # size change messaging
        from pyrticle.tools import NumberShiftForwarder
        self.particle_number_shift_signaller = NumberShiftForwarder()
        self.pic_algorithm.particle_number_shift_listener = self.particle_number_shift_signaller

        # visualization
        self.vis_info = {}
        self.pic_algorithm.vis_listener = MapStorageVisualizationListener(
                    self.vis_info,
                    self.particle_number_shift_signaller)

        # subsystem init
        self.reconstructor.initialize(self)
        self.pusher.initialize(self)

        self.derived_quantity_cache = {}

        # instrumentation 
        from pytools.log import IntervalTimer, EventCounter

        self.reconstruct_timer = IntervalTimer(
                "t_reconstruct",
                "Time spent reconstructing")
        self.reconstruct_counter = EventCounter(
                "n_reconstruct",
                "Number of reconstructions")

        self.find_el_timer = IntervalTimer(
                "t_find",
                "Time spent finding new elements")
        self.find_same_counter = EventCounter(
                "n_find_same",
                "#Particles found in same element")
        self.find_by_neighbor_counter = EventCounter(
                "n_find_neighbor",
                "#Particles found through neighbor")
        self.find_by_vertex_counter = EventCounter(
                "n_find_by_vertex",
                "#Particles found by vertex")
        self.find_global_counter = EventCounter(
                "n_find_global",
                "#Particles found by global search")

        self.force_timer = IntervalTimer(
                "t_force",
                "Time spent calculating forces")

    def __len__(self):
        return self.pic_algorithm.particle_count

    def add_instrumentation(self, mgr):
        mgr.add_quantity(self.reconstruct_timer)
        mgr.add_quantity(self.reconstruct_counter)

        mgr.add_quantity(self.find_el_timer)
        mgr.add_quantity(self.find_same_counter)
        mgr.add_quantity(self.find_by_neighbor_counter)
        mgr.add_quantity(self.find_by_vertex_counter)
        mgr.add_quantity(self.find_global_counter)

        mgr.add_quantity(self.force_timer)

        self.reconstructor.add_instrumentation(mgr)

        mgr.set_constant("reconstructor", self.reconstructor.name)
        mgr.set_constant("pusher", self.pusher.name)

    @property
    def positions(self):
        from hedge.tools import FixedSizeSliceAdapter
        return FixedSizeSliceAdapter(
                self.pic_algorithm.positions, 
                self.dimensions_pos,
                self.pic_algorithm.particle_count)

    @property
    def momenta(self):
        from hedge.tools import FixedSizeSliceAdapter
        return FixedSizeSliceAdapter(
                self.pic_algorithm.momenta, 
                self.dimensions_velocity,
                self.pic_algorithm.particle_count)

    @property
    def masses(self):
        return self.pic_algorithm.masses[:self.pic_algorithm.particle_count]

    @property
    def charges(self):
        return self.pic_algorithm.charges[:self.pic_algorithm.particle_count]

    def raw_velocities(self):
        if "raw_velocities" in self.derived_quantity_cache:
            return self.derived_quantity_cache["raw_velocities"]
        else:
            result = self.pic_algorithm.velocities()
            self.derived_quantity_cache["raw_velocities"] = result
            return result

    def velocities(self):
        from hedge.tools import FixedSizeSliceAdapter
        return FixedSizeSliceAdapter(
                self.raw_velocities(),
                self.dimensions_velocity)

    def add_particles(self, positions, velocities, charges, masses):
        """Add the particles with the given data to the cloud."""

        new_count = len(positions)

        # expand scalar masses and charges
        try:
            len(charges)
        except:
            charges = charges*num.ones((new_count,))

        try:
            len(masses)
        except:
            masses = new_count * [masses]

        # convert velocities to momenta
        momenta = [m*self.units.gamma(v)*v for m, v in zip(masses, velocities)]

        # find containing elements
        pic = self.pic_algorithm
        containing_elements = [pic.mesh_data.find_containing_element(p) 
                for p in positions]

        # weed out uncontained particles
        deathflags = [ce == MeshData.INVALID_ELEMENT 
                for ce in containing_elements]
        containing_elements = [ce for ce in containing_elements 
                if ce != MeshData.INVALID_ELEMENT]
        positions = [p for p, dead in zip(positions, deathflags) 
                if not dead]
        momenta = [v for v, dead in zip(momenta, deathflags) 
                if not dead]
        charges = num.array([c for c, dead in zip(charges, deathflags) 
            if not dead])
        masses = num.array([m for m, dead in zip(masses, deathflags) 
            if not dead])

        # check vector dimensionalities
        xdim = self.dimensions_pos
        vdim = self.dimensions_velocity
        for x in positions:
            assert len(x) == xdim
        for p in momenta:
            assert len(p) == vdim

        # add particles

        pic.containing_elements[pic.particle_count:] = containing_elements

        prev_count = pic.particle_count
        pic.positions = num.hstack(
                (pic.positions[:prev_count*xdim], num.hstack(positions))
                )
        pic.momenta = num.hstack(
                (pic.momenta[:prev_count*vdim], num.hstack(momenta)))

        pic.charges = num.hstack(
                (pic.charges[:prev_count], charges))
        pic.masses = num.hstack(
                (pic.masses[:prev_count], masses))

        pic.particle_count += len(containing_elements)

        self.check_containment()

        for pn in range(prev_count, pic.particle_count):
            self.reconstructor.add_particle_hook(pn)

    def check_containment(self):
        """Check that a containing element is known for each particle.

        This is a new invariant as of 1/17/08, and violations of this end
        segfaulting, which we should avoid.
        """
        for ce in self.pic_algorithm.containing_elements[:len(self)]:
            assert ce != MeshData.INVALID_ELEMENT
                
    def upkeep(self):
        """Perform any operations must fall in between timesteps,
        such as resampling or deleting particles.
        """
        self.vis_info.clear()
        self.pic_algorithm.perform_reconstructor_upkeep()

    def reconstruct_densities(self):
        """Return a tuple (charge_density, current_densities), where
        current_densities is an 
          ArithmeticList([[jx0,jx1,...],[jy0,jy1,...]])  
        of the densities in each direction.
        """

        if "j" in self.derived_quantity_cache:
            j = self.derived_quantity_cache["j"]
            if "rho" in self.derived_quantity_cache:
                rho = self.derived_quantity_cache["rho"]
            else:
                rho = self.reconstruct_rho()
            return rho, j
        else:
            if "rho" in self.derived_quantity_cache:
                rho = self.derived_quantity_cache["rho"]
                j = self.reconstruct_j(self.raw_velocities())
                return rho, j
            else:
                rho = self.discretization.volume_zeros()
                from hedge.tools import FixedSizeSliceAdapter
                j = FixedSizeSliceAdapter(num.zeros(
                    (self.dimensions_velocity*len(self.discretization),)),
                    self.dimensions_velocity)

                self.reconstruct_timer.start()
                self.reconstructor.reconstruct_hook()
                self.pic_algorithm.reconstruct_densities(rho, j, self.raw_velocities())
                j = j.get_alist_of_components()
                self.reconstruct_timer.stop()
                self.reconstruct_counter.add(self.dimensions_velocity+1)

                self.derived_quantity_cache["rho"] = rho
                self.derived_quantity_cache["j"] = j

                return rho, j

    def reconstruct_j(self):
        """Return a the current densities in the structure::
          ArithmeticList([[jx0,jx1,...],[jy0,jy1,...]])  
        """

        if "j" in self.derived_quantity_cache:
            return self.derived_quantity_cache["j"]

        from hedge.tools import FixedSizeSliceAdapter
        j = FixedSizeSliceAdapter(num.zeros(
            (self.dimensions_velocity*len(self.discretization),)),
            self.dimensions_velocity)

        self.reconstruct_timer.start()
        self.reconstructor.reconstruct_hook()
        self.pic_algorithm.reconstruct_j(j.adaptee, self.raw_velocities())
        j = j.get_alist_of_components()
        self.reconstruct_timer.stop()
        self.reconstruct_counter.add(self.dimensions_velocity)

        self.derived_quantity_cache["j"] = j

        return j

    def reconstruct_rho(self):
        """Return a the charge_density as a volume vector.
        """

        if "rho" in self.derived_quantity_cache:
            return self.derived_quantity_cache["rho"]

        rho = self.discretization.volume_zeros()

        self.reconstruct_timer.start()
        self.reconstructor.reconstruct_hook()
        self.pic_algorithm.reconstruct_rho(rho)
        self.reconstruct_timer.stop()
        self.reconstruct_counter.add(1)

        self.derived_quantity_cache["rho"] = rho

        return rho

    def rhs(self, t, e, b):
        """Return an ArithmeticList of velocities and forces on the particles.

        @arg e: triple of M{E_x}, M{E_y}, M{E_z}, each of which may be either 
          a Pylinear vector or a L{ZeroVector}.
        @arg b: triple of M{B_x}, M{B_y}, M{B_z}, each of which may be either 
          a Pylinear vector or a L{ZeroVector}. Caution: The hedge Maxwell operator
          deals with M{H}, not M{B}.
        """

        velocities = self.raw_velocities()

        field_args = tuple(e) + tuple(b)

        # compute forces
        self.force_timer.start()
        forces = self.pic_algorithm.forces(
                velocities=velocities,
                verbose_vis=self.verbose_vis,
                *field_args
                )
        self.force_timer.stop()

        from pyrticle.tools import NumberShiftableVector
        result = ArithmeticList([
            NumberShiftableVector(velocities, 
                multiplier=self.dimensions_pos,
                signaller=self.particle_number_shift_signaller),
            NumberShiftableVector(forces, 
                multiplier=self.dimensions_velocity,
                signaller=self.particle_number_shift_signaller),
            self.reconstructor.rhs()
            ])
        return result

    def __iadd__(self, rhs):
        self.derived_quantity_cache.clear()

        assert isinstance(rhs, ArithmeticList)

        from pyrticle.tools import NumberShiftableVector

        dx, dp, drecon = rhs
        dx = NumberShiftableVector.unwrap(dx)
        dp = NumberShiftableVector.unwrap(dp)
        assert len(dx) == self.dimensions_pos*len(self)
        assert len(dp) == self.dimensions_velocity*len(self)
        self.pic_algorithm.add_rhs(dx, dp)
        self.reconstructor.add_rhs(drecon)

        self.find_el_timer.start()
        self.pic_algorithm.update_containing_elements()
        self.find_el_timer.stop()
        self.find_same_counter.transfer(
                self.pic_algorithm.find_same)
        self.find_by_neighbor_counter.transfer(
                self.pic_algorithm.find_by_neighbor)
        self.find_by_vertex_counter.transfer(
                self.pic_algorithm.find_by_vertex)
        self.find_global_counter.transfer(
                self.pic_algorithm.find_global)

        return self

    def add_to_vis(self, visualizer, vis_file, time=None, step=None, verbose=None):
        if verbose is None:
            verbose = self.verbose_vis

        from hedge.visualization import VtkVisualizer, SiloVisualizer
        if isinstance(visualizer, VtkVisualizer):
            return self._add_to_vtk(visualizer, vis_file, time, step, verbose)
        elif isinstance(visualizer, SiloVisualizer):
            return self._add_to_silo(vis_file, time, step, verbose)
        else:
            raise ValueError, "unknown visualizer type `%s'" % type(visualizer)

    def _add_to_vtk(self, visualizer, vis_file, time, step, verbose):
        from hedge.vtk import \
                VTK_VERTEX, VF_INTERLEAVED, \
                DataArray, \
                UnstructuredGrid, \
                AppendedDataXMLGenerator

        dim = self.pic_algorithm.mesh_data.dimensions

        points = (len(self.pic_algorithm.containing_elements),
                DataArray("points", self.pic_algorithm.positions,
                    vector_format=VF_INTERLEAVED,
                    components=dim))
        grid = UnstructuredGrid(
                points, 
                cells=[[i] for i in range(points[0])],
                cell_types=[VTK_VERTEX] * points[0],
                )

        grid.add_pointdata(
                DataArray("velocity", 
                    self.pic_algorithm.velocities(),
                    vector_format=VF_INTERLEAVED,
                    components=dim)
                )

        def add_vis_vector(name):
            if name in self.vis_info:
                vec = self.vis_info[name]
            else:
                vec = num.zeros((len(self.pic_algorithm.containing_elements) * dim,))

            grid.add_pointdata(
                    DataArray(name, 
                        vec, vector_format=VF_INTERLEAVED,
                        components=dim)
                    )

        if verbose:
            add_vis_vector("pt_e")
            add_vis_vector("pt_b")
            add_vis_vector("el_force")
            add_vis_vector("mag_force")

        from os.path import splitext
        pathname = splitext(vis_file.pathname)[0] + "-particles.vtu"
        from pytools import assert_not_a_file
        assert_not_a_file(pathname)

        outf = open(pathname, "w")
        AppendedDataXMLGenerator()(grid).write(outf)
        outf.close()

        visualizer.register_pathname(time, pathname)

    def _add_to_silo(self, db, time, step, verbose):
        from pylo import DBOPT_DTIME, DBOPT_CYCLE
        from warnings import warn

        optlist = {}
        if time is not None:
            optlist[DBOPT_DTIME] = time
        if step is not None:
            optlist[DBOPT_CYCLE] = step

        pcount = len(self)

        if pcount:
            db.put_pointmesh("particles", self.dimensions_pos, 
                    num.hstack(self.positions.get_alist_of_components()), optlist)
            db.put_pointvar("momenta", "particles", 
                    self.momenta.get_alist_of_components())
            db.put_pointvar("velocity", "particles", 
                    self.velocities().get_alist_of_components())

            for name, value in self.vis_info.iteritems():
                from pyrticle.tools import NumberShiftableVector
                value = NumberShiftableVector.unwrap(value)
                dim, remainder = divmod(len(value), pcount)
                assert remainder == 0, (
                        "particle vis value '%s' had invalid number of entries: "
                        "%d (#particles=%d)" % (name, len(value), pcount))
                if dim == 1:
                    db.put_pointvar1(name, "particles", value)
                else:
                    db.put_pointvar(name, "particles", [value[i::dim] for i in range(dim)])




class FieldsAndCloud:
    def __init__(self, maxwell_op, e, h, cloud):
        from pytools.arithmetic_container import join_fields
        self.maxwell_op = maxwell_op
        self.e = e
        self.h = h
        self.cloud = cloud

        self.eh_components = maxwell_op.component_count()

        from pytools.log import IntervalTimer
        self.field_solve_timer = IntervalTimer(
                "t_field",
                "Time spent in field solver")

    def add_instrumentation(self, mgr):
        mgr.add_quantity(self.field_solve_timer)

        self.cloud.add_instrumentation(mgr)
        self.maxwell_op.discr.add_instrumentation(mgr)

    def __iadd__(self, other):
        assert len(other) == self.eh_components + 1
        rhs_e, rhs_h = self.maxwell_op.split_fields(other[:self.eh_components])
        rhs_cloud = other[self.eh_components]

        self.e += rhs_e
        self.h += rhs_h
        self.cloud += rhs_cloud

        return self

    def rhs(self, t, y):
        assert y is self

        from pytools.arithmetic_container import join_fields
        from pyrticle._internal import ZeroVector

        # assemble field_args of the form [ex,ey,ez] and [bx,by,bz],
        # inserting ZeroVectors where necessary.
        idx = 0
        e_arg = []
        for use_component in self.maxwell_op.get_subset()[0:3]:
            if use_component:
                e_arg.append(self.e[idx])
                idx += 1
            else:
                e_arg.append(ZeroVector())

        idx = 0
        b_arg = []
        for use_component in self.maxwell_op.get_subset()[3:6]:
            if use_component:
                b_arg.append(self.maxwell_op.mu * self.h[idx])
                idx += 1
            else:
                b_arg.append(ZeroVector())

        cloud_rhs = self.cloud.rhs(t, e_arg, b_arg)

        self.field_solve_timer.start()
        rhs_e, rhs_h = self.maxwell_op.split_fields(
                self.maxwell_op.rhs(t, join_fields(self.e, self.h))
                )
        self.field_solve_timer.stop()

        return join_fields(
                rhs_e - 1/self.maxwell_op.epsilon*self.cloud.reconstruct_j(),
                rhs_h,
                ).plus([cloud_rhs])




# shape bandwidth -------------------------------------------------------------
def guess_shape_bandwidth(cloud):
    self.cloud.reconstructor.set_radius(mesh_data.advisable_particle_radius())




def optimize_shape_bandwidth(cloud, discr, analytic_rho, rhovis=False, plot_l1_errors=False):
    adv_radius = cloud.mesh_data.advisable_particle_radius()
    radii = [adv_radius*2**exponent 
            for exponent in num.linspace(-4, 2, 50)]

    if rhovis:
        from hedge.visualization import SiloVisualizer
        vis = SiloVisualizer(discr)

    from hedge.discretization import integral

    tried_radii = []
    l1_errors = []
    for step, radius in enumerate(radii):
        cloud.reconstructor.set_radius(radius)
        cloud.derived_quantity_cache.clear()
        try:
            rec_rho = cloud.reconstruct_rho()
        except RuntimeError:
            continue

        tried_radii.append(radius)
        l1_errors.append(integral(discr, num.abs(rec_rho-analytic_rho)))

        if rhovis:
            visf = vis.make_file("rho-%04d" % step)
            cloud.add_to_vis(vis, visf, time=radius, step=step)
            vis.add_data(visf, [ 
                ("rho", rec_rho), 
                ("anarho", analytic_rho), 
                ],
                time=radius, step=step)
            visf.close()

    if rhovis:
        vis.close()

    from pytools import argmin
    good_rad = tried_radii[argmin(l1_errors)]

    print "radius: guessed optimum=%g, found optimum=%g" % (adv_radius, good_rad)

    if plot_l1_errors:
        from pylab import plot, show
        plot(tried_radii, l1_errors)
        show()

    cloud.reconstructor.set_radius(good_rad)




# initial condition -----------------------------------------------------------
def compute_initial_condition(pcon, discr, cloud, mean_beta, max_op,
        debug=False, force_zero=False):
    from hedge.operators import WeakPoissonOperator
    from hedge.mesh import TAG_ALL, TAG_NONE
    from hedge.data import ConstantGivenFunction, GivenVolumeInterpolant
    from hedge.tools import dot

    def l2_norm(field):
        return sqrt(dot(field, discr.mass_operator*field))
    def rel_l2_error(field, true):
        err = l2_norm(field-true)
        if err == 0:
            return 0
        else:
            return err/l2_norm(true)

    gamma = (1-comp.norm_2_squared(mean_beta))**(-0.5)

    # see doc/notes.tm for derivation of IC

    def make_scaling_matrix(beta_scale, other_scale):
        if comp.norm_2(mean_beta) < 1e-10:
            return other_scale*num.identity(discr.dimensions) 
        else:
            beta_unit = mean_beta/comp.norm_2(mean_beta)
            return (other_scale*num.identity(discr.dimensions) 
                    + (beta_scale-other_scale)*(beta_unit <<num.outer>> beta_unit))

    poisson_op = WeakPoissonOperator(discr, 
            diffusion_tensor=ConstantGivenFunction(
                make_scaling_matrix(1/gamma**2, 1)),
            dirichlet_tag=TAG_ALL,
            neumann_tag=TAG_NONE,
            )

    rho_prime = cloud.reconstruct_rho() 
    rho_tilde = rho_prime/gamma

    if force_zero:
        phi_tilde = discr.volume_zeros()
    else:
        from hedge.tools import parallel_cg
        phi_tilde = -parallel_cg(pcon, -poisson_op, 
                poisson_op.prepare_rhs(
                    GivenVolumeInterpolant(discr, rho_tilde/max_op.epsilon)), 
                debug=True, tol=1e-10)

    from pytools.arithmetic_container import ArithmeticListMatrix
    ALM = ArithmeticListMatrix
    e_tilde = ALM(make_scaling_matrix(1/gamma, 1))*poisson_op.grad(phi_tilde)
    e_prime = ALM(make_scaling_matrix(1, gamma))*e_tilde
    h_prime = (1/max_op.mu)*gamma/max_op.c * max_op.e_cross(mean_beta, e_tilde)
    if debug:
        from hedge.discretization import ones_on_volume
        reconstructed_charge = ones_on_volume(discr)*(discr.mass_operator*rho_prime)
        real_charge = sum(cloud.charges)
        print "charge: supposed=%g reconstructed=%g error=%g %%" % (
                real_charge,
                reconstructed_charge,
                100*abs(reconstructed_charge-real_charge)/abs(real_charge)
                )

        from hedge.operators import DivergenceOperator

        div_op = DivergenceOperator(discr)

        d_tilde = max_op.epsilon*e_tilde
        d_prime = max_op.epsilon*e_prime

        divD_prime_ldg = poisson_op.div(d_prime)
        divD_prime_ldg2 = poisson_op.div(d_prime, max_op.epsilon*gamma*phi_tilde)
        divD_prime_ldg3 = max_op.epsilon*\
                (discr.inverse_mass_operator*poisson_op.op(gamma*phi_tilde))
        divD_prime_central = div_op(d_prime)

        print "l2 div D_prime error central: %g" % \
                rel_l2_error(divD_prime_central, rho_prime)
        print "l2 div D_prime error ldg: %g" % \
                rel_l2_error(divD_prime_ldg, rho_prime)
        print "l2 div D_prime error ldg with phi: %g" % \
                rel_l2_error(divD_prime_ldg2, rho_prime)
        print "l2 div D_prime error ldg with phi 3: %g" % \
                rel_l2_error(divD_prime_ldg3, rho_prime)

        from hedge.visualization import SiloVisualizer
        vis = SiloVisualizer(discr)
        visf = vis.make_file("ic")
        vis.add_data(visf, [ 
            #("phi_tilde", phi_tilde),
            #("rho_tilde", rho_tilde), 
            #("e_tilde", e_tilde), 

            ("rho_prime", rho_prime), 
            ("divD_prime_ldg", divD_prime_ldg),
            ("divD_prime_ldg2", divD_prime_ldg2),
            ("divD_prime_ldg3", divD_prime_ldg3),
            ("divD_prime_central", divD_prime_central),
            ("e_prime", e_prime), 
            ("h_prime", h_prime), 
            ],
            )
        cloud.add_to_vis(vis, visf, verbose=False)
        visf.close()

    return FieldsAndCloud(max_op, e_prime, h_prime, cloud)

