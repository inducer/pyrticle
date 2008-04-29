"""Main state container Python interface"""

from __future__ import division

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




import numpy
import numpy.linalg as la
import pyrticle._internal as _internal

from pyrticle.meshdata import MeshData




class MapStorageVisualizationListener(_internal.VisualizationListener):
    def __init__(self, particle_number_shift_signaller):
        _internal.VisualizationListener.__init__(self)
        self.mesh_vis_map = {}
        self.particle_vis_map = {}
        self.particle_number_shift_signaller = particle_number_shift_signaller

    def store_mesh_vis_vector(self, name, vec, components=1):
        self.mesh_vis_map[name] = numpy.reshape(
                vec, (components, len(vec)//components))

    def store_particle_vis_vector(self, name, vec, entries_per_particle=1):
        from pyrticle.tools import NumberShiftableVector
        # FIXME is this still correct?
        self.particle_vis_map[name] = NumberShiftableVector(vec, 
                multiplier=entries_per_particle,
                signaller=self.particle_number_shift_signaller)

    def clear(self):
        self.mesh_vis_map.clear()
        self.particle_vis_map.clear()




class ElementFinder(object):
    pass
class FaceBasedElementFinder(ElementFinder):
    name = "FaceBasedFind"
class HeuristicElementFinder(ElementFinder):
    name = "HeuristicFind"




class ParticleCloud:
    """State container for a cloud of particles. Supports particle
    problems of any dimension, examples below are given for three
    dimensions for simplicity.
    """
    def __init__(self, discr, units, 
            reconstructor, pusher, finder,
            dimensions_pos, dimensions_velocity,
            verbose_vis=False):

        self.units = units
        self.discretization = discr
        self.verbose_vis = verbose_vis

        self.reconstructor = reconstructor
        self.pusher = pusher
        self.finder = finder

        self.dimensions_mesh = discr.dimensions
        self.dimensions_pos = dimensions_pos
        self.dimensions_velocity = dimensions_velocity

        dims = (dimensions_pos, dimensions_velocity)

        self.pic_algorithm = getattr(_internal, 
                "PIC%s%s%s%d%d" % (
                    self.reconstructor.name,
                    self.pusher.name,
                    self.finder.name,
                    self.dimensions_pos,
                    self.dimensions_velocity
                    ),
                )(discr.dimensions, units.VACUUM_LIGHT_SPEED)
        self.pic_algorithm.positions = numpy.zeros((0,), dtype=float)
        self.pic_algorithm.momenta = numpy.zeros((0,), dtype=float)
        self.pic_algorithm.charges = numpy.zeros((0,), dtype=float)
        self.pic_algorithm.masses = numpy.zeros((0,), dtype=float)

        # We need to retain this particular Python wrapper
        # of our mesh_data object because we write new stuff
        # to its __dict__. (And every time we access it via
        # pic_algorithm, a new wrapper--with a fresh __dict__--gets
        # created, which is not what we want.)
        self.mesh_data = self.pic_algorithm.mesh_data
        self.mesh_data.fill_from_hedge(discr)

        # size change messaging
        from pyrticle.tools import NumberShiftMultiplexer
        self.particle_number_shift_signaller = NumberShiftMultiplexer()
        self.pic_algorithm.particle_number_shift_listener = self.particle_number_shift_signaller

        # visualization
        self.vis_listener = \
                self.pic_algorithm.vis_listener = \
                MapStorageVisualizationListener(
                        self.particle_number_shift_signaller)

        # subsystem init
        self.reconstructor.initialize(self)
        self.pusher.initialize(self)

        self.derived_quantity_cache = {}

        # instrumentation 
        from pytools.log import IntervalTimer, EventCounter

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

    def set_ignore_core_warnings(self, ignore):
        from pyrticle.tools import WarningForwarder, WarningIgnorer
        import pyrticle.tools
        del pyrticle.tools.warning_forwarder
        if ignore:
            pyrticle.tools.warning_forwarder = WarningIgnorer()
        else:
            pyrticle.tools.warning_forwarder = WarningForwarder()

    def __len__(self):
        return self.pic_algorithm.particle_count

    def add_instrumentation(self, mgr):
        mgr.add_quantity(self.find_el_timer)
        mgr.add_quantity(self.find_same_counter)
        mgr.add_quantity(self.find_by_neighbor_counter)
        mgr.add_quantity(self.find_by_vertex_counter)
        mgr.add_quantity(self.find_global_counter)

        self.reconstructor.add_instrumentation(mgr)
        self.pusher.add_instrumentation(mgr)

    @property
    def positions(self):
        allocated_p_count = len(self.pic_algorithm.containing_elements)
        return numpy.reshape(self.pic_algorithm.positions,
                (allocated_p_count, self.dimensions_pos)
                )[:self.pic_algorithm.particle_count]

    @property
    def momenta(self):
        allocated_p_count = len(self.pic_algorithm.containing_elements)
        return numpy.reshape(self.pic_algorithm.momenta,
                (allocated_p_count, self.dimensions_velocity)
                )[:self.pic_algorithm.particle_count]

    @property
    def masses(self):
        return self.pic_algorithm.masses[:self.pic_algorithm.particle_count]

    @property
    def charges(self):
        return self.pic_algorithm.charges[:self.pic_algorithm.particle_count]

    def velocities(self):
        if "velocities" in self.derived_quantity_cache:
            return self.derived_quantity_cache["velocities"]
        else:
            result = numpy.reshape(
                    self.pic_algorithm.velocities(),
                    (self.pic_algorithm.particle_count, self.dimensions_velocity))
            self.derived_quantity_cache["velocities"] = result
            return result

    def add_particles(self, positions, velocities, charges, masses):
        """Add the particles with the given data to the cloud."""

        new_count = len(positions)

        # expand scalar masses and charges
        try:
            len(charges)
        except:
            charges = charges*numpy.ones((new_count,))

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
        charges = numpy.array([c for c, dead in zip(charges, deathflags) 
            if not dead])
        masses = numpy.array([m for m, dead in zip(masses, deathflags) 
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
        pic.positions = numpy.hstack(
                (pic.positions[:prev_count*xdim], numpy.hstack(positions))
                )
        pic.momenta = numpy.hstack(
                (pic.momenta[:prev_count*vdim], numpy.hstack(momenta)))

        pic.charges = numpy.hstack(
                (pic.charges[:prev_count], charges))
        pic.masses = numpy.hstack(
                (pic.masses[:prev_count], masses))

        pic.particle_count += len(containing_elements)

        self.check_containment()
        pic.note_change_particle_count(pic.particle_count)

        self.derived_quantity_cache.clear()

    def clear_particles(self):
        self.pic_algorithm.particle_count = 0
        self.reconstructor.clear_particles()
        self.derived_quantity_cache.clear()

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
        self.vis_listener.clear()
        self.pic_algorithm.perform_reconstructor_upkeep()

    def reconstruct_densities(self):
        """Return a tuple (charge_density, current_densities), where
        current_densities is an d-by-n array, where d is the number 
        of velocity dimensions, and n is the discretization nodes.
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
                j = self.reconstruct_j(self.velocities())
                return rho, j
            else:
                self.reconstruct_timer.start()
                self.reconstructor.reconstruct_hook()
                self.reconstruct_timer.stop()
                self.reconstruct_counter.add(self.dimensions_velocity+1)
                rho_j

                j = numpy.asarray(j.T, order="C")

                self.derived_quantity_cache["rho"] = rho
                self.derived_quantity_cache["j"] = j

                return rho, j

    def reconstruct_j(self):
        """Return a the current densities as an d-by-n array, where d 
        is the number of velocity dimensions, and n is the number of 
        discretization nodes.
        """

        if "j" in self.derived_quantity_cache:
            return self.derived_quantity_cache["j"]

        j = self.reconstructor.reconstruct_j(self.velocities())
        j = numpy.asarray(j.T, order="C")
        self.derived_quantity_cache["j"] = j

        return j

    def reconstruct_rho(self):
        """Return a the charge_density as a volume vector.
        """

        if "rho" in self.derived_quantity_cache:
            return self.derived_quantity_cache["rho"]

        rho = self.pic_algorithm.reconstruct_rho()
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

        velocities = self.velocities()

        field_args = tuple(e) + tuple(b)

        # compute forces
        forces = self.pusher.forces(
                velocities,
                self.verbose_vis,
                *field_args
                )
        forces = numpy.reshape(forces,
                (len(self), self.dimensions_velocity))

        from pyrticle.tools import NumberShiftableVector
        from hedge.tools import join_fields
        result = join_fields(
            NumberShiftableVector(velocities, 
                multiplier=self.dimensions_pos,
                signaller=self.particle_number_shift_signaller),
            NumberShiftableVector(forces, 
                multiplier=self.dimensions_velocity,
                signaller=self.particle_number_shift_signaller),
            self.reconstructor.rhs()
            )
        return result

    def __add__(self, rhs):
        self.derived_quantity_cache.clear()

        from pyrticle.tools import NumberShiftableVector

        dx, dp, drecon = rhs
        dx = NumberShiftableVector.unwrap(dx)
        dp = NumberShiftableVector.unwrap(dp)
        assert dx.shape == (len(self), self.dimensions_pos)
        assert dp.shape == (len(self), self.dimensions_velocity)
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

    def get_mesh_vis_vars(self):
        return self.vis_listener.mesh_vis_map.items()

    def add_to_vis(self, visualizer, vis_file, time=None, step=None, beamaxis=None):
        from hedge.visualization import VtkVisualizer, SiloVisualizer
        if isinstance(visualizer, VtkVisualizer):
            return self._add_to_vtk(visualizer, vis_file, time, step)
        elif isinstance(visualizer, SiloVisualizer):
            return self._add_to_silo(visualizer, vis_file, time, step, beamaxis)
        else:
            raise ValueError, "unknown visualizer type `%s'" % type(visualizer)

    def _add_to_vtk(self, visualizer, vis_file, time, step):
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

        def add_particle_vis_vector(name):
            if name in self.vis_listener.particle_vis_map:
                vec = self.vis_listener.particle_vis_map[name]
            else:
                vec = numpy.zeros((len(self.pic_algorithm.containing_elements) * dim,))

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

    def _add_to_silo(self, visualizer, db, time, step, beamaxis):
        from pylo import DBOPT_DTIME, DBOPT_CYCLE
        from warnings import warn

        optlist = {}
        if time is not None:
            optlist[DBOPT_DTIME] = time
        if step is not None:
            optlist[DBOPT_CYCLE] = step

        pcount = len(self)

        if pcount:
            # real-space ------------------------------------------------------
            db.put_pointmesh("particles", self.dimensions_pos, 
                    numpy.asarray(self.positions.T, order="C"), optlist)
            db.put_pointvar1("charge", "particles", self.charges)
            db.put_pointvar1("mass", "particles", self.masses)
            db.put_pointvar("momentum", "particles", 
                    numpy.asarray(self.momenta.T, order="C"))
            db.put_pointvar("velocity", "particles", 
                    numpy.asarray(self.velocities().T, order="C"))

            for name, value in self.vis_listener.particle_vis_map.iteritems():
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
            
            # phase-space -----------------------------------------------------
            axes_names = ["x", "y", "z"]

            if beamaxis is not None:
                for axis in range(min(self.dimensions_pos, self.dimensions_velocity)):
                    if axis == beamaxis:
                        continue 

                    axname = axes_names[axis]

                    db.put_defvars("phasespace_%s" % axname, 
                            [
                            ("part_%s" % axname, "coord(particles)[%d]" % axis), 
                            ("part_%s_prime" % axname, 
                                "momentum[%d]/momentum[%d]" % (axis, beamaxis)), 
                            ("part_%s_momentum" % axname, 
                                "momentum[%d]" % (axis)),
                            ])








class FieldsAndCloud:
    def __init__(self, maxwell_op, e, h, cloud):
        from pytools.arithmetic_container import join_fields
        self.maxwell_op = maxwell_op
        self.em_fields = maxwell_op.assemble_fields(e=e, h=h)
        self.cloud = cloud

        self.em_filters = []

        from pytools.log import IntervalTimer
        self.field_solve_timer = IntervalTimer(
                "t_field",
                "Time spent in field solver")

    def add_em_filter(self, filt):
        self.em_filters.append(filt)

    @property
    def e(self):
        e, h = self.maxwell_op.split_eh(self.em_fields)
        return e

    @property
    def h(self):
        e, h = self.maxwell_op.split_eh(self.em_fields)
        return h

    @property
    def phi(self):
        from pyrticle.hyperbolic import CleaningMaxwellOperator
        if isinstance(self.maxwell_op, CleaningMaxwellOperator):
            e, h, phi = self.maxwell_op.split_ehphi(self.em_fields)
            return phi
        else:
            return self.maxwell_op.discr.volume_zeros()

    def add_instrumentation(self, mgr):
        mgr.add_quantity(self.field_solve_timer)

        self.cloud.add_instrumentation(mgr)
        self.maxwell_op.discr.add_instrumentation(mgr)

    def __add__(self, other):
        d_em_fields, d_cloud = other

        self.em_fields += d_em_fields
        self.cloud += d_cloud

        return self

    def rhs(self, t, y):
        assert y is self

        from pyrticle._internal import ZeroVector

        # assemble field_args of the form [ex,ey,ez] and [bx,by,bz],
        # inserting ZeroVectors where necessary.
        idx = 0
        e_arg = []
        for use_component in self.maxwell_op.get_eh_subset()[0:3]:
            if use_component:
                e_arg.append(self.e[idx])
                idx += 1
            else:
                e_arg.append(ZeroVector())

        idx = 0
        b_arg = []
        for use_component in self.maxwell_op.get_eh_subset()[3:6]:
            if use_component:
                b_arg.append(self.maxwell_op.mu * self.h[idx])
                idx += 1
            else:
                b_arg.append(ZeroVector())

        rhs_cloud = self.cloud.rhs(t, e_arg, b_arg)

        # calculate EM right-hand side 
        self.field_solve_timer.start()
        from pyrticle.hyperbolic import CleaningMaxwellOperator
        if isinstance(self.maxwell_op, CleaningMaxwellOperator):
            rhs_em = self.maxwell_op.rhs(t, self.em_fields, 
                    self.cloud.reconstruct_rho())
        else:
            rhs_em = self.maxwell_op.rhs(t, self.em_fields)

        self.field_solve_timer.stop()

        # add current
        e_components = self.maxwell_op.count_subset(
                self.maxwell_op.get_eh_subset()[0:3])

        from hedge.tools import to_obj_array
        j = to_obj_array(self.cloud.reconstruct_j())
        rhs_em[:e_components] -= 1/self.maxwell_op.epsilon*j

        from hedge.tools import make_obj_array
        return make_obj_array([rhs_em, rhs_cloud])

    def upkeep(self):
        for f in self.em_filters:
            self.em_fields = f(self.em_fields)

        self.cloud.upkeep()




# shape bandwidth -------------------------------------------------------------
def guess_shape_bandwidth(cloud, exponent):
    from pyrticle._internal import ShapeFunction
    cloud.reconstructor.set_shape_function(
            ShapeFunction(
                cloud.mesh_data.advisable_particle_radius(),
                cloud.mesh_data.dimensions,
                exponent,
                ))




def optimize_shape_bandwidth(cloud, analytic_rho, exponent, 
        plot_l1_errors=False,
        visualize=False):
    discr = cloud.discretization

    adv_radius = cloud.mesh_data.advisable_particle_radius()
    radii = [adv_radius*2**i 
            for i in numpy.linspace(-4, 4, 50)]

    if visualize:
        from hedge.visualization import SiloVisualizer
        vis = SiloVisualizer(discr)

    from hedge.discretization import integral
    
    def set_radius(r):
        from pyrticle._internal import ShapeFunction
        cloud.reconstructor.set_shape_function(
                ShapeFunction(r, cloud.mesh_data.dimensions, exponent,))

    tried_radii = []
    l1_errors = []

    import sys

    sys.stdout.write("optimizing shape bw (%d attempts): " % len(radii))
    for step, radius in enumerate(radii):
        sys.stdout.write("%d." % step)
        sys.stdout.flush()

        try:
            cloud.set_ignore_core_warnings(True)
            set_radius(radius)
        except RuntimeError, re:
            if "particle mass is zero" in str(re):
                continue
            else:
                raise
        finally:
            cloud.set_ignore_core_warnings(False)

        cloud.derived_quantity_cache.clear()
        try:
            cloud.set_ignore_core_warnings(True)
            rec_rho = cloud.reconstruct_rho()
        except RuntimeError, re:
            if "particle mass is zero" in str(re):
                continue
            else:
                raise
        finally:
            cloud.set_ignore_core_warnings(False)

        tried_radii.append(radius)
        l1_errors.append(integral(discr, numpy.abs(rec_rho-analytic_rho)))

        if visualize:
            visf = vis.make_file("rho-%04d" % step)
            cloud.add_to_vis(vis, visf, time=radius, step=step)
            vis.add_data(visf, [ 
                ("rho", rec_rho), 
                ("anarho", analytic_rho), 
                ],
                time=radius, step=step)

            try:
                cloud.reconstructor.write_grid_quantities
            except AttributeError:
                pass
            else:
                cloud.reconstructor.write_grid_quantities(visf, 
                        ["rho", "usecount"])

            visf.close()

    sys.stdout.write("\n")
    sys.stdout.flush()

    if visualize:
        vis.close()

    if plot_l1_errors:
        from pylab import semilogx, show
        semilogx(tried_radii, l1_errors)
        show()

    from pytools import argmin
    min_idx = argmin(l1_errors)
    min_rad = tried_radii[min_idx]
    min_l1_error = l1_errors[min_idx]

    rel_l1_error = abs(min_l1_error / integral(discr, analytic_rho))
    if rel_l1_error > 0.1:
        from warnings import warn
        warn("Charge density is very poorly resolved (rel L1 error=%g)" % rel_l1_error)

    def is_local_minimum(list, i):
        if i == 0:
            return False
        elif i == len(list)-1:
            return False
        else:
            return list[i] < list[i-1] and list[i] < list[i+1]
        
    local_minima = [idx for idx in range(len(tried_radii)) 
            if is_local_minimum(l1_errors, idx)]

    chosen_idx = max(local_minima)
    chosen_rad = tried_radii[chosen_idx]
    chosen_l1_error = l1_errors[chosen_idx]

    print "radius: guessed optimum=%g, found optimum=%g, chosen=%g" % (
            adv_radius, min_rad, chosen_rad)

    print "radius: optimum l1 error=%g, chosen l1 error=%g" % (
            min_l1_error, chosen_l1_error)

    set_radius(chosen_rad)
    cloud.derived_quantity_cache.clear()




def set_shape_bandwidth(cloud, shape_bw, shape_exp, analytic_rho):
    if shape_bw.startswith("optimize"):
        optimize_shape_bandwidth(cloud, analytic_rho,
                shape_exp, 
                plot_l1_errors="plot" in shape_bw,
                visualize="visualize" in shape_bw,
                )
    elif shape_bw == "guess":
        guess_shape_bandwidth(cloud, shape_exp)
    else:
        from pyrticle._internal import ShapeFunction
        cloud.reconstructor.set_shape_function(
                ShapeFunction(
                    float(setup.shape_bandwidth),
                    cloud.mesh_data.dimensions,
                    shape_exp,
                    ))




# initial condition -----------------------------------------------------------
def compute_initial_condition(pcon, discr, cloud, mean_beta, max_op,
        debug=False, force_zero=False):
    from hedge.operators import WeakPoissonOperator
    from hedge.mesh import TAG_ALL, TAG_NONE
    from hedge.data import ConstantGivenFunction, GivenVolumeInterpolant
    from hedge.discretization import norm

    def rel_l2_error(field, true):
        err = norm(discr, field-true)
        if err == 0:
            return 0
        else:
            return err/norm(discr, true)

    gamma = (1-numpy.dot(mean_beta, mean_beta))**(-0.5)

    # see doc/notes.tm for derivation of IC

    def make_scaling_matrix(beta_scale, other_scale):
        if la.norm(mean_beta) < 1e-10:
            return other_scale*numpy.eye(discr.dimensions) 
        else:
            beta_unit = mean_beta/la.norm(mean_beta)
            return (other_scale*numpy.identity(discr.dimensions) 
                    + (beta_scale-other_scale)*numpy.outer(beta_unit, beta_unit))

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

    from hedge.tools import ptwise_dot
    e_tilde = ptwise_dot(make_scaling_matrix(1/gamma, 1), poisson_op.grad(phi_tilde))
    e_prime = ptwise_dot(make_scaling_matrix(1, gamma), e_tilde)
    h_prime = (1/max_op.mu)*gamma/max_op.c * max_op.e_cross(mean_beta, e_tilde)

    if debug:
        from hedge.discretization import integral
        reconstructed_charge = integral(discr, rho_prime)

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
        cloud.add_to_vis(vis, visf)
        visf.close()

    return FieldsAndCloud(max_op, e_prime, h_prime, cloud)

