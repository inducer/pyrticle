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
from pytools import memoize_method

from pyrticle.meshdata import MeshData




class PicState(object):
    def __init__(self, method,
            particle_count=None,
            containing_elements=None,
            positions=None,
            momenta=None,
            charges=None,
            masses=None,
            depositor_state=None,
            pusher_state=None,
            pnss=None,
            vis_listener=None,
            ):
        state_class = getattr(_internal, "ParticleState%s" % 
            method.get_dimensionality_suffix())

        pstate = self.particle_state = state_class()
        if depositor_state is None:
            self.depositor_state = method.depositor.make_state(self)
        else:
            self.depositor_state = depositor_state

        if pusher_state is None:
            self.pusher_state = method.pusher.make_state(self)
        else:
            self.pusher_state = pusher_state

        if particle_count is None:
            pstate.particle_count = 0
            pstate.containing_elements = numpy.zeros((0,), dtype=numpy.uint32)
            pstate.positions = numpy.zeros((0, pstate.xdim), dtype=float)
            pstate.momenta = numpy.zeros((0, pstate.vdim), dtype=float)
            pstate.charges = numpy.zeros((0,), dtype=float)
            pstate.masses = numpy.zeros((0,), dtype=float)
        else:
            pstate.particle_count = particle_count
            pstate.containing_elements = containing_elements
            pstate.positions = positions
            pstate.momenta = momenta
            pstate.charges = charges
            pstate.masses = masses

        self.derived_quantity_cache = {}

        if pnss is None:
            from pyrticle.tools import StatePassingNumberShiftMultiplexer
            self.particle_number_shift_signaller = \
                    StatePassingNumberShiftMultiplexer(self)
        else:
            self.particle_number_shift_signaller = pnss

        # visualization
        if vis_listener is None:
            from pyrticle.tools import MapStorageVisualizationListener
            self.vis_listener = MapStorageVisualizationListener(
                    self.particle_number_shift_signaller)
        else:
            self.vis_listener = vis_listener


    def __len__(self):
        return self.particle_state.particle_count

    @property
    def positions(self):
        return self.particle_state.positions[:self.particle_state.particle_count]

    @property
    def momenta(self):
        return self.particle_state.momenta[:self.particle_state.particle_count]

    @property
    def masses(self):
        return self.particle_state.masses[:self.particle_state.particle_count]

    @property
    def charges(self):
        return self.particle_state.charges[:self.particle_state.particle_count]

    def resize(self, newsize):
        pstate = self.particle_state

        pstate.containing_elements = numpy.resize(
                pstate.containing_elements, (newsize,))
        pstate.positions = numpy.resize(
                pstate.positions, (newsize, pstate.xdim))
        pstate.momenta = numpy.resize(
                pstate.momenta, (newsize, pstate.vdim))
        pstate.charges = numpy.resize(pstate.charges, (newsize,))
        pstate.masses = numpy.resize(pstate.masses, (newsize,))

    def clear(self):
        pstate = self.particle_state

        pstate.particle_count = 0
        self.depositor_state.clear()
        self.derived_quantity_cache.clear()

    # derived quantity cache --------------------------------------------------
    def get_derived_quantity_from_cache(self, name, getter):
        try:
            return self.derived_quantity_cache[name]
        except KeyError:
            self.derived_quantity_cache[name] = value = getter()
            return value

    def get_derived_quantities_from_cache(self, names, getters, all_getter=None):
        # all "getters" elements should update the cache themselves
        cached = tuple(self.derived_quantity_cache.get(name) for name in names)
        cached_number = sum(1 for v in cached if cached is not None)

        if cached_number == len(names):
            return cached
        elif cached_number == 0 and all_getter:
            all_values = all_getter()
            for name, value in zip(names, all_values):
                self.derived_quantity_cache[name] = value
            return all_values
        else:
            return tuple(c or g() for c, g in zip(cached, getters))




class TimesteppablePicState(object):
    def __init__(self, method, state):
        self.method = method
        self.state = state

    def __add__(self, update):
        from pyrticle.tools import NumberShiftableVector

        from pytools import typedump
        dx, dp, drecon = update
        dx = NumberShiftableVector.unwrap(dx)
        dp = NumberShiftableVector.unwrap(dp)

        new_state = self.method.advance_state(self.state, dx, dp, drecon)

        return TimesteppablePicState(self.method, new_state)




class FieldRhsCalculator(object):
    def __init__(self, method, maxwell_op):
        self.method = method
        self.maxwell_op = maxwell_op

        self.bound_maxwell_op = maxwell_op.bind(self.method.discretization)

        from pytools.log import IntervalTimer
        self.field_solve_timer = IntervalTimer(
                "t_field",
                "Time spent in field solver")

    def add_instrumentation(self, mgr):
        mgr.add_quantity(self.field_solve_timer)

    def __call__(self, t, fields_f, state_f):
        # calculate EM right-hand side 
        sub_timer = self.field_solve_timer.start_sub_timer()

        from pyrticle.hyperbolic import CleaningMaxwellOperator
        if isinstance(self.maxwell_op, CleaningMaxwellOperator):
            rhs_fields = self.bound_maxwell_op(t, fields_f(), 
                    self.method.deposit_rho(state_f()))
        else:
            rhs_fields = self.bound_maxwell_op(t, fields_f())

        sub_timer.stop().submit()

        return rhs_fields



class ParticleToFieldRhsCalculator(object):
    def __init__(self, method, maxwell_op):
        self.method = method
        self.maxwell_op = maxwell_op

    def __call__(self, t, fields_f, state_f):
        return self.maxwell_op.assemble_eh(
                e=-1/self.maxwell_op.epsilon
                *self.method.deposit_j(state_f()))




class FieldToParticleRhsCalculator(object):
    def __init__(self, method, maxwell_op):
        self.method = method
        self.maxwell_op = maxwell_op

    def __call__(self, t, fields_f, state_f):
        from hedge.tools import ZeroVector

        fields = fields_f()
        state = state_f()

        e, h = self.maxwell_op.split_eh(fields)

        # assemble field_args of the form [ex,ey,ez] and [bx,by,bz],
        # inserting ZeroVectors where necessary.
        idx = 0
        e_arg = []
        for use_component in self.maxwell_op.get_eh_subset()[0:3]:
            if use_component:
                e_arg.append(e[idx])
                idx += 1
            else:
                e_arg.append(ZeroVector())

        idx = 0
        b_arg = []
        for use_component in self.maxwell_op.get_eh_subset()[3:6]:
            if use_component:
                b_arg.append(self.maxwell_op.mu * h[idx])
                idx += 1
            else:
                b_arg.append(ZeroVector())

        field_args = tuple(e_arg) + tuple(b_arg)

        velocities = self.method.velocities(state)

        # compute forces
        forces = self.method.pusher.forces(
                state,
                velocities,
                *field_args)

        from pyrticle.tools import NumberShiftableVector
        from hedge.tools import make_obj_array
        result = make_obj_array([
            0,
            NumberShiftableVector(forces, 
                signaller=state.particle_number_shift_signaller),
            0])
        return result




class ParticleRhsCalculator(object):
    def __init__(self, method, maxwell_op):
        self.method = method
        self.maxwell_op = maxwell_op

    def __call__(self, t, fields_f, state_f):
        state = state_f()

        velocities = self.method.velocities(state)

        from pyrticle.tools import NumberShiftableVector
        from hedge.tools import make_obj_array
        result = make_obj_array([
            NumberShiftableVector(velocities, 
                signaller=state.particle_number_shift_signaller),
            0,
            self.method.depositor.rhs(state)
            ])
        return result




class PicMethod(object):
    """
    @arg debug: A set of strings telling what to debug. So far, the
      following debug flags are in use:
    
      - depositor: Debug the depositor.
      - verbose_vis: Generate E and B fields and force 
        visualizations at particle locations.
      - ic: Check the initial condition when it's generated.
      - no_ic: Start with zero fields.
      - discretization: (See driver.py.) Turn on debug mode for the
        discretization.
      - shape_bw: Debug the finding of the optimal shape bandwidth.

      - interactive: Allow debug measures that require user interaction.
      - vis_files: Allow debug measures that write extra visualization
        files.
    """

    def __init__(self, discr, units,
            depositor, pusher,
            dimensions_pos, dimensions_velocity,
            debug=set()):

        self.units = units
        self.discretization = discr
        self.debug = debug

        self.depositor = depositor
        self.pusher = pusher

        self.dimensions_mesh = discr.dimensions
        self.dimensions_pos = dimensions_pos
        self.dimensions_velocity = dimensions_velocity

        dims = (dimensions_pos, dimensions_velocity)

        self.mesh_data = _internal.MeshData(discr.dimensions)
        self.mesh_data.fill_from_hedge(discr)

        # subsystem init
        self.depositor.initialize(self)
        self.pusher.initialize(self)

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


    def make_state(self):
        state = PicState(self)

        state.particle_number_shift_signaller.subscribe_with_state(self.depositor)
        state.particle_number_shift_signaller.subscribe_with_state(self.pusher)

        return state

    def get_dimensionality_suffix(self):
        return "%dd%dv" % (self.dimensions_pos, self.dimensions_velocity)

    def get_shape_function_class(self):
        from pyrticle.tools import \
                PolynomialShapeFunction, \
                CInfinityShapeFunction

        from pyrticle._internal import get_shape_function_name
        name = get_shape_function_name()

        if name == "c_infinity":
            return CInfinityShapeFunction
        elif name == "polynomial":
            return PolynomialShapeFunction
        else:
            raise ValueError, "unknown shape function class"

    def set_ignore_core_warnings(self, ignore):
        from pyrticle.tools import WarningForwarder, WarningIgnorer
        import pyrticle.tools
        del pyrticle.tools.warning_forwarder
        if ignore:
            pyrticle.tools.warning_forwarder = WarningIgnorer()
        else:
            pyrticle.tools.warning_forwarder = WarningForwarder()

    def add_instrumentation(self, mgr, observer):
        mgr.add_quantity(self.find_el_timer)
        mgr.add_quantity(self.find_same_counter)
        mgr.add_quantity(self.find_by_neighbor_counter)
        mgr.add_quantity(self.find_by_vertex_counter)
        mgr.add_quantity(self.find_global_counter)

        self.depositor.add_instrumentation(mgr, observer)
        self.pusher.add_instrumentation(mgr, observer)

    def velocities(self, state):
        try:
            return state.derived_quantity_cache["velocities"]
        except KeyError:
            result = _internal.get_velocities(
                    state.particle_state, self.units.VACUUM_LIGHT_SPEED())
            state.derived_quantity_cache["velocities"] = result
            return result

    def mean_beta(self, state):
        if len(state):
            return numpy.average(self.velocities(state), axis=0) \
                    / self.units.VACUUM_LIGHT_SPEED()
        else:
            return numpy.zeros((self.dimensions_velocity,))

    def add_particles(self, state, iterable, maxcount=None):
        """Add the  particles from C{iterable} to the cloud.
        
        C{iterable} is expected to yield tuples 
        C{(position, velocity, charge, mass)}.
        
        If C{maxcount} is specified, maximally C{maxcount}
        particles are obtained from the iterable.
        """

        pstate = state.particle_state

        if maxcount is not None:
            if pstate.particle_count+maxcount >= len(pstate.containing_elements):
                state.resize(pstate.particle_count+maxcount)

        for pos, vel, charge, mass in iterable:
            if maxcount is not None:
                if maxcount == 0:
                    break
                maxcount -= 1

            assert len(pos) == self.dimensions_pos
            assert len(vel) == self.dimensions_velocity

            pos = numpy.asarray(pos)
            vel = numpy.asarray(vel)
            mom = mass*self.units.gamma_from_v(vel)*vel 

            cont_el = self.mesh_data.find_containing_element(pos) 
            if cont_el == MeshData.INVALID_ELEMENT:
                print "not in valid element"

                continue

            if pstate.particle_count >= len(pstate.containing_elements):
                state.resize(max(128, 2*pstate.particle_count))

            pstate.containing_elements[pstate.particle_count] = cont_el
            pstate.positions[pstate.particle_count] = pos
            pstate.momenta[pstate.particle_count] = mom
            pstate.charges[pstate.particle_count] = charge
            pstate.masses[pstate.particle_count] = mass

            pstate.particle_count += 1

        self.check_containment(state)
        state.particle_number_shift_signaller.note_change_size(
                pstate.particle_count)
        state.derived_quantity_cache.clear()

    def check_containment(self, state):
        """Check that a containing element is known for each particle.

        This is a new invariant as of 1/17/08, and violations of this end up
        segfaulting, which we should avoid.
        """
        assert (state.particle_state.containing_elements[:len(state)] 
                != MeshData.INVALID_ELEMENT).all()
                
    def upkeep(self, state):
        """Perform any operations must fall in between timesteps,
        such as resampling or deleting particles.
        """
        self.depositor.upkeep(state)
        self.pusher.upkeep(state)

        state.vis_listener.clear()



    # deposition ----------------------------------------------------------
    def deposit_densities(self, state):
        """Return a tuple (charge_density, current_densities), where
        current_densities is an d-by-n array, where d is the number 
        of velocity dimensions, and n is the discretization nodes.
        """
        def all_getter():
            rho, j = self.depositor.deposit_densities(state, self.velocities())
            j = numpy.asarray(j.T, order="C")
            return rho, j

        return state.get_derived_quantities_from_cache(
                ["rho", "j"],
                [self.deposit_rho, self.deposit_j],
                all_getter)

    def deposit_j(self, state):
        """Return a the current densities as an d-by-n array, where d 
        is the number of velocity dimensions, and n is the number of 
        discretization nodes.
        """

        def j_getter():
            return numpy.asarray(
                    self.depositor.deposit_j(
                        state, self.velocities(state))
                    .T, order="C")

        return state.get_derived_quantity_from_cache("j", j_getter)

    def deposit_rho(self, state):
        """Return a the charge_density as a volume vector."""

        def rho_getter():
            return self.depositor.deposit_rho(state)

        return state.get_derived_quantity_from_cache("rho", rho_getter)




    # time advance ------------------------------------------------------------
    def advance_state(self, state, dx, dp, drecon):
        pstate = state.particle_state
        cnt = pstate.particle_count

        positions = pstate.positions.copy()
        momenta = pstate.momenta.copy()
        positions[:cnt] += dx
        momenta[:cnt] += dp

        new_state = PicState(
                self,
                particle_count=cnt,
                containing_elements=pstate.containing_elements,
                positions=positions,
                momenta=momenta,
                charges=pstate.charges,
                masses=pstate.masses,
                depositor_state=self.depositor.advance_state(
                    state, drecon),
                pusher_state=self.pusher.advance_state(state),
                pnss=state.particle_number_shift_signaller,
                vis_listener=state.vis_listener,
                )

        from pyrticle._internal import FindEventCounters
        find_counters = FindEventCounters()

        class BHitListener(_internal.BoundaryHitListener):
            def note_boundary_hit(subself, pn):
                _internal.kill_particle(
                        new_state.particle_state, 
                        pn, new_state.particle_number_shift_signaller)

        from pyrticle._internal import update_containing_elements
        sub_timer = self.find_el_timer.start_sub_timer()
        update_containing_elements(
                self.mesh_data, new_state.particle_state, 
                BHitListener(), find_counters)
        sub_timer.stop().submit()

        self.find_same_counter.transfer(
                find_counters.find_same)
        self.find_by_neighbor_counter.transfer(
                find_counters.find_by_neighbor)
        self.find_by_vertex_counter.transfer(
                find_counters.find_by_vertex)
        self.find_global_counter.transfer(
                find_counters.find_global)

        return new_state

    # visualization -----------------------------------------------------------
    def get_mesh_vis_vars(self):
        return self.vis_listener.mesh_vis_map.items()

    def add_to_vis(self, visualizer, vis_file, state, time=None, step=None, beamaxis=None,
            vis_listener=None):
        from hedge.visualization import VtkVisualizer, SiloVisualizer
        if isinstance(visualizer, VtkVisualizer):
            return self._add_to_vtk(visualizer, vis_file, state, time, step)
        elif isinstance(visualizer, SiloVisualizer):
            return self._add_to_silo(visualizer, vis_file, state, time, step, beamaxis, vis_listener)
        else:
            raise ValueError, "unknown visualizer type `%s'" % type(visualizer)

    def _add_to_silo(self, visualizer, db, state, time, step, beamaxis, vis_listener=None):
        from pylo import DBOPT_DTIME, DBOPT_CYCLE
        from warnings import warn

        optlist = {}
        if time is not None:
            optlist[DBOPT_DTIME] = time
        if step is not None:
            optlist[DBOPT_CYCLE] = step

        pcount = len(state)

        if pcount:
            # real-space ------------------------------------------------------
            db.put_pointmesh("particles", 
                    numpy.asarray(state.positions.T, order="C"), optlist)
            db.put_pointvar1("charge", "particles", state.charges)
            db.put_pointvar1("mass", "particles", state.masses)
            db.put_pointvar("momentum", "particles", 
                    numpy.asarray(state.momenta.T, order="C"))
            db.put_pointvar("velocity", "particles", 
                    numpy.asarray(self.velocities(state).T, order="C"))

            if vis_listener is not None:
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








# shape bandwidth -------------------------------------------------------------
def guess_shape_bandwidth(method, state, exponent):
    method.depositor.set_shape_function(state,
            method.get_shape_function_class()(
                method.mesh_data.advisable_particle_radius(),
                method.mesh_data.dimensions,
                exponent,
                ))




def optimize_shape_bandwidth(method, state, analytic_rho, exponent):
    discr = method.discretization
    rec = method.depositor

    adv_radius = method.mesh_data.advisable_particle_radius()
    radii = [adv_radius*2**i 
            for i in numpy.linspace(-4, 2, 50)]

    def set_radius(r):
        method.depositor.set_shape_function(
                state,
                method.get_shape_function_class()
                (r, method.mesh_data.dimensions, exponent,))

    tried_radii = []
    l1_errors = []

    debug = "shape_bw" in method.debug
    visualize = set(["shape_bw", "vis_files"]) <= method.debug

    if visualize:
        from hedge.visualization import SiloVisualizer
        vis = SiloVisualizer(discr)

    import sys

    if debug:
        sys.stdout.write("optimizing shape bw (%d attempts): " % len(radii))
    for step, radius in enumerate(radii):
        if debug:
            sys.stdout.write("%d." % step)
            sys.stdout.flush()

        try:
            try:
                method.set_ignore_core_warnings(True)
                set_radius(radius)
            except RuntimeError, re:
                if "particle mass is zero" in str(re):
                    continue
                else:
                    raise
        finally:
            method.set_ignore_core_warnings(False)

        state.derived_quantity_cache.clear()

        try:
            try:
                method.set_ignore_core_warnings(True)
                rec_rho = method.deposit_rho(state)
            except RuntimeError, re:
                if "particle mass is zero" in str(re):
                    continue
                else:
                    raise
        finally:
            method.set_ignore_core_warnings(False)

        tried_radii.append(radius)
        l1_errors.append(discr.integral(numpy.abs(rec_rho-analytic_rho)))

        if visualize:
            visf = vis.make_file("bwopt-%04d" % step)
            method.add_to_vis(vis, visf, state, time=radius, step=step)
            vis.add_data(visf, [ 
                ("rho", rec_rho), 
                ("j", method.deposit_j(state)),
                ("rho_analytic", analytic_rho), 
                ],
                expressions=[("rho_diff", "rho-rho_analytic")],
                time=radius, step=step)

            try:
                rec.visualize_grid_quantities
            except AttributeError:
                pass
            else:
                rec.visualize_grid_quantities(visf, [
                        ("rho_grid", rec.deposit_grid_rho()),
                        ("j_grid", rec.deposit_grid_j(method.velocities(state))),
                        ("rho_resid", rec.remap_residual(rec.deposit_grid_rho())),
                        ])

            visf.close()

    if debug:
        sys.stdout.write("\n")
        sys.stdout.flush()

    if visualize:
        vis.close()

    from pytools import argmin
    min_idx = argmin(l1_errors)
    min_rad = tried_radii[min_idx]
    min_l1_error = l1_errors[min_idx]

    rel_l1_error = abs(min_l1_error / discr.integral(analytic_rho))
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

    if "shape_bw" in method.debug:
        chosen_l1_error = l1_errors[chosen_idx]
        print "radius: guessed optimum=%g, found optimum=%g, chosen=%g" % (
                adv_radius, min_rad, chosen_rad)
        print "radius: optimum l1 error=%g, chosen l1 error=%g" % (
                min_l1_error, chosen_l1_error)

    set_radius(chosen_rad)
    state.derived_quantity_cache.clear()

    if set(["interactive", "shape_bw"]) < method.debug:
        from pylab import semilogx, show
        semilogx(tried_radii, l1_errors)
        show()




# initial condition -----------------------------------------------------------
def compute_initial_condition(rcon, discr, method, state,
        maxwell_op, potential_bc,
        force_zero=False):
    from hedge.models.poisson import WeakPoissonOperator
    from hedge.mesh import TAG_ALL, TAG_NONE
    from hedge.data import ConstantGivenFunction, GivenVolumeInterpolant

    def rel_l2_error(field, true):
        from hedge.tools import relative_error
        return relative_error(
                discr.norm(field-true),
                discr.norm(true))

    mean_beta = method.mean_beta(state)
    gamma = method.units.gamma_from_beta(mean_beta)

    # see doc/notes.tm for derivation of IC

    def make_scaling_matrix(beta_scale, other_scale):
        if la.norm(mean_beta) < 1e-10:
            return other_scale*numpy.eye(discr.dimensions) 
        else:
            beta_unit = mean_beta/la.norm(mean_beta)
            return (other_scale*numpy.identity(discr.dimensions) 
                    + (beta_scale-other_scale)*numpy.outer(beta_unit, beta_unit))

    poisson_op = WeakPoissonOperator(
            discr.dimensions,
            diffusion_tensor=ConstantGivenFunction(
                make_scaling_matrix(1/gamma**2, 1)),
            dirichlet_tag=TAG_ALL,
            neumann_tag=TAG_NONE,
            dirichlet_bc=potential_bc)

    rho_prime = method.deposit_rho(state) 
    rho_tilde = rho_prime/gamma

    bound_poisson = poisson_op.bind(discr)
    if force_zero:
        phi_tilde = discr.volume_zeros()
    else:
        from hedge.tools import parallel_cg
        phi_tilde = -parallel_cg(rcon, -bound_poisson, 
                bound_poisson.prepare_rhs(
                    GivenVolumeInterpolant(discr, rho_tilde/maxwell_op.epsilon)), 
                debug=40 if "poisson" in method.debug else False, tol=1e-10)

    from hedge.tools import ptwise_dot
    from hedge.models.nd_calculus import GradientOperator
    #e_tilde = ptwise_dot(2, 1, make_scaling_matrix(1/gamma, 1), bound_poisson.grad(phi_tilde))

    from pyrticle.tools import make_cross_product
    v_e_to_h_cross = make_cross_product(method, maxwell_op, "v", "e", "h")

    e_tilde = ptwise_dot(2, 1, make_scaling_matrix(1/gamma, 1), 
            GradientOperator(discr.dimensions).bind(discr)(phi_tilde))
    e_prime = ptwise_dot(2, 1, make_scaling_matrix(1, gamma), e_tilde)
    h_prime = (1/maxwell_op.mu)*gamma/maxwell_op.c * v_e_to_h_cross(mean_beta, e_tilde)

    if "ic" in method.debug:
        deposited_charge = discr.integral(rho_prime)

        real_charge = numpy.sum(state.charges)
        print "charge: supposed=%g deposited=%g error=%g %%" % (
                real_charge,
                deposited_charge,
                100*abs(deposited_charge-real_charge)/abs(real_charge)
                )

        from hedge.models.nd_calculus import DivergenceOperator

        bound_div_op = DivergenceOperator(discr.dimensions).bind(discr)

        d_tilde = maxwell_op.epsilon*e_tilde
        d_prime = maxwell_op.epsilon*e_prime

        from hedge.optemplate import InverseMassOperator
        divD_prime_ldg = bound_poisson.div(d_prime)
        divD_prime_ldg2 = bound_poisson.div(d_prime, maxwell_op.epsilon*gamma*phi_tilde)
        divD_prime_ldg3 = maxwell_op.epsilon*\
                (InverseMassOperator().apply(discr, bound_poisson.op(gamma*phi_tilde)))
        divD_prime_central = bound_div_op(d_prime)

        print "l2 div D_prime error central: %g" % \
                rel_l2_error(divD_prime_central, rho_prime)
        print "l2 div D_prime error ldg: %g" % \
                rel_l2_error(divD_prime_ldg, rho_prime)
        print "l2 div D_prime error ldg with phi: %g" % \
                rel_l2_error(divD_prime_ldg2, rho_prime)
        print "l2 div D_prime error ldg with phi 3: %g" % \
                rel_l2_error(divD_prime_ldg3, rho_prime)

        if "vis_files" in method.debug:
            from hedge.visualization import SiloVisualizer
            vis = SiloVisualizer(discr)
            visf = vis.make_file("ic")
            vis.add_data(visf, [ 
                ("phi_moving", phi_tilde),
                ("rho_moving", rho_tilde), 
                ("e_moving", e_tilde), 

                ("rho_lab", rho_prime), 
                ("divD_lab_ldg", divD_prime_ldg),
                ("divD_lab_ldg2", divD_prime_ldg2),
                ("divD_lab_ldg3", divD_prime_ldg3),
                ("divD_lab_central", divD_prime_central),
                ("e_lab", e_prime), 
                ("h_lab", h_prime), 
                ],
                )
            method.add_to_vis(vis, visf, state)
            visf.close()

    return maxwell_op.assemble_eh(e=e_prime, h=h_prime)

