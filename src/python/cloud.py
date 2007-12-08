import pylinear.array as num
import pylinear.computation as comp
from pytools.arithmetic_container import work_with_arithmetic_containers
from pytools.arithmetic_container import ArithmeticList
from math import sqrt
import pyrticle._internal as _internal




ZeroVector = _internal.ZeroVector




def v_from_p(p, m, c):
    from math import sqrt
    value =  c*p*(comp.norm_2_squared(p)+c*c*m*m)**(-0.5)
    return value

def find_containing_element(discr, point):
    for el in discr.mesh.elements:
        if el.contains_point(point):
            return el
    return None




class MeshInfo(_internal.MeshInfo):
    pass




def add_mesh_info_methods():
    def add_local_discretizations(self, discr):
        from hedge.polynomial import generic_vandermonde

        ldis_indices = {}
        
        for i, eg in enumerate(discr.element_groups):
            ldis = eg.local_discretization
            ldis_indices[ldis] = i

            mon_basis = [_internal.MonomialBasisFunction(*idx)
                    for idx in ldis.node_tuples()]
            mon_vdm = generic_vandermonde(ldis.unit_nodes(), mon_basis).T

            l_vdm, u_vdm, perm, sign = comp.lu(mon_vdm)
            p_vdm = num.permutation_matrix(from_indices=perm)

            self.add_local_discretization(
                    mon_basis, l_vdm, u_vdm, p_vdm)

        return ldis_indices

    def add_elements(self, discr):
        ldis_indices = self.add_local_discretizations(discr)

        mesh = discr.mesh

        neighbor_map = {}
        for face, (e2, f2) in discr.mesh.both_interfaces():
            neighbor_map[face] = e2.id
        from hedge.mesh import TAG_ALL
        for face in discr.mesh.tag_to_boundary[TAG_ALL]:
            neighbor_map[face] = MeshInfo.INVALID_ELEMENT

        for el in mesh.elements:
            (estart, eend), ldis = discr.find_el_data(el.id)

            self.add_element(
                    el.inverse_map,
                    ldis_indices[ldis],
                    estart, eend,
                    [vi for vi in el.vertex_indices],
                    el.face_normals,
                    set(neighbor_map[el,fi] for fi in xrange(len(el.faces)))
                    )

    def add_vertices(self, discr):
        mesh = discr.mesh

        vertex_to_element_map = {}
        for el in mesh.elements:
            for vi in el.vertex_indices:
                vertex_to_element_map.setdefault(vi, set()).add(el.id)
                for other_vi in mesh.periodic_opposite_vertices.get(vi, []):
                    vertex_to_element_map.setdefault(other_vi, set()).add(el.id)

        for vi in xrange(len(mesh.points)):
            self.add_vertex(vi, 
                    mesh.points[vi],
                    vertex_to_element_map[vi])

    _internal.MeshInfo.add_local_discretizations = add_local_discretizations
    _internal.MeshInfo.add_elements = add_elements
    _internal.MeshInfo.add_vertices = add_vertices
add_mesh_info_methods()




class ParticleCloud:
    """State container for a cloud of particles. Supports particle
    problems of any dimension, examples below are given for three
    dimensions for simplicity.
    """
    def __init__(self, discr, units, 
            dimensions_pos, dimensions_velocity,
            verbose_vis=False):

        self.units = units
        self.discretization = discr
        self.verbose_vis = verbose_vis

        self.dimensions_mesh = discr.dimensions
        self.dimensions_pos = dimensions_pos
        self.dimensions_velocity = dimensions_velocity

        dims = (dimensions_pos, dimensions_velocity)

        constructor_args  = (
                self.dimensions_mesh,
                len(discr.mesh.points), len(discr.mesh.elements), len(discr.element_groups),
                units.VACUUM_LIGHT_SPEED)

        if dims == (3,3):
            self.icloud = _internal.ParticleCloud33(*constructor_args)
        elif dims == (2,2):
            self.icloud = _internal.ParticleCloud22(*constructor_args)
        else:
            raise ValueError, "unsupported combination of dimensions"

        self.icloud.mesh_info.add_elements(discr)
        self.icloud.mesh_info.add_vertices(discr)
        self.icloud.mesh_info.add_nodes(len(discr.nodes), discr.nodes)

        def min_vertex_distance(el):
            vertices = [discr.mesh.points[vi] 
                    for vi in el.vertex_indices]

            return min(min(comp.norm_2(vi-vj)
                    for i, vi in enumerate(vertices)
                    if i != j)
                    for j, vj in enumerate(vertices))

        self.particle_radius = 0.5*min(min_vertex_distance(el) 
                for el in discr.mesh.elements)

        for axis, ((ax_min, ax_max), periodicity_tags) in enumerate(zip(
                zip(*discr.mesh.bounding_box), discr.mesh.periodicity)):
            if periodicity_tags is not None:
                self.icloud.mesh_info.add_periodicity(
                        axis, ax_min, ax_max)

    def __len__(self):
        return len(self.icloud.containing_elements) - len(self.icloud.deadlist)

    @property
    def positions(self):
        return self.icloud.positions

    @property
    def momenta(self):
        return self.icloud.momenta

    @property
    def masses(self):
        return self.icloud.masses

    @property
    def charges(self):
        return self.icloud.charges

    def velocities(self):
        return self.icloud.velocities()

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
        containing_elements = [self.icloud.mesh_info.find_containing_element(p) 
                for p in positions]

        # weed out uncontained particles
        deathflags = [ce == MeshInfo.INVALID_ELEMENT 
                for ce in containing_elements]
        containing_elements = [ce for ce in containing_elements 
                if ce != MeshInfo.INVALID_ELEMENT]
        positions = [p for p, dead in zip(positions, deathflags) 
                if not dead]
        momenta = [v for v, dead in zip(momenta, deathflags) 
                if not dead]
        charges = num.array([c for c, dead in zip(charges, deathflags) 
            if not dead])
        masses = num.array([m for m, dead in zip(masses, deathflags) 
            if not dead])

        # check vector dimensionalities
        for x in positions:
            assert len(x) == self.dimensions_pos
        for p in momenta:
            assert len(p) == self.dimensions_velocity

        # first, fill up the spots of formerly dead particles
        already_placed = 0
        while len(self.icloud.deadlist):
            i = self.icloud.deadlist.pop()
            x_pstart = i*self.dimensions_pos
            x_pend = (i+1)*self.dimensions_pos
            v_pstart = i*self.dimensions_velocity
            v_pend = (i+1)*self.dimensions_velocity

            self.icloud.containing_elements[i] = containing_elements[alreay_placed]
            self.icloud.positions[x_pstart:x_pend] = positions[already_placed]
            self.icloud.momenta[v_pstart:v_pend] = momenta[already_placed]
            self.icloud.charges[i] = charges[already_placed]
            self.icloud.masses[i] = masses[already_placed]

            already_placed += 1

        # next, append
        self.icloud.containing_elements[-1:] = \
                containing_elements[already_placed:]

        self.icloud.positions = num.hstack((self.icloud.positions, 
            num.hstack(positions[already_placed:])))
        self.icloud.momenta = num.hstack((self.icloud.momenta, 
            num.hstack(momenta[already_placed:])))

        self.icloud.charges = num.hstack((self.icloud.charges, 
            charges[already_placed:]))
        self.icloud.masses = num.hstack((self.icloud.masses, 
            masses[already_placed:]))

    def upkeep(self):
        """Perform any operations must fall in between timesteps,
        such as resampling or deleting particles.
        """
        self.icloud.vis_info.clear()

    def reconstruct_densities(self, velocities=None):
        """Return a tuple (charge_density, current_densities), where
        current_densities is an 
          ArithmeticList([[jx0,jx1,...],[jy0,jy1,...]])  
        of the densities in each direction.
        """

        if velocities is None:
            velocities = self.icloud.velocities()

        rho = self.discretization.volume_zeros()
        j = ArithmeticList([self.discretization.volume_zeros()
            for axis in range(self.dimensions_velocity)])

        self.icloud.reconstruct_densities(rho, j, self.particle_radius,
                velocities)

        if self.verbose_vis:
            self.icloud.vis_info["rho"] = rho
            self.icloud.vis_info["j"] = j

        return rho, j

    def reconstruct_j(self, velocities=None):
        """Return a the current densities in the structure::
          ArithmeticList([[jx0,jx1,...],[jy0,jy1,...]])  
        """

        if velocities is None:
            velocities = self.icloud.velocities()

        j = ArithmeticList([self.discretization.volume_zeros()
            for axis in range(self.dimensions_velocity)])

        self.icloud.reconstruct_j(j, self.particle_radius, velocities)

        if self.verbose_vis:
            self.icloud.vis_info["j"] = j

        return j

    def reconstruct_rho(self):
        """Return a the charge_density as a volume vector.
        """

        rho = self.discretization.volume_zeros()

        self.icloud.reconstruct_rho(rho, self.particle_radius)

        if self.verbose_vis:
            self.icloud.vis_info["rho"] = rho

        return rho

    def rhs(self, t, e, b, velocities=None):
        """Return an ArithmeticList of velocities and forces on the particles.

        @arg e: triple of M{E_x}, M{E_y}, M{E_z}, each of which may be either 
          a Pylinear vector or a L{ZeroVector}.
        @arg b: triple of M{B_x}, M{B_y}, M{B_z}, each of which may be either 
          a Pylinear vector or a L{ZeroVector}. Caution: The hedge Maxwell operator
          deals with M{H}, not M{B}.
        """

        if velocities is None:
            velocities = self.icloud.velocities()

        field_args = tuple(e) + tuple(b)
        # compute forces
        return ArithmeticList([
            velocities, 
            self.icloud.forces(
                velocities=velocities,
                verbose_vis=self.verbose_vis,
                *field_args
                )
            ])

    def __iadd__(self, rhs):
        from pytools import argmin, argmax

        assert isinstance(rhs, ArithmeticList)

        dx, dp = rhs
        self.icloud.positions += dx
        self.icloud.momenta += dp

        self.icloud.update_containing_elements()

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

        dim = self.icloud.mesh_info.dimensions

        points = (len(self.icloud.containing_elements),
                DataArray("points", self.icloud.positions,
                    vector_format=VF_INTERLEAVED,
                    components=dim))
        grid = UnstructuredGrid(
                points, 
                cells=[[i] for i in range(points[0])],
                cell_types=[VTK_VERTEX] * points[0],
                )

        grid.add_pointdata(
                DataArray("velocity", 
                    self.icloud.velocities(),
                    vector_format=VF_INTERLEAVED,
                    components=dim)
                )

        def add_vis_vector(name):
            if name in self.icloud.vis_info:
                vec = self.icloud.vis_info[name]
            else:
                vec = num.zeros((len(self.icloud.containing_elements) * dim,))

            grid.add_pointdata(
                    DataArray(name, 
                        vec, vector_format=VF_INTERLEAVED,
                        components=dim)
                    )

        mesh_scalars = []
        mesh_vectors = []

        if verbose:
            add_vis_vector("pt_e")
            add_vis_vector("pt_b")
            add_vis_vector("el_force")
            add_vis_vector("lorentz_force")

            mesh_scalars.append(("rho", self.vis_info["rho"]))
            mesh_vectors.append(("j", self.vis_info["j"]))

        from os.path import splitext
        pathname = splitext(vis_file.pathname)[0] + "-particles.vtu"
        from pytools import assert_not_a_file
        assert_not_a_file(pathname)

        outf = open(pathname, "w")
        AppendedDataXMLGenerator()(grid).write(outf)
        outf.close()

        visualizer.register_pathname(time, pathname)

        return mesh_scalars, mesh_vectors

    def _add_to_silo(self, db, time, step, verbose):
        from pylo import DBOPT_DTIME, DBOPT_CYCLE
        from warnings import warn

        optlist = {}
        if time is not None:
            optlist[DBOPT_DTIME] = time
        if step is not None:
            optlist[DBOPT_CYCLE] = step

        coords = num.hstack([
            self.icloud.positions[i::self.dimensions_pos] 
            for i in range(self.dimensions_pos)])
        db.put_pointmesh("particles", self.dimensions_pos, coords, optlist)

        db.put_pointvar("momenta", "particles", 
                [self.icloud.momenta[i::self.dimensions_velocity] 
                    for i in range(self.dimensions_velocity)])

        velocities = self.icloud.velocities()
        db.put_pointvar("velocity", "particles", 
                [velocities[i::self.dimensions_velocity] 
                    for i in range(self.dimensions_velocity)])

        pcount = len(self.icloud.containing_elements)
        def add_particle_vis(name, dim):
            if name in self.icloud.vis_info:
                db.put_pointvar(name, "particles", 
                        [self.icloud.vis_info[name][i::dim] for i in range(dim)])
            else:
                warn("writing zero for particle visualization variable '%s'" % name)
                db.put_pointvar(name, "particles", 
                        [num.zeros((pcount,)) for i in range(dim)])

        def add_mesh_vis(name, unavail_ok=False):
            if name in self.icloud.vis_info:
                mesh_scalars.append((name, self.icloud.vis_info[name]))
            else:
                if not unavail_ok:
                    warn("visualization of unavailable mesh variable '%s' requested" % name)

        mesh_scalars = []
        mesh_vectors = []

        if verbose:
            add_particle_vis("pt_e", 3)
            add_particle_vis("pt_b", 3)
            add_particle_vis("el_force", 3)
            add_particle_vis("lorentz_force", 3)

            add_mesh_vis("rho", unavail_ok=True)
            add_mesh_vis("j")

        return mesh_scalars, mesh_vectors




class FieldsAndCloud:
    def __init__(self, maxwell_op, e, h, cloud):
        from pytools.arithmetic_container import join_fields
        self.maxwell_op = maxwell_op
        self.e = e
        self.h = h
        self.cloud = cloud

        self.eh_components = maxwell_op.component_count()

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

        velocities = self.cloud.velocities()

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

        cloud_rhs = self.cloud.rhs(t, e_arg, b_arg, velocities)

        rhs_e, rhs_h = self.maxwell_op.split_fields(
                self.maxwell_op.rhs(t, join_fields(self.e, self.h))
                )

        return join_fields(
                rhs_e - 1/self.maxwell_op.epsilon*self.cloud.reconstruct_j(velocities),
                rhs_h,
                ).plus([cloud_rhs])




def compute_initial_condition(pcon, discr, cloud, mean_beta, max_op,
        debug=False):
    from hedge.operators import WeakPoissonOperator
    from hedge.mesh import TAG_ALL, TAG_NONE
    from hedge.data import ConstantGivenFunction, GivenVolumeInterpolant
    from hedge.tools import dot

    def l2_norm(field):
        return sqrt(dot(field, discr.mass_operator*field))
    def l2_error(field, true):
        return l2_norm(field-true)/l2_norm(true)

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

        divD_tilde_ldg = poisson_op.div(d_tilde)
        divD_tilde_ldg2 = poisson_op.div(d_tilde, max_op.epsilon*phi_tilde)
        divD_tilde_central = div_op(d_tilde)

        divD_prime_ldg = poisson_op.div(d_prime)
        divD_prime_ldg2 = poisson_op.div(d_prime, max_op.epsilon*gamma*phi_tilde)
        divD_prime_ldg3 = max_op.epsilon*\
                (discr.inverse_mass_operator*poisson_op.op(gamma*phi_tilde))
        divD_prime_central = div_op(d_prime)

        print "l2 div D_tilde error central: %g" % \
                l2_error(divD_tilde_central, rho_tilde)
        print "l2 div D_tilde error ldg: %g" % \
                l2_error(divD_tilde_ldg, rho_tilde)
        print "l2 div D_tilde error ldg2: %g" % \
                l2_error(divD_tilde_ldg2, rho_tilde)

        print "l2 div D_prime error central: %g" % \
                l2_error(divD_prime_central, rho_prime)
        print "l2 div D_prime error ldg: %g" % \
                l2_error(divD_prime_ldg, rho_prime)
        print "l2 div D_prime error ldg with phi: %g" % \
                l2_error(divD_prime_ldg2, rho_prime)
        print "l2 div D_prime error ldg with phi 3: %g" % \
                l2_error(divD_prime_ldg3, rho_prime)

        from hedge.visualization import SiloVisualizer
        vis = SiloVisualizer(discr)
        visf = vis.make_file("ic")
        vis.add_data(visf, [ 
            ("phi_tilde", phi_tilde),
            ("rho_tilde", rho_tilde), 
            ("divD_tilde_central", divD_tilde_central),
            ("divD_tilde_ldg", divD_tilde_ldg),
            ("divD_tilde_ldg2", divD_tilde_ldg2),
            ("e_tilde", e_tilde), 

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

