import pylinear.array as num
import pylinear.computation as comp
from pytools.arithmetic_container import work_with_arithmetic_containers
from pytools.arithmetic_container import ArithmeticList
from math import sqrt
import pyrticle._internal as _internal




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
                    [neighbor_map[el,fi] for fi in xrange(len(el.faces))]
                    )

    def add_vertices(self, discr):
        mesh = discr.mesh

        vertex_to_element_map = {}
        for el in mesh.elements:
            for vi in el.vertex_indices:
                vertex_to_element_map.setdefault(vi, []).append(el.id)

        for vi in xrange(len(mesh.points)):
            self.add_vertex(vi, 
                    mesh.points[vi],
                    vertex_to_element_map[vi])

    _internal.MeshInfo.add_local_discretizations = add_local_discretizations
    _internal.MeshInfo.add_elements = add_elements
    _internal.MeshInfo.add_vertices = add_vertices
add_mesh_info_methods()




class ParticleCloud(_internal.ParticleCloud):
    """State container for a cloud of particles. Supports particle
    problems of any dimension, examples below are given for three
    dimensions for simplicity.

    It contains the following data elements:
    @ivar positions: layout as [x0,y0,z0,x1,y1,z1,...]
    @ivar momenta: layout as [x0,y0,z0,x1,y1,z1,...]
    @ivar charges: single scalar per particle
    @ivar masses: single scalar per particle
    @ivar containing_elements: a number of a mesh element
      (or MeshInfo.INVALID_ELEMENT, indicating this particle index is dead--
      remember, it's called particle-*in*-cell :-)
     """
    def __init__(self, maxwell_op, 
            dimensions_pos, dimensions_velocity,
            verbose_vis=False):
        discr = maxwell_op.discr

        _internal.ParticleCloud.__init__(self, 
                discr.dimensions, dimensions_pos, dimensions_velocity,
                len(discr.mesh.points), len(discr.mesh.elements),
                len(discr.element_groups),
                maxwell_op.epsilon, maxwell_op.mu)

        self.verbose_vis = verbose_vis
        self.discretization = discr
        self.maxwell_op = maxwell_op

        self.mesh_info.add_elements(discr)
        self.mesh_info.add_vertices(discr)
        self.mesh_info.add_nodes(len(discr.nodes), discr.nodes)
               
        def min_vertex_distance(el):
            vertices = [discr.mesh.points[vi] 
                    for vi in el.vertex_indices]

            return min(min(comp.norm_2(vi-vj)
                    for i, vi in enumerate(vertices)
                    if i != j)
                    for j, vj in enumerate(vertices))

        self.particle_radius = 0.5*min(min_vertex_distance(el) 
                for el in discr.mesh.elements)

        self.periodicity = []
        for axis_interval, periodicity_tags in zip(
                zip(*discr.mesh.bounding_box), discr.mesh.periodicity):
            if periodicity_tags is None:
                self.periodicity.append(None)
            else:
                self.periodicity.append(axis_interval)

    def __len__(self):
        return len(self.containing_elements) - len(self.deadlist)

    def gamma(self, v):
        value = (1-comp.norm_2_squared(v)/self.c**2)**(-0.5)
        if value < 0:
            raise RuntimeError, "particle velocity > speed of light"
        return value

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
        momenta = [m*self.gamma(v)*v for m, v in zip(masses, velocities)]

        # find containing elements
        containing_elements = [self.mesh_info.find_containing_element(p) 
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
        while len(self.deadlist):
            i = self.deadlist.pop()
            x_pstart = i*self.dimensions_pos
            x_pend = (i+1)*self.dimensions_pos
            v_pstart = i*self.dimensions_velocity
            v_pend = (i+1)*self.dimensions_velocity

            self.containing_elements[i] = containing_elements[alreay_placed]
            self.positions[x_pstart:x_pend] = positions[already_placed]
            self.momenta[v_pstart:v_pend] = momenta[already_placed]
            self.charges[i] = charges[already_placed]
            self.masses[i] = masses[already_placed]

            already_placed += 1

        # next, append
        self.containing_elements[-1:] = \
                containing_elements[already_placed:]

        self.positions = num.hstack((self.positions, 
            num.hstack(positions[already_placed:])))
        self.momenta = num.hstack((self.momenta, 
            num.hstack(momenta[already_placed:])))

        self.charges = num.hstack((self.charges, 
            charges[already_placed:]))
        self.masses = num.hstack((self.masses, 
            masses[already_placed:]))

    def kill_particle(self, pn):
        dim = self.mesh_info.dimensions

        x_pstart = pn*self.dimensions_pos
        x_pend = (pn+1)*self.dimensions_pos

        periodicity_trip = False

        pt = self.positions[x_pstart:x_pend]
        for i, (axis_interval, xi) in enumerate(zip(self.periodicity, pt)):
            if axis_interval is not None:
                xmin, xmax = axis_interval
                if not (xmin <= xi <= xmax):
                    pt[i] = (xi-xmin) % (xmax-xmin) + xmin;
                    periodicity_trip = True

        if periodicity_trip:
            self.positions[x_pstart:x_pend] = pt
            ce = self.find_new_containing_element(
                    pn, self.containing_elements[pn])
            if ce != MeshInfo.INVALID_ELEMENT:
                self.containing_elements[pn] = ce
                self.periodic_hits.tick()
                return

        print "KILL %d" % pn
        
        self.containing_elements[pn] = MeshInfo.INVALID_ELEMENT
        self.deadlist.append(pn)

        v_pstart = pn*self.dimensions_velocity
        v_pend = (pn+1)*self.dimensions_velocity

        self.positions[x_pstart:x_pend] = num.zeros((dim,))
        self.momenta[v_pstart:v_pend] = num.zeros((dim,))

    def upkeep(self):
        """Perform any operations must fall in between timesteps,
        such as resampling or deleting particles.
        """
        pass

    def reconstruct_densities(self, velocities=None):
        """Return a tuple (charge_density, current_densities), where
        current_densities is an 
          ArithmeticList([[jx0,jx1,...],[jy0,jy1,...]])  
        of the densities in each direction.
        """

        rho = self.discretization.volume_zeros()
        j = ArithmeticList([self.discretization.volume_zeros()
            for axis in range(self.dimensions_velocity)])

        if velocities is None:
            velocities = self.velocities()

        self._reconstruct_densities(rho, j, self.particle_radius, velocities)

        if self.verbose_vis:
            self.vis_info["rho"] = rho
            self.vis_info["j"] = j

        return rho, j

    def reconstruct_rho(self):
        """Return a the charge_density as a volume vector.
        """

        rho = self.discretization.volume_zeros()

        self._reconstruct_rho(rho, self.particle_radius)

        if self.verbose_vis:
            self.vis_info["rho"] = rho

        return rho

    def rhs(self, t, e, h, velocities=None):
        from pyrticle._internal import ZeroVector
        
        if velocities is None:
            velocities = self.velocities()

        # assemble field_args of the form [ex,ey,ez,hx,hy,hz],
        # inserting ZeroVectors where necessary.
        field_args = []
        idx = 0
        for use_component in self.maxwell_op.get_subset()[0:3]:
            if use_component:
                field_args.append(e[idx])
                idx += 1
            else:
                field_args.append(ZeroVector())
        idx = 0
        for use_component in self.maxwell_op.get_subset()[3:6]:
            if use_component:
                field_args.append(h[idx])
                idx += 1
            else:
                field_args.append(ZeroVector())

        # compute forces
        return ArithmeticList([
            velocities, 
            self.forces(
                velocities=velocities,
                verbose_vis=self.verbose_vis,
                *field_args
                )
            ])

    def __iadd__(self, rhs):
        from pytools import argmin, argmax

        assert isinstance(rhs, ArithmeticList)

        dx, dp = rhs
        self.positions += dx
        self.momenta += dp

        self.update_containing_elements()

        return self

    def add_to_vis(self, visualizer, vis_file, time=None, step=None):
        from hedge.visualization import VtkVisualizer, SiloVisualizer
        if isinstance(visualizer, VtkVisualizer):
            return self._add_to_vtk(visualizer, vis_file, time, step)
        elif isinstance(visualizer, SiloVisualizer):
            return self._add_to_silo(vis_file, time, step)
        else:
            raise ValueError, "unknown visualizer type `%s'" % type(visualizer)

    def _add_to_vtk(self, visualizer, vis_file, time, step):
        from hedge.vtk import \
                VTK_VERTEX, VF_INTERLEAVED, \
                DataArray, \
                UnstructuredGrid, \
                AppendedDataXMLGenerator

        dim = self.mesh_info.dimensions

        points = (len(self.containing_elements),
                DataArray("points", self.positions,
                    vector_format=VF_INTERLEAVED,
                    components=dim))
        grid = UnstructuredGrid(
                points, 
                cells=[[i] for i in range(points[0])],
                cell_types=[VTK_VERTEX] * points[0],
                )

        grid.add_pointdata(
                DataArray("velocity", 
                    self.velocities(),
                    vector_format=VF_INTERLEAVED,
                    components=dim)
                )

        def add_vis_vector(name):
            if name in self.vis_info:
                vec = self.vis_info[name]
            else:
                vec = num.zeros((len(self.containing_elements) * dim,))

            grid.add_pointdata(
                    DataArray(name, 
                        vec, vector_format=VF_INTERLEAVED,
                        components=dim)
                    )

        mesh_scalars = []
        mesh_vectors = []

        if self.verbose_vis:
            add_vis_vector("pt_e")
            add_vis_vector("pt_h")
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

    def _add_to_silo(self, db, time, step):
        from pylo import DBOPT_DTIME, DBOPT_CYCLE
        dim = self.mesh_info.dimensions

        optlist = {}
        if time is not None:
            optlist[DBOPT_DTIME] = time
        if step is not None:
            optlist[DBOPT_CYCLE] = step

        coords = num.hstack([self.positions[i::dim] for i in range(dim)])
        db.put_pointmesh("particles", dim, coords, optlist)

        db.put_pointvar("momenta", "particles", 
                [self.momenta[i::dim] for i in range(dim)])
        velocities = self.velocities()
        db.put_pointvar("velocity", "particles", 
                [velocities[i::dim] for i in range(dim)])

        pcount = len(self.containing_elements)
        def add_vis_vector(name):
            if name in self.vis_info:
                db.put_pointvar(name, "particles", 
                        [self.vis_info[name][i::dim] for i in range(dim)])
            else:
                db.put_pointvar(name, "particles", 
                        [num.zeros((pcount,)) for i in range(dim)])

        mesh_scalars = []
        mesh_vectors = []

        if self.verbose_vis:
            add_vis_vector("pt_e")
            add_vis_vector("pt_h")
            add_vis_vector("el_force")
            add_vis_vector("lorentz_force")

            if "rho" in self.vis_info:
                mesh_scalars.append(("rho", self.vis_info["rho"]))
            if "j" in self.vis_info:
                mesh_vectors.append(("j", self.vis_info["j"]))

        return mesh_scalars, mesh_vectors
