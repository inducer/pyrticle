import pylinear.array as num
import pylinear.computation as comp
from pytools.arithmetic_container import work_with_arithmetic_containers
from pytools.arithmetic_container import ArithmeticList
from math import sqrt
import pyrticle._internal as _internal




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
        for face in discr.mesh.tag_to_boundary[None]:
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
    - positions and velocities, layout as
      [x0,y0,z0,x1,y1,z1,...]
    - charges, masses
      single scalar per particle
    - containing_elements
      a reference to a hedge.mesh.Element containing the particle
      (or None, indicating this particle index is dead--
      remember, it's called particle-*in*-cell :-)
     """
    def __init__(self, discr, epsilon, mu, verbose_vis=False):
        _internal.ParticleCloud.__init__(self, 
                discr.dimensions,
                len(discr.mesh.points),
                len(discr.mesh.elements),
                len(discr.element_groups))

        self.verbose_vis = verbose_vis
        self.discretization = discr

        self.mesh_info.add_elements(discr)
        self.mesh_info.add_vertices(discr)
        self.mesh_info.add_nodes(len(discr.nodes), discr.nodes)
               
        self.epsilon = epsilon
        self.mu = mu
        self.c = 1/(epsilon*mu)**0.5

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

    def add_particles(self, positions, velocities, charges, masses):
        """Add the particles with the given data to the cloud."""

        dim = self.mesh_info.dimensions

        new_count = len(positions)

        for v in velocities:
            assert comp.norm_2(v) < self.c

        containing_elements = [self.mesh_info.find_containing_element(p) 
                for p in positions]

        try:
            len(charges)
        except:
            charges = charges*num.ones((new_count,))

        try:
            len(masses)
        except:
            masses = new_count * [masses]

        deathflags = [ce == MeshInfo.INVALID_ELEMENT 
                for ce in containing_elements]
        containing_elements = [ce for ce in containing_elements 
                if ce != MeshInfo.INVALID_ELEMENT]
        positions = [p for p, dead in zip(positions, deathflags) 
                if not dead]
        velocities = [v for v, dead in zip(velocities, deathflags) 
                if not dead]
        charges = num.array([c for c, dead in zip(charges, deathflags) 
            if not dead])
        masses = num.array([m for m, dead in zip(masses, deathflags) 
            if not dead])

        for p in positions:
            assert len(p) == dim
        for v in velocities:
            assert len(v) == dim

        # first, fill up the spots of formerly dead particles
        already_placed = 0
        while len(self.deadlist):
            i = self.deadlist.pop()
            pstart = i*dim
            pend = (i+1)*dim

            self.containing_elements[i] = containing_elements[alreay_placed]
            self.positions[pstart:pend] = positions[already_placed]
            self.velocities[pstart:pend] = velocities[already_placed]
            self.charges[i] = charges[already_placed]
            self.masses[i] = masses[already_placed]

            already_placed += 1

        # next, append
        self.containing_elements[-1:] = \
                containing_elements[already_placed:]

        self.positions = num.vstack((self.positions, 
            num.vstack(positions[already_placed:])))
        self.velocities = num.vstack((self.velocities, 
            num.vstack(velocities[already_placed:])))

        self.charges = num.vstack((self.charges, 
            charges[already_placed:]))
        self.masses = num.vstack((self.masses, 
            masses[already_placed:]))

    def kill_particle(self, pn):
        dim = self.mesh_info.dimensions

        pstart = pn*dim
        pend = (pn+1)*dim

        periodicity_trip = False

        pt = self.positions[pstart:pend]
        for i, (axis_interval, xi) in enumerate(zip(self.periodicity, pt)):
            if axis_interval is not None:
                xmin, xmax = axis_interval
                if not (xmin <= xi <= xmax):
                    pt[i] = (xi-xmin) % (xmax-xmin) + xmin;
                    periodicity_trip = True

        if periodicity_trip:
            self.positions[pstart:pend] = pt
            ce = self.find_new_containing_element(
                    pn, self.containing_elements[pn])
            if ce != MeshInfo.INVALID_ELEMENT:
                self.containing_elements[pn] = ce
                self.periodic_hits.tick()
                return

        print "KILL %d" % pn
        
        self.containing_elements[pn] = MeshInfo.INVALID_ELEMENT
        self.deadlist.append(pn)

        self.positions[pstart:pend] = num.zeros((dim,))
        self.velocities[pstart:pend] = num.zeros((dim,))

    def upkeep(self):
        """Perform any operations must fall in between timesteps,
        such as resampling or deleting particles.
        """
        pass

    def reconstruct_densities(self):
        """Return a tuple (charge_density, current_densities), where
        current_densities is an 
          ArithmeticList([[jx0,jx1,...],[jy0,jy1,...]])  
        of the densities in each direction.
        """

        rho = self.discretization.volume_zeros()
        j = ArithmeticList([self.discretization.volume_zeros()
            for axis in range(self.discretization.dimensions)])

        self._reconstruct_densities(rho, j[0], j[1], j[2], self.particle_radius)

        if self.verbose_vis:
            self.vis_info["rho"] = rho
            self.vis_info["j"] = j

        return rho, j

    def rhs(self, t, e, h):
        return ArithmeticList([
            self.velocities, 
            self.accelerations(
                e[0], e[1], e[2], 
                h[0], h[1], h[2], 
                self.c, self.mu,
                self.verbose_vis)
            ])

    def __iadd__(self, rhs):
        from pytools import argmin, argmax

        assert isinstance(rhs, ArithmeticList)

        dx, dv = rhs
        self.positions += dx
        self.velocities += dv

        self.update_containing_elements()

        return self

    def add_to_vis(self, visualizer, vis_file, time, step):
        from hedge.visualization import VtkVisualizer, SiloVisualizer
        if isinstance(visualizer, VtkVisualizer):
            return self.add_to_vtk(visualizer, vis_file, time, step)
        elif isinstance(visualizer, SiloVisualizer):
            return self.add_to_silo(vis_file, time, step)
        else:
            raise ValueError, "unknown visualizer type `%s'" % type(visualizer)

    def add_to_vtk(self, visualizer, vis_file, time, step):
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
                    self.velocities,
                    vector_format=VF_INTERLEAVED,
                    components=dim)
                )

        def add_vis_vector(name):
            if name in self.vis_info:
                vec = self.vis_info[name]
            else:
                vec = num.zeros((len(self.containing_elements) * dim,))

            grid.add_pointdata(
                    DataArray("velocity", 
                        vec, vector_format=VF_INTERLEAVED,
                        components=dim)
                    )

        mesh_scalars = []
        mesh_vectors = []

        if self.verbose_vis:
            add_vis_vector("pt_e")
            add_vis_vector("pt_h")
            add_vis_vector("el_acc")
            add_vis_vector("lorentz_acc")

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

    def add_to_silo(self, db, time, step):
        from pylo import DBOPT_DTIME, DBOPT_CYCLE
        dim = self.mesh_info.dimensions

        coords = num.vstack([self.positions[i::dim] for i in range(dim)])
        db.put_pointmesh("particles", dim, coords, 
                {DBOPT_DTIME: time, DBOPT_CYCLE:step})

        db.put_pointvar("velocity", "particles", 
                [self.velocities[i::dim] for i in range(dim)])

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
            add_vis_vector("el_acc")
            add_vis_vector("lorentz_acc")

            mesh_scalars.append(("rho", self.vis_info["rho"]))
            mesh_vectors.append(("j", self.vis_info["j"]))

        return mesh_scalars, mesh_vectors
