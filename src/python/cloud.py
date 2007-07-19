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




def enum_subvectors(x, subdim):
    for i in range(len(x)//subdim):
        yield x[i*subdim:(i+1)*subdim]




class Interpolator:
    def __init__(self, discr, el_id, x):
        unit_coords = discr.mesh.elements[el_id].inverse_map(x)

        (self.el_start, self.el_end), ldis = discr.find_el_data(el_id)

        point_vdm = num.array([f(unit_coords) for f in ldis.basis_functions()])
        self.interp_coeff = ldis.vandermonde() <<num.leftsolve>> point_vdm

    @work_with_arithmetic_containers
    def __call__(self, field):
            return self.interp_coeff * field[self.el_start:self.el_end]




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
            mon_vdm = generic_vandermonde( ldis.unit_nodes(), mon_basis)

            l_vdm, u_vdm, perm, sign = comp.lu(mon_vdm.T)
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
                    [neighbor_map[el,fi] for fi in range(len(el.faces))]
                    )

    def add_vertices(self, discr):
        mesh = discr.mesh

        vertex_to_element_map = {}
        for el in mesh.elements:
            for vi in el.vertex_indices:
                vertex_to_element_map.setdefault(vi, []).append(el.id)

        for vi in range(len(mesh.points)):
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
    def __init__(self, discr, epsilon, mu):
        _internal.ParticleCloud.__init__(self, 
                discr.dimensions,
                len(discr.mesh.points),
                len(discr.mesh.elements),
                len(discr.element_groups))

        self.discretization = discr

        self.mesh_info.add_elements(discr)
        self.mesh_info.add_vertices(discr)
               
        self.epsilon = epsilon
        self.mu = mu
        self.c = 1/sqrt(epsilon*mu)

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

        return (self.discretization.volume_zeros(),
                ArithmeticList([self.discretization.volume_zeros() 
                    for i in range(self.discretization.dimensions)]))

    def rhs(self, t, e, h):
        accelerations = num.zeros(self.velocities.shape)

        dim = self.mesh_info.dimensions

        #self.vis_info = {}

        for i in range(len(self.charges)):
            if self.containing_elements[i] == MeshInfo.INVALID_ELEMENT:
                continue

            pstart = dim*i
            pend = dim*(i+1)
           
            interp = Interpolator(self.discretization,
                    self.containing_elements[i],
                    self.positions[pstart:pend])
            e_at_pt = num.array(interp(e))
            h_at_pt = num.array(interp(h))

            v = self.velocities[pstart:pend]
            v_scalar = comp.norm_2(v)
            q = self.charges[i]

            el_force = q*e_at_pt
            lorentz_force = q*(v <<num.cross>> (self.mu*h_at_pt))
            force = el_force + lorentz_force

            rel_mass = self.masses[i] / sqrt(1-(v_scalar/self.c)**2)

            #self.vis_info[i, "pt_e"] = e_at_pt
            #self.vis_info[i, "pt_h"] = h_at_pt
            #self.vis_info[i, "el_force"] = el_force
            #self.vis_info[i, "lorentz_force"] = lorentz_force

            accelerations[pstart:pend] = force/rel_mass

        return ArithmeticList([self.velocities, accelerations])

    def __iadd__(self, rhs):
        from pytools import argmin, argmax

        assert isinstance(rhs, ArithmeticList)

        dx, dv = rhs
        self.positions += dx
        self.velocities += dv

        self.update_containing_elements()

        return self

    def add_to_silo_fast(self, db, e, h):
        d = self.dimensions
        coords = num.vstack([self.positions[i::d] for i in range(self.dimensions)])
        velocities = [self.velocities[i::d] for i in range(self.dimensions)]
        
        db.put_pointmesh("particles", d, coords)
        db.put_pointvar("velocities", "particles", velocities)

    def add_to_silo(self, db, e, h):
        dim = self.mesh_info.dimensions

        def make_empty_vis_field():
            return [[] for d in range(dim)]

        coords = make_empty_vis_field()
        pt_v = make_empty_vis_field()

        def add_vector_field(name, var):
            db.put_pointvar(name, "particles", 
                    [num.array(vi) for vi in var])

        for i, ce in enumerate(self.containing_elements):
            if ce == MeshInfo.INVALID_ELEMENT:
                continue

            pstart = dim*i
            pend = dim*(i+1)
           
            for d in range(dim):
                coords[d].append(self.positions[pstart + d])
                pt_v[d].append(self.velocities[pstart + d])

        from operator import add
        coords = reduce(add, coords)

        if not coords:
            return

        db.put_pointmesh("particles", dim, coords)
        add_vector_field("velocities", pt_v)

        def add_vis_info_vector(name):
            vf = make_empty_vis_field()

            try:
                for i, ce in enumerate(self.containing_elements):
                    if ce is None:
                        continue

                    for d in range(dim):
                        vf[d].append(self.vis_info[i, name][d])
            except KeyError:
                vf = [num.zeros((len(self),)) for d in range(dim)]

            db.put_pointvar(name, "particles", 
                    [num.array(vfi) for vfi in vf])

        #add_vis_info_vector("pt_e")
        #add_vis_info_vector("pt_h")
        #add_vis_info_vector("el_force")
        #add_vis_info_vector("lorentz_force")
