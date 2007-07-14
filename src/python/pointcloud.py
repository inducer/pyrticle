import pylinear.array as num
import pylinear.computation as comp
from pytools.arithmetic_container import work_with_arithmetic_containers
from pytools.arithmetic_container import ArithmeticList
from math import sqrt




def find_containing_element(discr, point):
    for el in discr.mesh.elements:
        if el.contains_point(point):
            return el
    return None




def enum_subvectors(x, subdim):
    for i in range(len(x)//subdim):
        yield x[i*subdim:(i+1)*subdim]




class Interpolator:
    def __init__(self, discr, el, x):
        unit_coords = el.inverse_map(x)

        (self.el_start, self.el_end), ldis = discr.find_el_data(el.id)

        point_vdm = num.array([f(unit_coords) for f in ldis.basis_functions()])
        self.interp_coeff = ldis.vandermonde() <<num.leftsolve>> point_vdm

    @work_with_arithmetic_containers
    def __call__(self, field):
            return self.interp_coeff * field[self.el_start:self.el_end]




class PointCloud:
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
        self.discretization = discr
        self.dimensions = discr.dimensions

        self.containing_elements = []
        self.positions = num.zeros((0,))
        self.velocities = num.zeros((0,))
        self.charges = num.zeros((0,))
        self.masses = num.zeros((0,))

        self.deadlist = []

        self.epsilon = epsilon
        self.mu = mu
        self.c = 1/sqrt(epsilon*mu)

        self.vis_info = {}

        self._build_vertex_to_element_map()
        self._build_neighbor_map()

    def _build_vertex_to_element_map(self):
        self.vertex_to_element_map = {}
        for el in self.discretization.mesh.elements:
            for vi in el.vertex_indices:
                self.vertex_to_element_map \
                        .setdefault(vi, []).append(el)

    def _build_neighbor_map(self):
        self.neighbor_map = {}
        for face, (e2, f2) in self.discretization.mesh.both_interfaces():
            self.neighbor_map[face] = e2

    def __len__(self):
        return len(self.charges) - len(self.deadlist)

    def add_points(self, positions, velocities, charges, masses):
        """Adds the points with the given data to the cloud."""

        dim = self.dimensions

        new_count = len(positions)

        for v in velocities:
            assert comp.norm_2(v) < self.c

        containing_elements = [find_containing_element(self.discretization, p) 
                for p in positions]

        try:
            len(charges)
        except:
            charges = charges*num.ones((new_count,))

        try:
            len(masses)
        except:
            masses = new_count * [masses]

        deathflags = [ce is None for ce in containing_elements]
        containing_elements = [ce for ce in containing_elements 
                if ce is not None]
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
        self.containing_elements.extend(
                containing_elements[already_placed:])

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

        dim = self.dimensions

        self.vis_info = {}

        for i in range(len(self.charges)):
            if self.containing_elements[i] is None:
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

            self.vis_info[i, "pt_e"] = e_at_pt
            self.vis_info[i, "pt_h"] = h_at_pt
            self.vis_info[i, "el_force"] = el_force
            self.vis_info[i, "lorentz_force"] = lorentz_force

            accelerations[pstart:pend] = force/rel_mass

        return ArithmeticList([self.velocities, accelerations])

    def __iadd__(self, rhs):
        from pytools import argmin, argmax

        assert isinstance(rhs, ArithmeticList)

        dx, dv = rhs
        self.positions += dx
        self.velocities += dv

        discr = self.discretization

        def find_new_containing_element(i, p, v, prev_el):
            if prev_el:
                if prev_el.contains_point(p):
                    return prev_el
                else:
                    # look via normal -----------------------------------------
                    best_normal_fi = argmax(v*n for n in prev_el.face_normals)
                    el_candidate = self.neighbor_map[prev_el, best_normal_fi]
                    if el_candidate.contains_point(p):
                        return el_candidate
                   
                    # look via closest vertex ---------------------------------
                    closest_vi = argmin(comp.norm_2(discr.mesh.points[vi]-p)
                            for vi in prev_el.vertex_indices)

                    el_candidates = self.vertex_to_element_map[
                            prev_el.vertex_indices[closest_vi]]

                    for el in el_candidates:
                        if el.contains_point(p):
                            return el

                    # look globally -------------------------------------------
                    result = find_containing_element(discr, p)
                    if not result:
                        self.deadlist.append(i)
                    return result
            else:
                return None

        self.containing_elements = [
                find_new_containing_element(i, p, v, prev_el) 
                for i, (p, v, prev_el) in enumerate(zip(
                        enum_subvectors(self.positions, self.dimensions),
                        enum_subvectors(self.velocities, self.dimensions),
                        self.containing_elements))]

        return self

    def add_to_silo_fast(self, db, e, h):
        d = self.dimensions
        coords = num.vstack([self.positions[i::d] for i in range(self.dimensions)])
        velocities = [self.velocities[i::d] for i in range(self.dimensions)]
        
        db.put_pointmesh("particles", d, coords)
        db.put_pointvar("velocities", "particles", velocities)

    def add_to_silo(self, db, e, h):
        dim = self.dimensions

        def make_empty_vis_field():
            return [[] for d in range(dim)]

        coords = make_empty_vis_field()
        pt_v = make_empty_vis_field()

        def add_vector_field(name, var):
            db.put_pointvar(name, "particles", 
                    [num.array(vi) for vi in var])

        for i, ce in enumerate(self.containing_elements):
            if ce is None:
                continue

            pstart = dim*i
            pend = dim*(i+1)
           
            for d in range(dim):
                coords[d].append(self.positions[pstart + d])
                pt_v[d].append(self.velocities[pstart + d])

        from operator import add
        db.put_pointmesh("particles", dim, reduce(add, coords))
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

        add_vis_info_vector("pt_e")
        add_vis_info_vector("pt_h")
        add_vis_info_vector("el_force")
        add_vis_info_vector("lorentz_force")
