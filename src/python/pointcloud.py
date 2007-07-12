import pylinear.array as num
from pytools.arithmetic_container import work_with_arithmetic_containers




def find_containing_element(discr, point):
    for el in discr.mesh.elements:
        unit_coords = el.contains_point(point)
        if unit_coors:
            return el, unit_coords




class Interpolator:
    def __init__(self, discr, el, x):
        unit_coords = el.inverse_map(x)

        (self.el_start, self.el_end), ldis = discr.find_el_data(el.id)

        point_vdm = num.array([f(unit_coords) for f in ldis.basis_functions()])
        self.interp_coeff = ldis.vandermonde() <<num.leftsolve>> point_vdm

    @work_with_arithmetic_containers
    def __call__(self, field):
            return self.interp_coeff * field[el_start:el_end]




class PointCloud:
    """State container for a cloud of particles. Supports particle
    problems of any dimension, examples below are given for three
    dimensions for simplicity.

    It contains the following data elements:
    - positions and velocities, layout as
      [x0,y0,z0,x1,y1,z1,...]
    - charges, masses
      single scalar per particle
     """
    def __init__(self, discr, epsilon, mu):
        self.discretization = discr
        self.dimensions = discr.dimensions
        self.positions = num.zeros((0,))
        self.velocities = num.zeros((0,))
        self.charges = num.zeros((0,))
        self.masses = num.zeros((0,))
        self.containing_elements = []

        self.epsilon = epsilon
        self.mu = mu

    def __len__(self):
        return len(self.charges)

    def add_points(self, positions, velocities, charges, masses):
        """Adds the points with the given data to the cloud."""

        new_count = len(self.positions)
        self.positions = num.vstack((self.positions, num.vstack(positions)))
        self.velocities = num.vstack((self.velocities, num.vstack(velocities)))

        try:
            len(charges)
        except:
            charges = charges*num.ones((new_count,))

        try:
            len(masses)
        except:
            masses = masses*num.ones((new_count,))

        self.charges = num.vstack((self.charges, num.array(charges)))
        self.masses = num.vstack((self.masses, num.array(masses)))
        self.containing_elements.extend(
                find_containing_element(p) for p in positions)

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
        pass

    def rhs(self, e, h):
        from pytools.arithmetic_container import ArithmeticList

        accelerations = num.zeros(self.velocities.shape)

        d = self.dimensions

        for i in range(len(self)):
            interp = Interpolator(self.discretization,
                    self.containing_element[i],
                    self.positions[d*i:d*i+3])
            e_at_pt = num.array(interp(e))
            h_at_pt = num.array(interp(h))
            v = self.velocities[d*i:d*i+3]
            q = self.charges[i]
            force = q*(e_at_pt + v <<num.cross>> (self.mu*h_at_pt))
            accelerations[d*i:d*i+3] = force/self.masses[i]

        return ArithmeticList([self.velocities, accelerations])

    def __iadd__(self, rhs):
        assert isinstance(rhs, ArithmeticList)
        dx, dv = rhs
        self.positions += dx
        self.velocities += dv
        self.containing_elements = [
                find_containing_element(p) for p in self.positions]

    def add_to_silo_db(self, db):
        d = self.dimensions
        coords = num.vstack([self.positions[i::d] for i in range(self.dimensions)])
        velocities = [self.velocities[i::d] for i in range(self.dimensions)]

        db.put_pointmesh("particles", d, coords)
        db.put_pointvar("velocities", "particles", velocities)

