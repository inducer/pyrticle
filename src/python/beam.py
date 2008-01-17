"""Beam physics functions"""

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




def calculate_rms_beam_size(cloud, axis):
    from pytools import average
    from math import sqrt

    xdim = cloud.dimensions_pos

    return sqrt(average(x**2 for x in cloud.pic_algorithm.positions[axis::xdim]))




def calculate_rms_emittance(cloud, axis, beam_axis):
    from pytools import average
    from math import sqrt

    xdim = cloud.dimensions_pos
    vdim = cloud.dimensions_velocity

    def calc_xprime(px, pz):
        if pz:
            return px/pz
        else:
            return 0

    xprime = [calc_xprime(px, pz) for px, pz in zip(
        cloud.pic_algorithm.momenta[axis::vdim], 
        cloud.pic_algorithm.momenta[beam_axis::vdim])]
    mean_x_squared = average(x**2 for x in cloud.pic_algorithm.positions[axis::xdim])
    mean_p_squared = average(xp**2 for xp in xprime)
    squared_mean_xp = average(
            x*xp 
            for x, xp in zip(cloud.pic_algorithm.positions[axis::xdim], xprime)
            )**2

    return sqrt(mean_x_squared*mean_p_squared - squared_mean_xp)




def calculate_rms_energy_spread(cloud):
    from pytools import average
    from math import sqrt

    xdim = cloud.dimensions_pos
    vdim = cloud.dimensions_velocity

    velocities = cloud.velocities()

    from pytools import average

    gammas = [cloud.units.gamma(v) for v in cloud.velocities()]
    energies = [(gamma-1)*m*cloud.units.VACUUM_LIGHT_SPEED**2
            for gamma, m in zip(gammas, cloud.masses)]
    mean_energy = average(energies)
    squared_mean_energy = average(energies)**2

    return sqrt(average(energy**2 for energy in energies)
            -squared_mean_energy)/mean_energy





