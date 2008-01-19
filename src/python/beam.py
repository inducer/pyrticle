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




def rms_beam_size(cloud, axis):
    from pyrticle._internal import rms_beam_size
    return rms_beam_size(cloud.pic_algorithm, axis)




def rms_emittance(cloud, axis, beam_axis):
    from pyrticle._internal import rms_beam_emittance
    return rms_beam_emittance(cloud.pic_algorithm, axis, beam_axis)




def rms_energy_spread(cloud):
    from pyrticle._internal import rms_energy_spread
    return rms_energy_spread(cloud.pic_algorithm)





