"""Little bits of usefulness for Pyrticle"""

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



import pyrticle._internal as _internal
import pylinear.array as num
import pylinear.computation as comp




def v_from_p(p, m, c):
    from math import sqrt
    value =  c*p*(comp.norm_2_squared(p)+c*c*m*m)**(-0.5)
    return value




ZeroVector = _internal.ZeroVector




class DOFShiftSignaller:
    def __init__(self):
        from weakref import WeakKeyDictionary
        self.subscribers = WeakKeyDictionary()

    def subscribe(self, shiftable):
        self.subscribers[shiftable] = None

    def change_size(self, new_size):
        for subscriber in self.subscribers.iterkeys():
            subscriber.change_size(new_size)

    def move_dof(self, orig, dest):
        for subscriber in self.subscribers.iterkeys():
            subscriber.move_dof(orig, dest)



class DOFShiftForwarder(_internal.DOFShiftListener, DOFShiftSignaller):
    def __init__(self):
        _internal.DOFShiftListener.__init__(self)
        DOFShiftSignaller.__init__(self)

    def note_change_size(self, new_size):
        self.change_size(new_size)

    def note_move_dof(self, orig, dest):
        self.move_dof(orig, dest)




class DOFShiftableVector(object):
    """A vector that may be notified of shifts in the DOFs it contains.

    This solves the following problem: Particles in 
    L{pyrticle.cloud.ParticleCloud} may die at any time, including
    during a Runge-Kutta timestep. The state vector inside Runge-Kutta,
    needs to be notified that degrees of freedom may have shifted and/or
    been deleted. This class, together with the corresponding
    L{DOFShiftSignaller}, fulfills that purpose.
    """

    def __init__(self, vector, signaller):
        self.vector = vector
        self.signaller = signaller
        signaller.subscribe(self)

    @staticmethod
    def unwrap(instance):
        if isinstance(instance, DOFShiftableVector):
            return instance.vector
        else:
            return instance

    # arithmetic --------------------------------------------------------------
    def __add__(self, other):
        if len(self.vector) != len(self.unwrap(other)):
            print len(self.vector), len(self.unwrap(other)), type(other)
            print self.signaller, other.signaller
            print other in other.signaller.subscribers
            print "---------------------------------------------"
        return DOFShiftableVector(
                self.vector + self.unwrap(other),
                self.signaller)

    __radd__ = __add__
        
    def __iadd__(self, other):
        self.vector += self.unwrap(other)
        return self

    def __mul__(self, other):
        result = DOFShiftableVector(
                self.vector * self.unwrap(other),
                self.signaller)
        return result

    __rmul__ = __mul__

    # shiftiness --------------------------------------------------------------
    def change_size(self, new_size):
        old_size = len(self.vector)

        if new_size > old_size:
            self.vector = num.hstack((
                    self.vector, 
                    num.zeros((new_size-len(self.vector),), 
                        dtype=self.vector.dtype)))
        elif new_size < old_size:
            self.vector = self.vector[:new_size]

    def move_dof(self, orig, dest):
        self.vector[dest] = self.vector[orig]
