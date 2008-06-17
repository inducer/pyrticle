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
import numpy
import numpy.linalg as la
import pytools




# warnings --------------------------------------------------------------------
class WarningForwarder(_internal.WarningListener):
    def note_warning(self, message, filename, lineno):
        from warnings import warn_explicit
        warn_explicit(message, UserWarning, filename, lineno)




class WarningIgnorer(_internal.WarningListener):
    def note_warning(self, message, filename, lineno):
        pass




warning_forwarder = WarningForwarder()




# number-shifting vectors -----------------------------------------------------
class NumberShiftMultiplexer(_internal.NumberShiftListener):
    def __init__(self):
        _internal.NumberShiftListener.__init__(self)
        from weakref import WeakKeyDictionary
        self.subscribers = WeakKeyDictionary()

    def subscribe(self, subscriber):
        assert isinstance(subscriber, _internal.NumberShiftListener)
        self.subscribers[subscriber] = None

    def note_change_size(self, new_size):
        for subscriber in self.subscribers.iterkeys():
            subscriber.note_change_size(new_size)

    def note_move(self, orig, dest, size):
        for subscriber in self.subscribers.iterkeys():
            subscriber.note_move(orig, dest, size)

    def note_reset(self, start, size):
        for subscriber in self.subscribers.iterkeys():
            subscriber.note_reset(start, size)




class NumberShiftableVector(_internal.NumberShiftListener):
    """A vector that may be notified of shifts in the DOFs it contains.

    This solves the following problem: Particles in 
    L{pyrticle.cloud.ParticleCloud} may die at any time, including
    during a Runge-Kutta timestep. The state vector inside Runge-Kutta,
    needs to be notified that degrees of freedom may have shifted and/or
    been deleted. This class, together with the corresponding
    L{NumberShiftSignaller}, fulfills that purpose.
    """

    def __init__(self, vector, signaller):
        _internal.NumberShiftListener.__init__(self)
        self.vector = vector
        self.signaller = signaller
        signaller.subscribe(self)

    @staticmethod
    def unwrap(instance):
        if isinstance(instance, NumberShiftableVector):
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
        return NumberShiftableVector(
                self.vector + self.unwrap(other),
                self.signaller)

    __radd__ = __add__
        
    def __iadd__(self, other):
        self.vector += self.unwrap(other)
        return self

    def __mul__(self, other):
        result = NumberShiftableVector(
                self.vector * self.unwrap(other),
                self.signaller)
        return result

    __rmul__ = __mul__

    # shiftiness --------------------------------------------------------------
    def note_change_size(self, new_size):
        old_size = len(self.vector)

        if new_size > old_size:
            raise NotImplementedError, "vector growing unimplemented"
        elif new_size < old_size:
            self.vector = self.vector[:new_size]

    def note_move(self, orig, dest, size):
        self.vector[dest] = self.vector[orig]

    def note_reset(self, start, size):
        self.vector[start:(start+size)] = 0




# math stuff ------------------------------------------------------------------
def uniform_on_unit_sphere(dim):
    from random import gauss

    # cf.
    # http://www-alg.ist.hokudai.ac.jp/~jan/randsphere.pdf
    # Algorith due to Knuth

    pt = numpy.array([gauss(0,1) for i in range(dim)])
    n2 = la.norm(pt)
    return pt/n2




class ODEDefinedFunction:
    def __init__(self, t0, y0, dt):
        self.t = [t0]
        self.y = [y0]
        self.dt = dt

        from hedge.timestep import RK4TimeStepper
        self.forward_stepper = RK4TimeStepper()
        self.backward_stepper = RK4TimeStepper()

    def __call__(self, t):
        def copy_if_necessary(x):
            try:
                return x[:]
            except TypeError:
                return x

        if t < self.t[0]:
            steps = int((self.t[0]-t)/self.dt)+1
            t_list = [self.t[0]]
            y_list = [self.y[0]]
            for n in range(steps):
                y_list.append(self.backward_stepper(
                    copy_if_necessary(y_list[-1]), 
                    t_list[-1], -self.dt, self.rhs))
                t_list.append(t_list[-1]-self.dt)

            self.t = t_list[:0:-1] + self.t
            self.y = y_list[:0:-1] + self.y
        elif t >= self.t[-1]:
            steps = int((t-self.t[-1])/self.dt)+1
            t_list = [self.t[-1]]
            y_list = [self.y[-1]]
            for n in range(steps):
                y_list.append(self.forward_stepper(
                    copy_if_necessary(y_list[-1]), 
                    t_list[-1], self.dt, self.rhs))
                t_list.append(t_list[-1]+self.dt)

            self.t = self.t + t_list[1:]
            self.y = self.y + y_list[1:]

        from bisect import bisect_right
        below_idx = bisect_right(self.t, t)-1
        assert below_idx >= 0
        above_idx = below_idx + 1

        assert above_idx < len(self.t)
        assert self.t[below_idx] <= t <= self.t[above_idx]

        # FIXME linear interpolation, bad
        slope = ((self.y[above_idx]-self.y[below_idx]) 
                /
                (self.t[above_idx]-self.t[below_idx]))

        return self.y[below_idx] + (t-self.t[below_idx]) * slope

    def rhs(self, t, y):
        raise NotImplementedError
