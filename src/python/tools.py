import pyrticle._internal as _internal
import pylinear.computation as comp




def v_from_p(p, m, c):
    from math import sqrt
    value =  c*p*(comp.norm_2_squared(p)+c*c*m*m)**(-0.5)
    return value




ZeroVector = _internal.ZeroVector
