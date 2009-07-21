print "enabling multirate..."
from hedge.timestep.multirate_ab import \
        TwoRateAdamsBashforthTimeStepper as _TwoRateAB
_enable_multirate = True
if _enable_multirate:
    _substep_count = 20
    timestepper_maker = lambda dt: _TwoRateAB(
            "fastest_first_1a",
            dt, substep_count=_substep_count, order=3)
    _rk4_stability = 2
    #dt_scale = 0.36/_rk4_stability*_substep_count # ab4

    # be conservative for jan-data run
    dt_scale = 0.30/_rk4_stability*_substep_count # ab4
    #dt_scale = 0.18/_rk4_stability*_substep_count # ab5

vis_interval = 10
