import random as _random
_random.seed(0)

pusher = PushMonomial()
#pusher = PushAverage()
if True:
    depositor = DepGrid(
            el_tolerance=0.1,
            submethod="simplex_reduce",
            #submethod="simplex_extra",
            #submethod="simplex_enlarge",
            #submethod="brick",
            #jiggle_radius=0.001,
            enforce_continuity=True
            )
#depositor = DepAdv()
#depositor = DepShape()
#depositor = DepNormShape()
#depositor = DepGridFind()

debug.add("shape_bw")
#debug.add("no_ic")
#debug.add("ic")

dimensions_pos = 2
dimensions_velocity = 2

beam_axis = 0
beam_diag_axis = 1
tube_length = 2

# chi = 4

shape_bandwidth = "optimize"
shape_bandwidth = 0.1

_cloud_charge = 10e-9 * units.C
#nparticles = 1
element_order = 4
final_time = 10*units.M/units.VACUUM_LIGHT_SPEED()
_electrons_per_particle = abs(_cloud_charge/nparticles/units.EL_CHARGE)

shape_exponent = 2

_tube_width = 1
import hedge.mesh.generator as _meshgen
mesh = _meshgen.make_rect_mesh(
        a=(-0.5, -_tube_width/2),
        b=(-0.5+tube_length, _tube_width/2),
        periodicity=(True, False),
        subdivisions=(10,5),
        max_area=0.02)

_c0 = units.VACUUM_LIGHT_SPEED()

_mean_v = numpy.array([_c0*0.9,0])
_sigma_v = numpy.array([_c0*0.9*1e-3, 0.1*_c0])

_mean_beta = _mean_v/_c0
_gamma = units.gamma_from_v(_mean_v)
_pmass = _electrons_per_particle*units.EL_MASS
_mean_p = _gamma*_pmass*_mean_v

distribution = pyrticle.distribution.JointParticleDistribution([
    pyrticle.distribution.GaussianPos([0.0,0.0], [0.01, 0.01]),
    pyrticle.distribution.GaussianMomentum(
        #_mean_p, _sigma_v*_gamma*_pmass, 
        _mean_p, _sigma_v*_gamma*_pmass, 
        units,
        pyrticle.distribution.DeltaChargeMass(
            _cloud_charge/nparticles,
            _pmass))
    ])

vis_interval = 1
#vis_order = 8

if isinstance(depositor, DepGrid):
    def hook_visualize(runner, vis, visf, observer):
        meth = runner.method
        dep = meth.depositor
        state = observer.state

        dep.visualize_grid_quantities(visf, [
                ("rho_grid", dep.deposit_grid_rho(state)),
                ("j_grid", dep.deposit_grid_j(state, meth.velocities(state))),
                ("ones_resid", dep.remap_residual(dep.ones_on_grid())),
                ("rho_resid", dep.remap_residual(dep.deposit_grid_rho(state))),
                ("usecount", dep.grid_usecount()),
                ])

def hook_vis_quantities(observer):
    return [
                    ("e", observer.e), 
                    ("h", observer.h), 
                    ("j", observer.method.deposit_j(observer.state)), 
                    ("rho", observer.method.deposit_rho(observer.state)), 
                    #("phi", observer.phi), 
                    ]
