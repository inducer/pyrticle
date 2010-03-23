import random as _random
_random.seed(0)

dimensions_pos = 2
dimensions_velocity = 2

beam_axis = 0
beam_diag_axis = 1
tube_length = 40*units.MM

chi = None

#shape_bandwidth = "optimize,visualize,plot"
shape_bandwidth = "optimize"

pusher = PushMonomial()
if False:
    depositor = RecGrid(
            #FineCoreBrickGenerator(core_axis=0, core_fraction=0.08),
            el_tolerance=0.1,
            method="simplex_reduce",
            #filter_min_amplification=0.1,
            #filter_order=6,
            )
else:
    depositor = DepShape()

_cloud_charge = -10e-9 * units.C
nparticles = 1
element_order = 3
final_time = 0.1*units.M/units.VACUUM_LIGHT_SPEED()
_electrons_per_particle = abs(_cloud_charge/nparticles/units.EL_CHARGE)

_el_energy = units.EL_REST_ENERGY()*10
_gamma = _el_energy/units.EL_REST_ENERGY()
_mean_beta = (1-1/_gamma**2)**0.5

_tube_width = 33*units.MM
mesh = pyrticle.geometry.make_fine_center_rect_mesh(
        a=(-5*units.MM, -_tube_width/2),
        b=(-5*units.MM+tube_length, _tube_width/2),
        periodicity=(True, False),
        subdivisions=(10,5), inner_subdivisions=(16,2),
        inner_width=2.5*units.MM,
        max_area=1e-5, 
        inner_max_area=0.15*1e-5
        )

distribution = pyrticle.distribution.KVZIntervalBeam(
        units, total_charge=_cloud_charge, 
        p_charge=_cloud_charge/nparticles, 
        p_mass=_electrons_per_particle*units.EL_MASS,
        radii=[1.7*units.MM],
        emittances=[5*units.MM*units.MRAD], 
        z_length=5*units.MM,
        z_pos=10*units.MM,
        beta=_mean_beta,
        axis_first=True)

vis_interval = 10

def hook_vis_quantities(runner):
    return [
            ("e", runner.fields.e), 
            ("h", runner.fields.h), 
            ("j", runner.cloud.deposit_j()), 
            ("rho", runner.cloud.deposit_rho()), 
            ]

if isinstance(depositor, DepGrid):
    def hook_visualize(runner, vis, visf):
        dep = runner.cloud.depositor
        dep.visualize_grid_quantities(visf, [
                ("rho_grid", dep.deposit_grid_rho()),
                ("j_grid", dep.deposit_grid_j(runner.cloud.velocities())),
                ("ones_resid", dep.remap_residual(dep.ones_on_grid())),
                ("rho_resid", dep.remap_residual(dep.deposit_grid_rho())),
                ("usecount", dep.grid_usecount()),
                ])
