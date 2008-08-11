Running your first :mod:`pyrticle` Simulation
---------------------------------------------

To run your first :mod:`pyrticle` simulation, change to the
:file:`examples/pic/` directory. To keep things simple, we will run a 2D
simulation of a Gaussian beam in a beam tube.  A parameter file for this
simulation comes in the above-mentioned directory :mod:`pyrticle` distribution,
it is called :file:`simple.py`.

To first give you an idea of what is going to happen, run this simulation now. 
Type the following command into your shell::

    $ python driver.py simple.cpy

You will see :mod:`pyrticle` perform a number of preparatory steps, and then, once
the screen starts looking like this::

    step=0 | t_sim=0 | W_field=2.94401e-06 | t_step=0.932571 | t_eta=0 | n_part=2000
    step=1 | t_sim=4.02369e-11 | W_field=2.94297e-06 | t_step=0.934013 | t_eta=820.506 | n_part=2000
    step=2 | t_sim=8.04738e-11 | W_field=2.94239e-06 | t_step=0.569611 | t_eta=643.318 | n_part=2000
    step=3 | t_sim=1.20711e-10 | W_field=2.94206e-06 | t_step=0.66167 | t_eta=611.747 | n_part=2000
    step=4 | t_sim=1.60948e-10 | W_field=2.94188e-06 | t_step=0.704858 | t_eta=604.52 | n_part=2000
    step=5 | t_sim=2.01185e-10 | W_field=2.94188e-06 | t_step=0.693117 | t_eta=614.649 | n_part=2000
    step=6 | t_sim=2.41422e-10 | W_field=2.94205e-06 | t_step=0.860489 | t_eta=618.254 | n_part=2000
    ...

:mod:`pyrticle` is performing the actual timestepping. :file:`simple.cpy` tells
:mod:`pyrticle` to write out a visualization file every 10 timesteps, so after
a hundred or so timesteps (which should take less than a minute on a modern
computer), you will have a few files to look at. We suggest that you download
`VisIt <https://wci.llnl.gov/codes/visit/>`_, and play around with visualizing
them.

Next, let us look at the anatomy of the simulation control file. It is broken
up into a few sections. We will look at each section in turn. Please open
:file:`simple.cpy` in an editor.

The first section of the control file is labeled "PIC setup", and contains just that::

    # -----------------------------------------------------------------------------
    # pic setup
    # -----------------------------------------------------------------------------
    pusher = PushMonomial()
    reconstructor = RecShape()

    dimensions_pos = 2
    dimensions_velocity = 2

    final_time = 10*units.M/units.VACUUM_LIGHT_SPEED

    vis_interval = 10

It sets up the solver's dimensionality, which deposition and push methods it
should use, how long it should run, and how often it should dump the
visualization output to disk.

Next up is the mesh setup::

    # -----------------------------------------------------------------------------
    # geometry and field discretization
    # -----------------------------------------------------------------------------
    element_order = 3

    _tube_width = 1*units.M
    _tube_length = 2*units.M

    import hedge.mesh as _mesh
    mesh = _mesh.make_rect_mesh(
            [0, -_tube_width/2],
            [_tube_length, _tube_width/2],
            periodicity=(True, False),
            subdivisions=(10,5),
            max_area=0.02)

Here we observe the first use of *user variables*. :mod:`pyrticle` wants to
interpret every setting in a control file as a simulation parameter. If it does
not something you are setting, it will stop with an error.  If a setting
starts with an underscore (like `_tube_length` above), however, this setting 
becomes a *user variable* and is ignored by :mod:`pyrticle`. This way, you can
perform your own symbolic calculations in your control files.

Apart from user variables, this section is reasonably straightforward--it sets
up an x-periodic rectangular domain.

The final section in the control file is concerned with setting up a particle
distribution::

    # -----------------------------------------------------------------------------
    # particle setup
    # -----------------------------------------------------------------------------
    _cloud_charge = 10e-9 * units.C
    nparticles = 2000

    _electrons_per_particle = abs(_cloud_charge/nparticles/units.EL_CHARGE)
    _pmass = _electrons_per_particle*units.EL_MASS

    _c0 = units.VACUUM_LIGHT_SPEED

    _mean_v = numpy.array([_c0*0.9,0])
    _sigma_v = numpy.array([_c0*0.9*1e-3, _c0*1e-5])

    _mean_beta = _mean_v/units.VACUUM_LIGHT_SPEED
    _gamma = units.gamma_from_v(_mean_v)
    _mean_p = _gamma*_pmass*_mean_v

    distribution = pyrticle.distribution.JointParticleDistribution([
        pyrticle.distribution.GaussianPos(
            [_tube_length*0.25,0.0], 
            [0.1, 0.1]),
        pyrticle.distribution.GaussianMomentum(
            _mean_p, _sigma_v*_gamma*_pmass, 
            units,
            pyrticle.distribution.DeltaChargeMass(
                _cloud_charge/nparticles,
                _pmass))
        ])

As you can see, it makes extensive use of user variables and finally sets
up a Gaussian distribution in position and momentum, coupled with delta 
distributions in particle charge and mass.

More documentation on the details of how control files are written will
become available as time permits, for now we invite you to tinker and 
ask questions. If something is unclear, please do not hesitate to get in
touch.

Thanks for trying Pyrticle!
