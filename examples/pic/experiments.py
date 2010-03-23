from pytools.batchjob import \
        guess_job_class, \
        ConstructorPlaceholder, \
        get_timestamp

BatchJob = guess_job_class()

def cn(placeholder):
    """Chop the first bit off a CamelCasedClassName, and convert to lower case."""

    result = placeholder.classname
    assert result[0].isupper()

    for i in range(1, len(result)):
        if result[i].isupper():
            return result[i:].lower()
    else:
        return result.lower()

def cn_with_args(placeholder):
    result = cn(placeholder)
    if placeholder.args:
        result += "-" + ("-".join(str(a) for a in args))

    def mogrify_kwarg(k, v):
        result = k[0]
        
        i = 1
        while i < len(k):
            if k[i] == "_":
                result += k[i+1]
                i += 2
            else:
                i += 1

        return result+str(v)

    if placeholder.kwargs:
        result += "-" + ("-".join(
            mogrify_kwarg(k, v) 
            for k, v in placeholder.kwargs.iteritems()))

    return result

def multiline_to_setup(mls):
    lines = mls.split("\n")

    assert lines[0] == ""
    lines.pop(0)

    first_line = lines[0]
    indent = 0
    while first_line[indent] == " " and indent < len(first_line):
        indent += 1

    indent_str = indent*" "

    def remove_indent(s):
        assert s.startswith(indent_str) or len(s) == 0
        return s[indent:]

    return [remove_indent(line) for line in lines]




# basic setups ----------------------------------------------------------------
def basic_2d_gauss_setup():
    return multiline_to_setup("""
    import random as _random
    _random.seed(0)

    dimensions_pos = 2
    dimensions_velocity = 2

    beam_axis = 0
    beam_diag_axis = 1
    tube_length = 2

    _cloud_charge = -10e-9 * units.C
    final_time = 10*units.M/units.VACUUM_LIGHT_SPEED
    _electrons_per_particle = abs(_cloud_charge/nparticles/units.EL_CHARGE)

    _tube_width = 1
    import hedge.mesh as _mesh
    mesh = _mesh.make_rect_mesh(
            a=(-0.5, -_tube_width/2),
            b=(-0.5+tube_length, _tube_width/2),
            periodicity=(True, False),
            subdivisions=(10,5),
            max_area=0.02)

    _c0 = units.VACUUM_LIGHT_SPEED

    _mean_v = numpy.array([_c0*0.9,0])
    _sigma_v = numpy.array([_c0*0.9*1e-3, _c0*1e-5])

    _mean_beta = _mean_v/units.VACUUM_LIGHT_SPEED
    _gamma = units.gamma_from_v(_mean_v)
    _pmass = _electrons_per_particle*units.EL_MASS
    _mean_p = _gamma*_pmass*_mean_v

    distribution = pyrticle.distribution.JointParticleDistribution([
        pyrticle.distribution.GaussianPos([0,0], [0.1, 0.1]),
        pyrticle.distribution.GaussianMomentum(
            _mean_p, _sigma_v*_gamma*_pmass, units,
            pyrticle.distribution.DeltaChargeMass(
                _cloud_charge/nparticles,
                _pmass))
        ])
    """)




# experiments -----------------------------------------------------------------
def compare_methods():
    """Submit jobs to compare deposition/pushing methods."""

    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    for rec in [
        #O("DepGrid", jiggle_radius=0),
        #O("DepGrid"),
        #O("DepGridFind"),
        O("DepAdv"), 
        #O("DepNormShape"), 
        #O("DepShape"), 
        ]:
        for eorder in [3]:
            for sexp in [3]:
                if "Grid" in rec.classname:
                    pushers = [O("PushMonomial")]
                else:
                    pushers = [
                            O("PushMonomial"), 
                            #O("PushAverage")
                            ]
                for pusher in pushers:
                    job = BatchJob(
                            "compmeth-$DATE/eo%d-se%d-%s-%s" % (
                                eorder, sexp, cn_with_args(rec), cn(pusher)),
                            "driver.py",
                            timestamp=timestamp,
                            )
                    job.write_setup([
                        "pusher = %s" % pusher,
                        "depositor = %s" % rec,
                        "element_order = %d" % eorder,
                        "shape_exponent = %d" % sexp,
                        ]
                        +basic_2d_gauss_setup()
                        )
                    job.submit()

def study_rec_grid(output_path=None):
    """Submit jobs to study the behavior of grid deposition."""

    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    pusher = O("PushMonomial")
    eorder = 3

    nparticles = 30000


    for method in ["simplex_enlarge", "simplex_extra", "simplex_reduce"]:
        #for el_tolerance in [0, 0.1, 0.15, 0.2]:
        for el_tolerance in [0.1, 0.15]:
            #for enf_cont in [True, False]:
            for enf_cont in [False]:
                #for overres in [0.8, 1.0, 1.2, 1.3, 1.4, 1.6, 2]:
                for overres in [1.2, 1.4, 1.6, 2]:
                    for mesh_margin in [0]:
                    #for mesh_margin in [0, 0.05, 0.1, 0.2]:
                        job = BatchJob(
                                "recgrid-$DATE/%s-tol%g-cont%s-or%g-mm%g" % (
                                    method, el_tolerance, enf_cont, overres, mesh_margin),
                                "driver.py",
                                timestamp=timestamp,
                                )
                        brick_gen = O("SingleBrickGenerator",
                                overresolve=overres,
                                mesh_margin=mesh_margin)
                        rec = O("DepGrid", brick_gen,
                                el_tolerance=el_tolerance,
                                enforce_continuity=enf_cont,
                                method=method)

                        setup = [
                            "pusher = %s" % pusher,
                            "depositor = %s" % rec,
                            "element_order = %d" % eorder,
                            "nparticles = %d" % nparticles,
                            ]+basic_2d_gauss_setup()

                        if output_path is not None:
                            import os
                            job_out_path = os.path.join(output_path, job.subdir)
                            os.makedirs(job_out_path)

                            setup.append(
                                    "vis_path = %s" % repr(job_out_path),
                                    )
                        job.write_setup(setup)
                        job.submit()




def compare_element_finders():
    """Submit jobs to compare element finders."""

    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    depositor = O("DepShape")
    pusher = O("PushMonomial")
    for xpos in [0, 0.5, 1]:
        for finder in [O("FindFaceBased"), O("FindHeuristic")]:
            job = BatchJob(
                    "compelfind-$DATE/%s-%s" % (xpos, cn(classname)),
                    "rec-by-area.py",
                    ["special_meshes.py"],
                    timestamp=timestamp,
                    )
            job.write_setup([
                "pusher = %s" % pusher,
                "depositor = %s" % depositor,
                "finder = %s" % finder,
                "mean_x = num.array([%g*tube_length,0])" % xpos,
                "nparticles = 1000",
                "mesh = 'glued'",
                ])
            job.submit()

def study_cleaning():
    """Submit jobs to see the effect of hyperbolic cleaning."""

    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    def filt_desc(f):
        if isinstance(f, tuple):
            return "-".join(str(i) for i in f)
        else:
            return str(f)

    for rec in [
      #O("DepAdv"), 
      #O("DepNormShape"), 
      O("DepShape"),
      O("DepGrid", jiggle_radius=0),
      ]:
        for chi in [None, 5]:
            if chi is None:
                filters = [None]
            else:
                filters = [
                  None,
                  #(0.6,6),
                  #(0.8,6),
                  #(0.95,6),
                  #(0.6,3),
                  #(0.8,3),
                  #(0.95,3),
                  ]
            for filter in filters:
                for pusher in [O("PushMonomial")]:
                        job = BatchJob(
                                "cleanstudy-$DATE/%s-chi%s-filt%s" % (cn(rec), chi, filt_desc(filter)),
                                "driver.py",
                                timestamp=timestamp,
                                )
                        job.write_setup([
                            "pusher = %s" % pusher,
                            "depositor = %s" % rec,
                            "chi = %s" % chi,
                            "phi_filter = %s" % str(filter),
                            "element_order = 4",
                            ]+basic_2d_gauss_setup())
                        job.submit()

def study_advec_filter():
    """Submit jobs to see whether filtering really improves advection."""

    O = ConstructorPlaceholder

    def get_filter_orders(amp):
        if filter_amp is None: 
            return [None]
        else:
            return [3, 6, 7, 9]

    timestamp = get_timestamp()

    for filter_amp in [None, 0.97, 0.93, 0.88, 0.8, 0.7, 0.5]:
        for filter_order in get_filter_orders(filter_amp):
            job = BatchJob(
                    "filtstudy-$DATE/amp%s-ord%s" % (filter_amp, filter_order),
                    "with-charge.py",
                    timestamp=timestamp,
                    )
            job.write_setup([
                "pusher = %s" % O("PushAverage"),
                "depositor = %s" % O("DepAdv", 
                    filter_order=filter_order, 
                    filter_amp=filter_amp),
                ])
            job.submit()

def study_blob_exponent():
    """Submit jobs to study the effect of the exponent in the shape blob."""

    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    for exponent in [1,2,3,4,5,6]:
        #for rec in [O("DepAdv"), O("DepNormShape"), O("DepNormShape")]:
        #for rec in [O("DepShape"), ]:
        for eorder in [2,3,4,5]:
            for rec in [O("DepShape"), ]:
                for push in [O("PushMonomial"), O("PushAverage")]:
                    job = BatchJob(
                            "expstudy-$DATE/exp%d-eo%d-%s-%s" % (exponent, eorder, cn(rec), cn(push)),
                            "with-charge.py",
                            timestamp=timestamp,
                            )
                    job.write_setup([
                        "pusher = %s" % push,
                        "depositor = %s" % rec,
                        "shape_exponent = %s" % exponent,
                        "element_order = %s" % eorder,
                        ])
                    job.submit()

def run_apsgun():
    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    job = BatchJob(
            "apsgun/$DATE",
            "driver.py",
            aux_files=["apsgun.cpy"],
            timestamp=timestamp,
            )
    job.write_setup([ "execfile('apsgun.cpy')" ])
    job.submit()

def run_kv3d():
    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    for rec in [
            O("DepShape"), 
            O("DepGrid", O("FineCoreBrickGenerator", core_axis=2),
                el_tolerance=0.1,
                method="simplex_reduce")
            ]:
        job = BatchJob(
                "kv3d-$DATE/%s" % cn(rec),
                "driver.py",
                aux_files=["kv3d.cpy"],
                timestamp=timestamp,
                )
        job.write_setup([ 
            "execfile('kv3d.cpy')",
            "depositor = %s" % rec,
            ])
        job.submit()

import sys
exec sys.argv[1]
