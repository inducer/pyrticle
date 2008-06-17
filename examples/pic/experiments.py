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

def compare_methods():
    """Submit jobs to compare reconstruction/pushing methods."""

    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    for rec in [
        O("RecGrid", el_tolerance=0.2),
        O("RecAdv"), 
        O("RecNormShape"), 
        O("RecShape"), 
        ]:
        for eorder in [3,4,5,7]:
            if rec.classname == "RecGrid":
                pushers = [O("PushMonomial")]
            else:
                pushers = [O("PushMonomial"), O("PushAverage")]
            for pusher in pushers:
                for finder in [
                        O("FindFaceBased"),
                        #O("FindHeuristic"),
                        ]:
                    job = BatchJob(
                            "compmeth-$DATE/eo%d-%s-%s-%s" % (eorder, cn(rec), cn(pusher), cn(finder)),
                            "with-charge.py",
                            timestamp=timestamp,
                            )
                    job.write_setup([
                        "pusher = %s" % pusher,
                        "reconstructor = %s" % rec,
                        "finder = %s" % finder,
                        "element_order = %d" % eorder,
                        ])
                    job.submit()

def study_rec_grid(output_path=None):
    """Submit jobs to study the behavior of grid reconstruction."""

    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    method = "simplex_enlarge"
    pusher = O("PushMonomial")
    eorder = 3

    nparticles = 30000

    for method in ["simplex_enlarge", "simplex_extra", "brick"]:
        for el_tolerance in [0, 0.1, 0.15, 0.2, 0.3]:
            for enf_cont in [True, False]:
                for overres in [0.8, 1.0, 1.2, 1.3, 1.4, 1.6, 2]:
                    for mesh_margin in [0, 0.05, 0.1, 0.2]:
                        job = BatchJob(
                                "recgrid-$DATE/%s-tol%g-cont%s-or%g-mm%g" % (
                                    method, el_tolerance, enf_cont, overres, mesh_margin),
                                "with-charge.py",
                                timestamp=timestamp,
                                )
                        brick_gen = O("SingleBrickGenerator",
                                overresolve=overres,
                                mesh_margin=mesh_margin)
                        rec = O("RecGrid", brick_gen,
                                el_tolerance=el_tolerance,
                                enforce_continuity=enf_cont,
                                method=method)

                        if output_path is not None:
                            import os
                            job_out_path = os.path.join(output_path, job.subdir)
                            os.makedirs(job_out_path)

                        setup = [
                            "pusher = %s" % pusher,
                            "reconstructor = %s" % rec,
                            "element_order = %d" % eorder,
                            "nparticles = %d" % nparticles,
                            "vis_path = %s" % repr(job_out_path),
                            ]
                        job.write_setup(setup)
                        job.submit()




def compare_element_finders():
    """Submit jobs to compare element finders."""

    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    reconstructor = O("RecShape")
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
                "reconstructor = %s" % reconstructor,
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
      #O("RecAdv"), 
      #O("RecNormShape"), 
      O("RecShape")
      ]:
        for chi in [None, 1, 2, 5]:
            if chi is None:
                filters = [None]
            else:
                filters = [
                  None,
                  (0.6,6),
                  (0.8,6),
                  (0.95,6),
                  (0.6,3),
                  (0.8,3),
                  (0.95,3),
                  ]
            for filter in filters:
                for pusher in [O("PushAverage")]:
                        job = BatchJob(
                                "cleanstudy-$DATE/%s-chi%s-filt%s" % (cn(rec), chi, filt_desc(filter)),
                                "with-charge.py",
                                timestamp=timestamp,
                                )
                        job.write_setup([
                            "pusher = %s" % pusher,
                            "reconstructor = %s" % rec,
                            "chi = %s" % chi,
                            "phi_filter = %s" % str(filter),
                            "element_order = 4",
                            ])
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
                "reconstructor = %s" % O("RecAdv", 
                    filter_order=filter_order, 
                    filter_amp=filter_amp),
                ])
            job.submit()

def study_blob_exponent():
    """Submit jobs to study the effect of the exponent in the shape blob."""

    O = ConstructorPlaceholder

    timestamp = get_timestamp()

    for exponent in [1,2,3,4,5,6]:
        #for rec in [O("RecAdv"), O("RecNormShape"), O("RecNormShape")]:
        #for rec in [O("RecShape"), ]:
        for eorder in [2,3,4,5]:
            for rec in [O("RecShape"), ]:
                for push in [O("PushMonomial"), O("PushAverage")]:
                    job = BatchJob(
                            "expstudy-$DATE/exp%d-eo%d-%s-%s" % (exponent, eorder, cn(rec), cn(push)),
                            "with-charge.py",
                            timestamp=timestamp,
                            )
                    job.write_setup([
                        "pusher = %s" % push,
                        "reconstructor = %s" % rec,
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

import sys
exec sys.argv[1]
