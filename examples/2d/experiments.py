from pytools.batchjob import guess_job_class, ConstructorPlaceholder

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

    for rec in [O("RecAdv"), O("RecNormShape"), O("RecShape")]:
        for pusher in [O("PushMonomial"), O("PushAverage")]:
            for finder in [
                    O("FindFaceBased"),
                    #O("FindHeuristic"),
                    ]:
                job = BatchJob(
                        "compmeth-$DATE/%s-%s-%s" % (cn(rec), cn(pusher), cn(finder)),
                        "with-charge.py",
                        )
                job.write_setup([
                    "pusher = %s" % pusher,
                    "reconstructor = %s" % rec,
                    "finder = %s" % finder,
                    ])
                job.submit()

def compare_element_finders():
    """Submit jobs to compare element finders."""

    O = ConstructorPlaceholder

    reconstructor = O("RecShape")
    pusher = O("PushMonomial")
    for xpos in [0, 0.5, 1]:
        for finder in [O("FindFaceBased"), O("FindHeuristic")]:
            job = BatchJob(
                    "compelfind-$DATE/%s-%s" % (xpos, cn(classname)),
                    "rec-by-area.py",
                    ["special_meshes.py"]
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

    for rec in [O("RecAdv"), O("RecNormShape"), O("RecShape")]:
        for pusher in [O("PushAverage")]:
            for chi in [None, 1, 2, 5]:
                job = BatchJob(
                        "cleanstudy-$DATE/%s-chi%s" % (cn(rec), chi),
                        "with-charge.py",
                        )
                job.write_setup([
                    "pusher = %s" % pusher,
                    "reconstructor = %s" % rec,
                    "chi = %s" % chi,
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

    for filter_amp in [None, 0.97, 0.93, 0.88, 0.8, 0.7, 0.5]:
        for filter_order in get_filter_orders(filter_amp):
            job = BatchJob(
                    "filtstudy-$DATE/amp%s-ord%s" % (filter_amp, filter_order),
                    "with-charge.py",
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

    for exponent in [1,2,3,4]:
        for rec in [O("RecAdv"), O("RecNormShape")]:
            job = BatchJob(
                    "expstudy-$DATE/exp%d-%s" % (exponent, cn(rec)),
                    "with-charge.py",
                    )
            job.write_setup([
                "pusher = PushAverage()",
                "reconstructor = %s" % rec,
                "shape_exponent = %s" % exponent,
                ])
            job.submit()

import sys
exec sys.argv[1]
