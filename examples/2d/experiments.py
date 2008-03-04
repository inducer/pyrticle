from pytools.batchjob import GridEngineJob, ConstructorPlaceholder

BatchJob = GridEngineJob

def compare_methods():
    """Submit jobs to compare reconstruction/pushing methods."""

    O = ConstructorPlaceholder

    for rec in [O("RecAdv"), O("RecNormShape"), O("RecShape")]:
        for pusher in [O("PushMonomial"), O("PushAverage")]:
            job = BatchJob(
                    "compmeth-%s-%s" % (rec.classname, pusher.classname),
                    "with-charge.py",
                    )
            job.write_setup([
                "pusher = %s" % pusher,
                "reconstructor = %s" % rec,
                ])

def study_cleaning():
    """Submit jobs to see the effect of hyperbolic cleaning."""

    O = ConstructorPlaceholder

    for rec in [O("RecAdv"), O("RecNormShape")]:
        for pusher in [O("PushAverage")]:
            for chi in [None, 1, 2, 5]:
                job = BatchJob(
                        "cleanstudy-%s-chi%s" % (rec.classname, chi),
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
            return [3,6]
    for filter_amp in [None, 0.97, 0.93]:
        for filter_order in get_filter_orders(filter_amp):
            job = BatchJob(
                    "filtstudy-amp%s-ord%s" % (filter_amp, filter_order),
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
                    "expstudy-exp%d-%s" % (exponent, rec.classname),
                    "with-charge.py",
                    )
            job.write_setup([
                "pusher = PushAverage()",
                "reconstructor = %s" % rec,
                "shape_exponent = %s" % exponent,
                ])
            job.submit()

study_blob_exponent()
study_advec_filter()
study_cleaning()
