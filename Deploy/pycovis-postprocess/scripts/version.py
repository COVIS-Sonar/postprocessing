from pycovis.postprocess import runtime

with runtime.Runtime() as covis:

    print(" * Matlab version:               %s" % covis.matlab_version() )
    print(" * COVIS postprocessing version: %s" % covis.covis_version() )
