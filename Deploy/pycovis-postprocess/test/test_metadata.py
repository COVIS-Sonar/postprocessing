
import pycovis.matlab as pymatlab
from pycovis.postprocess import runtime

def test_imaging_sweep():
    with runtime.Runtime() as pp:
        metadata = pp.postproc_metadata()

        for key in ["matlab_version", "verstr", "postprocessing_gitrev", "postprocessing_gittags"]:
            assert key in metadata
