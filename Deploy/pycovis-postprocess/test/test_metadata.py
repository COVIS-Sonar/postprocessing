
import pycovis.matlab as pymatlab
from pycovis.postprocess import runtime

def test_imaging_sweep():
    with runtime.Runtime() as pp:
        metadata = pp.postproc_metadata()

        for key in ["gitrev", "gittags"]:
            assert key in metadata
