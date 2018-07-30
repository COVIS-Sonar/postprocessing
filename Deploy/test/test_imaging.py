import os.path
import pytest
from pycovis.postprocess import runtime

## Process an imaging file through the specific function imaging_sweep
def test_imaging_sweep():
    testFile = "/input/APLUWCOVISMBSONAR001_20111001T210757.973Z-IMAGING.tar.gz"
    outputDir = "/output/test_imaging_sweep/"

    with runtime.Runtime() as pp:
        matfile = pp.covis_imaging_sweep(testFile, outputDir)
        assert( os.path.isfile(matfile) )

        plotfile = pp.covis_imaging_plot(matfile, outputDir)


## Process an imaging file through the generic function process_sweep
def test_imaging_process_sweep():
    testFile = "/input/APLUWCOVISMBSONAR001_20111001T210757.973Z-IMAGING.tar.gz"
    outputDir = "/output/test_imaging_process_sweep/"

    with runtime.Runtime() as pp:
        matfile = pp.covis_process_sweep(testFile, outputDir)
        assert( os.path.isfile(matfile) )

        plotfile = pp.covis_imaging_plot(matfile, outputDir)
