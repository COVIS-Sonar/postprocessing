from pycovis.postprocess import runtime

from pathlib import Path
import argparse


parser = argparse.ArgumentParser(description='Process one COVIS file.')
parser.add_argument('inputFile', help='an integer for the accumulator')
# parser.add_argument('--sum', dest='accumulate', action='store_const',
#                     const=sum, default=max,
#                     help='sum the integers (default: find the max)')

args = parser.parse_args()


inputPath = "../covis-test-data/"
#inputFile = "COVIS-20191024T003002-diffuse1.tar.gz"

with runtime.Runtime() as covis:

    result = covis.covis_raw_sweep( inputPath, args.inputFile, 0 )
    # matfile = pp.covis_imaging_sweep(input, outdir, '')
    # imgfile = pp.covis_imaging_plot(matfile, outdir, '')
