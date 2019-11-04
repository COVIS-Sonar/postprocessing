from pycovis.postprocess import runtime

import logging
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description='Process one COVIS file.')
parser.add_argument('inputFile', help='COVIS file to parse')
parser.add_argument("--output",  help="Filename for output", default="")

args = parser.parse_args()


inputPath = Path(args.inputFile).resolve()

if not inputPath.exists():
    sys.exit( "Input file %s doesn't exist" % inputPath )

print("Processing input file: %s" % inputPath )

with runtime.Runtime() as covis:

    result = covis.covis_raw_sweep( str(inputPath.parent), str(inputPath.name), 0 )

    if result:
        print("Saving results to %s" % result)
        covis.save( args.output, result )

    # matfile = pp.covis_imaging_sweep(input, outdir, '')
    # imgfile = pp.covis_imaging_plot(matfile, outdir, '')
