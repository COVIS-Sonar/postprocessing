from pycovis.postprocess import runtime

import logging
from pathlib import Path
import argparse

parser = argparse.ArgumentParser(description='Process one COVIS file.')
parser.add_argument('inputFile', help='COVIS file to parse')
parser.add_argument("--output",  help="Filename for output",
                        dest="outputMat", default="")

args = parser.parse_args()


inputPath = Path(args.inputFile).resolve()
outputPath = Path(args.outputMat).resolve()

if not inputPath.exists():
    sys.exit( "Input file %s doesn't exist" % inputPath )

print("Processing input file: %s" % inputPath )

with runtime.Runtime() as covis:

    result = covis.process( str(inputPath), str(outputPath) )

    # if result:
    #     print("Saving results to %s" % result)
    #     covis.save( args.output, result )

    # matfile = pp.covis_imaging_sweep(input, outdir, '')
    # imgfile = pp.covis_imaging_plot(matfile, outdir, '')
