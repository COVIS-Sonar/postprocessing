from . import runtime

from pathlib import Path
import argparse
import re
import sys

import os

from urllib.parse import urlparse
import tempfile

from minio import Minio
from minio.error import ResponseError

from decouple import config
from pyunpack import Archive


def postprocessing_metadata():
    with runtime.Runtime() as covis:
        return covis.postproc_metadata()

def process( inputFile, outputDir ):

    ## Parse the input file,
    inputPath = Path( inputFile ).resolve()
    outputPath = Path( outputDir ).resolve()

    tempdirObj = tempfile.TemporaryDirectory(prefix = 'covis')
    tempdir = Path(tempdirObj.name)

    if not inputPath.exists():
        sys.exit( "Input file %s doesn't exist" % inputPath )

    basename = inputPath.stem
    basename = re.sub(r'\.tar','',basename)   ## Pathlib::stem only removes the .gz from .tar.gz

    matOutput = basename + ".mat"
    matOutputPath = outputPath / matOutput

    print("Processing input file: %s" % inputPath )
    print("            to output: %s" % matOutputPath )

    if inputPath.is_dir():
        # Pass through if it is already unpacked
        unpackedInput = inputPath
    else:
        ## \TODO.  More robust checking if it's an archive
        Archive( inputPath ).extractall( tempdir )

        for p in Path(tempdir).rglob("*/"):
            print(p)
            if p.is_dir():
                unpackedInput = p
                break

    if not unpackedInput:
        logging.error("Could not find unpacked input")
        return

    print("Unpacked input file at: %s" % unpackedInput )

    ## Tough to diagnose what's happening with the unpacked archive
    # for x in os.listdir( tempdir.name ):
    #     print(x)
    #
    # for x in os.listdir( unpackedInput ):
    #     print(x)


    with runtime.Runtime() as covis:
        result = covis.process( str(unpackedInput), str(matOutputPath) )
