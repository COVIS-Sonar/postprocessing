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


def process( inputFile, outputDir ):

    ## Parse the input file,
    inputPath = Path( inputFile ).resolve()
    outputPath = Path( outputDir ).resolve()

    tempdir = tempfile.TemporaryDirectory(prefix = 'covis')

    if not inputPath.exists():
        sys.exit( "Input file %s doesn't exist" % inputPath )

    basename = inputPath.stem
    basename = re.sub(r'\.tar','',basename)   ## Pathlib::stem only removes the .gz from .tar.gz

    matOutput = basename + ".mat"
    matOutputPath = outputPath / matOutput

    print("Processing input file: %s" % inputPath )
    print("            to output: %s" % matOutputPath )

    ## Unpack archive in Python
    if inputPath.is_dir():
        unpackedInput = inputPath
    else:
        ## \TODO.  More robust checking if it's an archive
        Archive(inputPath).extractall( tempdir.name )

        ## \TODO  Need a way to find the directory name found in the archive rather
        ## than this open-loop version
        unpackedInput = Dir.glob( tempdir.name + "/*" )[0]  / basename

    ## Tough to diagnose what's happening with the unpacked archive
    # for x in os.listdir( tempdir.name ):
    #     print(x)
    #
    # for x in os.listdir( unpackedInput ):
    #     print(x)


    with runtime.Runtime() as covis:
        result = covis.process( str(unpackedInput), str(matOutputPath) )
