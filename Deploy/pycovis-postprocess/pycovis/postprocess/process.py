from . import runtime

from pathlib import Path
import argparse
import re
import sys

import io
import os
import logging

from urllib.parse import urlparse
import tempfile

from minio import Minio
from minio.error import ResponseError

from decouple import config
from pyunpack import Archive


class ProcessResult:

    def __init__(self, stdout):
        self.stdout = stdout


def postprocessing_metadata():
    with runtime.Runtime() as covis:
        return covis.postproc_metadata()

def process( inputFile, outputDir ):

    logger = logging.getLogger(__name__)

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

    logger.info("Processing input file: %s" % inputPath )
    logger.info("            to output: %s" % matOutputPath )

    if inputPath.is_dir():
        # Pass through if it is already unpacked
        unpackedInput = inputPath
    else:
        ## \TODO.  More robust checking if it's an archive
        Archive( inputPath ).extractall( tempdir )

        for p in Path(tempdir).rglob("*/"):
            if p.is_dir():
                unpackedInput = p
                break

    if not unpackedInput:
        logger.error("Could not find unpacked input")
        return

    logger.debug("Unpacked input file at: %s" % unpackedInput )

    out = io.StringIO()
    #err = io.StringIO()

    with runtime.Runtime() as covis:
        result = covis.process( str(unpackedInput), str(matOutputPath), stdout=out, stderr=out )

    out.seek(0)
    return ProcessResult( out )
