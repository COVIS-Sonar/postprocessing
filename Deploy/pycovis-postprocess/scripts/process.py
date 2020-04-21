from pycovis.postprocess import process
import logging
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


parser = argparse.ArgumentParser(description='Process one COVIS file.')

parser.add_argument('inputFile', help='COVIS file to parse')
parser.add_argument("--output",  help="Directory for output",
                        dest="outputDir", default="/output")

logging.basicConfig(level=logging.DEBUG)

args = parser.parse_args()

result = process.process( args.inputFile, args.outputDir )

for line in result.stdout:
    logging.info(line.rstrip())
