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

#parser = argparse.ArgumentParser(description='Process one COVIS file.')

# parser.add_argument('inputFile', help='COVIS file to parse')
# parser.add_argument("--output",  help="Directory for output",
#                         dest="outputDir", default="/output")
#
# parser.add_argument("--auto-output-path", action='store_true',
#                     help="Automatically add the YYYY/MM/DD/ path to the output")
#
# parser.add_argument('--input-s3-host', default=config('INPUT_S3_HOST', default=""),
#                     dest="inputS3Host",
#                     help='Hostname for input S3 server (if s3:// URL is provided for inputFile)')
# parser.add_argument('--input-s3-access-key', default=config('INPUT_S3_ACCESS_KEY', default=""),
#                     dest="inputS3AccessKey",
#                     help='Access key for input S3 server (if s3:// URL is provided for inputFile)')
# parser.add_argument('--input-s3-secret-key', default=config('INPUT_S3_SECRET_KEY', default=""),
#                     dest="inputS3SecretKey",
#                     help='Secret key for input S3 server (if s3:// URL is provided for inputFile)')
#
# parser.add_argument('--output-s3-host', default=config('OUTPUT_S3_HOST', default=""),
#                     dest="outputS3Host",
#                     help='Hostname for output S3 server (if s3:// URL is provided for outputFile)')
# parser.add_argument('--output-s3-access-key', default=config('OUTPUT_S3_ACCESS_KEY', default=""),
#                     dest="outputS3AccessKey",
#                     help='Access key for output S3 server (if s3:// URL is provided for outputFile)')
# parser.add_argument('--output-s3-secret-key', default=config('OUTPUT_S3_SECRET_KEY', default=""),
#                     dest="outputS3SecretKey",
#                     help='Secret key for output S3 server (if s3:// URL is provided for outputFile)')
#
# args = parser.parse_args()


    ## Parse the input file,
    input = urlparse( inputFile )
    output = urlparse( outputDir )

    outputMinioClient = None

    ## Validate the output first, rather than catch it at the end
    if output.scheme == '':
        ## If writing to file, work directly in the final destination
        workingOutputPath = Path(output.path)

    elif output.scheme == 's3':
        s3host = config('OUTPUT_S3_HOST', default="")
        if not s3host:
            sys.exit("s3:// output URL provided but --output-s3-host not provided")

        outputMinioClient = Minio(s3host,
                      access_key=config('OUTPUT_S3_ACCESS_KEY', default=""),
                      secret_key=config('OUTPUT_S3_SECRET_KEY', default=""),
                      secure=False)

        ## Otherwise, set a default working directory
        workingOutputPath = Path("/output")

        ## ...


    tempdir = tempfile.TemporaryDirectory(prefix = 'covis')

    ## Now deal with the input
    ## If a local file
    if input.scheme == '':
        inputPath = Path(input.path).resolve()
    elif input.scheme == 's3':

        s3host = config('INPUT_S3_HOST', default="")
        if not s3host:
            sys.exit("s3:// input URL provided but --input-s3-host not provided")

        bucket = input.netloc
        path   = input.path

        print("Retrieving bucket: %s, path %s from host %s" % (bucket, path, s3host))

        minioClient = Minio(s3host,
                      access_key=config('INPUT_S3_ACCESS_KEY', default=""),
                      secret_key=config('INPUT_S3_SECRET_KEY', default=""),
                      secure=False)

        basename = Path(path).name
        inputPath = Path( tempdir.name ) / basename

        print("Retrieving path %s from bucket %s to %s" % (path, bucket, inputPath))

        minioClient.fget_object(bucket_name=bucket, object_name=path, file_path=str(inputPath) )

    if not inputPath.exists():
        sys.exit( "Input file %s doesn't exist" % inputPath )

    basename = inputPath.stem
    basename = re.sub(r'\.tar','',basename)   ## Pathlib::stem only removes the .gz from .tar.gz

    matOutput = basename + ".mat"
    matOutputPath = workingOutputPath / matOutput

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
        unpackedInput = Path(tempdir.name) / basename

    ## Tough to diagnose what's happening with the unpacked archive
    # for x in os.listdir( tempdir.name ):
    #     print(x)
    #
    # for x in os.listdir( unpackedInput ):
    #     print(x)


    with runtime.Runtime() as covis:
        result = covis.process( str(unpackedInput), str(matOutputPath) )

    ## Reuse from above
    if outputMinioClient:

        outbucket = output.netloc

        if not outputMinioClient.bucket_exists( outbucket ):
            outputMinioClient.make_bucket(outbucket)

        ## \TODO  generate correct outputName
        outpath = Path(output.path) / matOutput
        outpath = str(outpath).strip("/")  ## Minio client doesn't like leading slash

        outputMinioClient.fput_object(file_path=matOutputPath, bucket_name=outbucket, object_name=outpath)

        print("Uploaded to bucket %s, path %s" % (outbucket,outpath) )

        # if result:
        #     print("Saving results to %s" % result)
        #     covis.save( args.output, result )

        # matfile = pp.covis_imaging_sweep(input, outdir, '')
        # imgfile = pp.covis_imaging_plot(matfile, outdir, '')
