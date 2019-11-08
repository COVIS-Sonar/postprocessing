from pycovis.postprocess import runtime

import logging
from pathlib import Path
import argparse
import re
import sys

from urllib.parse import urlparse
import tempfile

from minio import Minio
from minio.error import ResponseError

from decouple import config


parser = argparse.ArgumentParser(description='Process one COVIS file.')
parser.add_argument('inputFile', help='COVIS file to parse')
parser.add_argument("--output",  help="Directory for output",
                        dest="outputDir", default="")

parser.add_argument('--input-s3-host', default=config('INPUT_S3_HOST', default=""),
                    dest="inputS3Host",
                    help='Hostname for input S3 server (if s3:// URL is provided for inputFile)')
parser.add_argument('--input-s3-access-key', default=config('INPUT_S3_ACCESS_KEY', default=""),
                    dest="inputS3AccessKey",
                    help='Access key for input S3 server (if s3:// URL is provided for inputFile)')
parser.add_argument('--input-s3-secret-key', default=config('INPUT_S3_SECRET_KEY', default=""),
                    dest="inputS3SecretKey",
                    help='Secret key for input S3 server (if s3:// URL is provided for inputFile)')

parser.add_argument('--output-s3-host', default=config('OUTPUT_S3_HOST', default=""),
                    dest="outputS3Host",
                    help='Hostname for output S3 server (if s3:// URL is provided for outputFile)')
parser.add_argument('--output-s3-access-key', default=config('OUTPUT_S3_ACCESS_KEY', default=""),
                    dest="outputS3AccessKey",
                    help='Access key for output S3 server (if s3:// URL is provided for outputFile)')
parser.add_argument('--output-s3-secret-key', default=config('OUTPUT_S3_SECRET_KEY', default=""),
                    dest="outputS3SecretKey",
                    help='Secret key for output S3 server (if s3:// URL is provided for outputFile)')


args = parser.parse_args()


## Parse the input file,
input = urlparse( args.inputFile )
output = urlparse( args.outputDir )

## Validate the output first, rather than catch it at the end
if output.scheme == '':
    outputPath = args.outputDir
    if not outputPath.exists():
        sys.exit("Specified output path %s does not exist" % outputPath)
elif output.scheme == 's3':
    s3host = args.outputS3Host
    if not s3host:
        sys.exit("s3:// output URL provided but --output-s3-host not provided")

    outputMinioClient = Minio(s3host,
                  access_key=args.outputS3AccessKey,
                  secret_key=args.outputS3SecretKey,
                  secure=False)

    ## ...


## Now deal with the input
## If a local file
if input.scheme == '':
    inputPath = Path(args.inputFile).resolve()
elif input.scheme == 's3':

    s3host = args.inputS3Host
    if not s3host:
        sys.exit("s3:// input URL provided but --input-s3-host not provided")

    bucket = input.netloc
    path   = input.path

    print("Retrieving bucket: %s, path %s from host %s" % (bucket, path, s3host))

    tempdir = tempfile.TemporaryDirectory(prefix = 'covis')

    minioClient = Minio(s3host,
                  access_key=args.inputS3AccessKey,
                  secret_key=args.inputS3SecretKey,
                  secure=False)


    basename = Path(path).name
    inputPath = Path( tempdir.name ) / basename

    print("Retrieving path %s from bucket %s to %s" % (path, bucket, inputPath))

    minioClient.fget_object(bucket_name=bucket, object_name=path, file_path=str(inputPath) )

# if args.outputMat:
#     outputPath = Path(args.outputMat).resolve()
# else:
if True:

    inputName = inputPath.name

    # Ugly for now
    inputName = re.sub(r'\.tar\.gz','.mat',inputName)
    outputName = re.sub(r'\.7z','.mat',inputName)

    outputPath = Path("/output") / outputName

if not inputPath.exists():
    sys.exit( "Input file %s doesn't exist" % inputPath )

print("Processing input file: %s" % inputPath )
print("            to output: %s" % outputPath)

with runtime.Runtime() as covis:

    result = covis.process( str(inputPath), str(outputPath) )

## Reuse from above
if outputMinioClient:

    outbucket = output.netloc

    if not outputMinioClient.bucket_exists( outbucket ):
        outputMinioClient.make_bucket(outbucket)

    ## \TODO  generate correct outputName
    outpath = Path(output.path) / outputName
    outpath = str(outpath).strip("/")  ## Minio client doesn't like leading slash
    outputMinioClient.fput_object(file_path=outputPath,bucket_name=outbucket,object_name=outpath)

    print("Uploaded to bucket %s, path %s" % (outbucket,outpath) )

    # if result:
    #     print("Saving results to %s" % result)
    #     covis.save( args.output, result )

    # matfile = pp.covis_imaging_sweep(input, outdir, '')
    # imgfile = pp.covis_imaging_plot(matfile, outdir, '')
