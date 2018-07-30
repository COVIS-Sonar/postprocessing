This directory contains scripts and tools related to packaging this Matlab code
into a Python package, then building this package into a Docker
image `amarburg/covis-postprocess:latest` for distributed execution.

Tasks are documented in the `Makefile`.   Run `make help` to get a
list of commands.

`pycovis-matlab` is Matlab-geneated helper code.

`pycovis-postprocess` is "normal" (non-Matlab-generated) Python helper code.
It should be lightweight.
