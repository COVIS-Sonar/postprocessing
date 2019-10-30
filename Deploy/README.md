This directory contains scripts and tools related to packaging the `master_program` Matlab code into a Python package using the Matlab compiler, then building this package into a Docker image `amarburg/covis-postprocess` for distributed execution.

The general steps are:

1. `make matlab` calls the Matlab compiler (which must be installed locally) to convert `master_program` into the Python library `pycovis-matlab`.   The (normal, not auto-generated) `pycovis-postprocess` provides more Pythonic wrappers around the Matlab-generated code.

1. `make docker` then installs `pycovis-matlab` and  `pycovis-postprocess` into a Docker image built on the base image `amarburg/matlab_runtime`, a Debian-based image which includes the Matlab Runtime.

1.  The resulting `amarburg/covis-postprocess` Docker image can then be run "anywhere."    

    * Running `make pytest` will execute the Python test suite found in `test/` _in the Docker image._   
    * Similarly, `make version` will run the simple Python+Matlab script `scripts/version.py` _in the Docker image._

Alternatively, once the `pycovis-matlab` library has been built, it can be executed on the local machine _if Matlab or the Matlab runtime are installed._   The LOCAL_MATLAB_INSTALL env variable my need to be set (it defaults to `/usr/local/MATLAB/R{current Matlab version}`)    The Makefile jobs `make local_pytest` and `make local_version` are the local equivalents to `mae pytest / version` above.


All tasks are scripted in the `Makefile`.   Run `make help` to get a list of all commands.

## Python modules and directories

`pycovis-matlab` is the Matlab compiler-generated helper code.

`pycovis-postprocess` is "normal" (non-Matlab-generated) Python wrapper around `pycovis-matlab`.  It should be lightweight.  It contains:

   * `scripts/` which are Python script entrypoints for common tasks.
   * `test/` is a unit test suite.
