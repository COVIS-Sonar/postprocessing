# COVIS Post-Processing

Matlab code for parsing and post-processing data from the COVIS instrument as
installed on the ONC Neptune Cabled Array at the Endeavour vent field from
2010-2015 and to be installed on the OOI Cabled Array at the Ashes vent field 2018-2022.

The `Diffuse/` `Doppler/` and `Imaging/` directories contain Matlab functions
related to those three modes of COVIS.  `Common/` contains Matlab code shared
by all the three modes.

`Test/` contains a Matlab-language test suite.    On a (Linux or Mac?) machine
with Matlab and Gnu `make` installed:

  * `make test` or `make unittest` will run a (relatively fast) unit test suite.
  * `make integrationtest` will run a more thorough integration test suite.  This will
        fully process sample COVIS files from each mode and
        may take some time to run.


The `Deploy/` directory contains code related to packaging this code as a
Python library for distributed processing.
