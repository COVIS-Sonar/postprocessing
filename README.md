# COVIS Post-Processing

Matlab code for parsing and post-processing data from the COVIS instrument as installed on the ONC Neptune Cabled Array at the Endeavour vent field from 2010-2015 and on the OOI Cabled Array at the Ashes vent field 2018-2022.

`master_program/` ... contains Matlab code for processing COVIS data.  See [`master_program/README.md`](README) in that directory for instructions on running that code.

`Test/` contains a Matlab-language test suite.    On a (Linux or Mac?) machine with Matlab and Gnu make installed:

  * `make test` or `make unittest` will run a (relatively fast) unit test suite.
  * `make integrationtest` will run a more thorough integration test suite.  This will
        fully process multiple COVIS files and may take some time to run.

The `Deploy/` directory contains code related to packaging this code as a Python library for distribution.   See the [`Deploy/README.md`](README) in that directory for further details.

This repo includes 
