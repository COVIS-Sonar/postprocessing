help:
	@echo "make unittest         Runs Matlab-based unit test suite (from command line) -- runs quickly"
	@echo "make integrationtest         Runs Matlab-based integration test suite (from command line) -- may run slowly"


test: unittest

unittest:
	cd Test/Unit/ && matlab -nodisplay -nosplash -r "  addpath('..','../../Common'); result = runtests(); disp(result); exit()"

integrationtest:
	cd Test/Integration/ && matlab -nodisplay -nosplash -r "  addpath('..','../../Common'); result = runtests(); disp(result); exit()"


.PHONY: test unittest integrationtest
