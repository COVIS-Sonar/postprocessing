help:
	@echo "make unittest (or make test)        Runs Matlab-based unit test suite (from command line) -- runs quickly"
	@echo "make integrationtest                Runs Matlab-based integration test suite (from command line) -- may run slowly"

## By default run the shorter unit test
test: unittest

ALL_PATHS='..','../../Imaging','../../Diffuse','../../Doppler','../../Common'

unittest: covis_test_data
	git tag -d test && git tag test
	cd Test/Unit/ && matlab -nodisplay -nosplash -r "  addpath(${ALL_PATHS}); result = runtests(); disp(result); exit()"

integrationtest: covis_test_data
	cd Test/Integration/ && matlab -nodisplay -nosplash -r "  addpath(${ALL_PATHS}); result = runtests(); disp(result); exit()"

deps: gsw

gsw:
	mkdir -p master_program/ThirdParty/GSW/
	cd master_program/ThirdParty/GSW/ && \
		wget -nc http://www.teos-10.org/software/gsw_matlab_v3_06_11.zip && \
		unzip -f gsw_matlab_v3_06_11.zip

## Rule to retrieve git test data if it doesn't exist
covis_test_data: covis-test-data/old_covis_nas1.txt
	git submodule init && git submodule update


.PHONY: test unittest integrationtest covis_test_data
