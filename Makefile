help:
	@echo "make unittest                Runs Matlab-based unit test suite (from command line) -- runs quickly"
	@echo "make integrationtest         Runs Matlab-based integration test suite (from command line) -- may run slowly"
	@echo "make test                    Run both of the above tests"

## By default run the shorter unit test
test: unittest integrationtest

ALL_PATHS='..','../../Imaging','../../Diffuse','../../Doppler','../../Common'

unittest: covis_test_data
	git tag -d test && git tag test
	cd Test/Unit/ && matlab -nodisplay -nosplash -r "  addpath(${ALL_PATHS}); result = runtests(); disp(result); exit()"


## By default, every test will use a new temporary directory, so there will be
# a lot of time spent unpacking and reprocessing.  You can cirumvent this
# behavior by setting the COVIS_TEST_TEMP, which will be used as the tempdir
# location in every test -- note this will cause some tests to see
# existing data from previous tests and perhaps run differently
integrationtest: covis_test_data
	cd Test/Integration/ && matlab -nodisplay -nosplash -r "  addpath(${ALL_PATHS}); result = runtests(); disp(result); exit()"



## Rule to retrieve git test data if it doesn't exist
covis_test_data: covis-test-data/old_covis_nas1.txt
	git submodule init && git submodule update


.PHONY: test unittest integrationtest covis_test_data
