help:
	@echo "make unittest (or make test)        Runs Matlab-based unit test suite (from command line) -- runs quickly"
	@echo "make integrationtest                Runs Matlab-based integration test suite (from command line) -- may run slowly"
	@echo "make deps													 Install necessary dependencies to ThirdParty/"
	@echo "make testdata  										 Run \'make download\" in TestData/"


## By default run the shorter unit test
test: unittest

## Well, this is unpleasant and ugly, isn't it
ALL_MATLAB_PATHS= ${shell cd Test/Unit/ && find ../../master_program -type d}
EMPTY :=
SPACE := $(EMPTY) $(EMPTY)
COMMA := ,


unittest: covis_test_data
	git tag -f test
	cd Test/Unit/ && matlab -nodisplay -nosplash -r "  addpath('..',${subst $(SPACE),$(COMMA),$(patsubst %,'%',$(ALL_MATLAB_PATHS))}); result = runtests(); disp(result); exit()"

integrationtest: covis_test_data
	cd Test/Integration/ && matlab -nodisplay -nosplash -r "  addpath('..',${subst $(SPACE),$(COMMA),$(ALL_MATLAB_PATHS)}); result = runtests(); disp(result); exit()"

##
deps: gsw

gsw:
	mkdir -p master_program/ThirdParty/GSW/
	cd master_program/ThirdParty/GSW && \
			wget http://www.teos-10.org/software/gsw_matlab_v3_06_11.zip && \
			unzip gsw_matlab_v3_06_11.zip && \
			rm gsw_matlab_v3_06_11.zip

## Rule to retrieve git test data if it doesn't exist
testdata:
	cd TestData && make download


.PHONY: test unittest integrationtest covis_test_data deps gsw testdata
