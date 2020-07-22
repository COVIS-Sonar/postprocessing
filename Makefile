help:
	@echo "make unittest (or make test)        Runs Matlab-based unit test suite (from command line) -- runs quickly"
	@echo "make integrationtest                Runs Matlab-based integration test suite (from command line) -- may run slowly"
	@echo "make deps                           Install necessary dependencies to ThirdParty/"
	@echo "make testdata                       Run 'make download' in TestData/"
	@echo "make fullimaging1                   Process the \"fullimaging1\" test data set(s)"

process:
	matlab -nodisplay -nosplash -r "initial_matlab_process()"


## By default run the shorter unit test
test: unittest

## Well, this is unpleasant and ugly, isn't it
ALL_MATLAB_PATHS= ${shell cd Test/Unit/ && find ../../master_program -type d | grep -v .git }
EMPTY :=
SPACE := $(EMPTY) $(EMPTY)
COMMA := ,


unittest: covis_test_data
	## postproc_metadata test assumes this tag has been set!
	git tag -f test
	cd Test/Unit/ && matlab -nodisplay -nosplash -r "  addpath('..',${subst $(SPACE),$(COMMA),$(patsubst %,'%',$(ALL_MATLAB_PATHS))}); result = runtests(); disp(result); exit()"

integrationtest: covis_test_data
	cd Test/Integration/ && matlab -nodisplay -nosplash -r "  addpath('..',${subst $(SPACE),$(COMMA),$(ALL_MATLAB_PATHS)}); result = runtests(); disp(result); exit()"


## Run on test datasets
FULLIMAGING1 = 'TestData/2020/07/19/COVIS-20200719T000001-fullimaging1/'
fullimaging1:
	if [ ! -d $(FULLIMAGING1) ]; then cd TestData/2020/07/19/ && 7z x COVIS-20200719T000001-fullimaging1.7z; fi
	mkdir -p tmp/fullimaging1/
	cd tmp/fullimaging1 && matlab -nodisplay -nosplash -r "addpath('..',${subst $(SPACE),$(COMMA),$(patsubst %,'%',$(ALL_MATLAB_PATHS))}); result = covis_raw_sweep('../../TestData/2020/07/19/','COVIS-20200719T000001-fullimaging1',1); disp(result); exit()"




#== Handle thirdparty dependencies
deps: gsw matlab-axis-label-alignment

matlab-axis-label-alignment:
	cd master_program/ThirdParty && git clone https://github.com/phymhan/matlab-axis-label-alignment.git

gsw:
	mkdir -p master_program/ThirdParty/GSW/
	cd master_program/ThirdParty/GSW && \
			wget http://www.teos-10.org/software/gsw_matlab_v3_06_11.zip && \
			unzip gsw_matlab_v3_06_11.zip && \
			rm gsw_matlab_v3_06_11.zip

## Rule to retrieve git test data if it doesn't exist
testdata:
	cd TestData && make download


.PHONY: test unittest integrationtest covis_test_data deps \
				gsw testdata process \
				fullimaging1
