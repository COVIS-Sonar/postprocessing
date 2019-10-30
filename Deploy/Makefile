

IMAGE_TAG = amarburg/covis-postprocess:latest


help:
	@echo "make matlab          Use Matlab compiler to build Python package 'pycovis-postprocess'"
	@echo
	@echo "== Python tasks run in the Docker image =="
	@echo "make docker           Build the docker image"
	@echo "make force_docker     Build the docker image with --no-cache"
	@echo "make push            Push the docker image to DockerHub"
	@echo "make test            Run pytest in the docker image"
	@echo "make version         Demo which runs scripts/version.py in the Docker image"
	@echo
	@echo "== Python tasks run locally (requires Matlab or Matlab runtime) =="
	@echo
	@echo "make local_pytest"
	@echo "make local_version"
	@echo "make local_test_process"



#== Tasks related to building and testing Docker scripts ==
docker:  #pycovis-matlab/setup.py pycovis-postprocess/pycovis/postprocess/*.py
	docker build --tag ${IMAGE_TAG} .

force_docker: pycovis-matlab/setup.py
	docker build --tag ${IMAGE_TAG} --no-cache .

push: build
	docker push ${IMAGE_TAG}

pytest: build
	mkdir -p pytest_output/
	docker run --rm \
						-v $(shell pwd)/../TestData/:/input \
						-v $(shell pwd)/pytest_output:/output \
						${IMAGE_TAG} --workdir pycovis-postprocess -m pytest


#== Some examples of running python scripts in the Docker image... ==
version: build
	docker run -i --rm ${IMAGE} covis_version.py

imaging: build
	mkdir -p output/
	docker run --rm \
					-v $(shell pwd)/../TestData/:/input \
					-v $(shell pwd)/pytest_output:/output \
 					${IMAGE} covis_process.py

##== Tasks which run the Python library locally ==

LOCAL_MATLAB_INSTALL ?= /usr/local/MATLAB/R2019b
MATLAB_LD_LIBRARY_PATH=$(LOCAL_MATLAB_INSTALL)/runtime/glnxa64:$(LOCAL_MATLAB_INSTALL)/bin/glnxa64:$(LOCAL_MATLAB_INSTALL)/sys/os/glnxa64:$(LOCAL_MATLAB_INSTALL)/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:$(LOCAL_MATLAB_INSTALL)/sys/java/jre/glnxa64/jre/lib/amd64/server:$(LOCAL_MATLAB_INSTALL)/sys/java/jre/glnxa64/jre/lib/amd64

LOCAL_MATLAB_ENV = 	LD_LIBRARY_PATH=$(MATLAB_LD_LIBRARY_PATH):$(LD_LIBRARY_PATH) \
		PYTHONPATH=$(shell pwd)/pycovis-matlab:$(shell pwd)/pycovis-postprocess

local_pytest:
		${LOCAL_MATLAB_ENV} cd pycovis-postprocess && python3 -m pytest

local_version:
		${LOCAL_MATLAB_ENV} python3 covis_version.py

local_test_process:
		${LOCAL_MATLAB_ENV} python3 covis_process.py /input/COVIS-20191024T003346-diffuse3.tar.gz



##== Matlab packaging rules. ==
## Must be run on a machine with Matlab Compiler SDK installed

## TODO.  How to ensure the version of Matlab / mcc used matches the
## runtime in the image?

COVIS_REPO=..
MCC=mcc

matlab: pycovis-matlab/setup.py

GITREV=${shell git rev-parse HEAD }
GITTAG=${shell git describe --tags }

COVIS_CODE = $(shell find ${COVIS_REPO}/master_program/ -name *.m) \
					 	 $(shell find ${COVIS_REPO}/ThirdParty/ -name *.m)

tmp/static_git_info.m:
	mkdir -p tmp/
	echo "function [gitinfo] = static_git_info()" > $@
	printf "   gitinfo = struct;\n" >> $@
	printf "   gitinfo.gitrev = '%s';\n" ${GITREV} >> $@
	printf "   gitinfo.gittags = '%s';\n" ${GITTAG} >> $@
	printf "   gitinfo.dockertag = '%s';\n" ${IMAGE} >> $@


pycovis-matlab/setup.py: $(COVIS_CODE) tmp/static_git_info.m
	mkdir -p pycovis-matlab/Input
	${MCC} -v -d pycovis-matlab/ -W python:pycovis.matlab -T link:lib $^
	cp -r ${COVIS_REPO}/master_program/Inputs pycovis-matlab/


clean:
		rm -rf pycovis-matlab/ tmp/


.PHONY: build force_build test help matlab version imaging clean tmp/static_git_info.m
