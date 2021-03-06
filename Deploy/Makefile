

IMAGE_TAG = amarburg/covis-postprocess:latest

ifdef http_proxy
	HTTP_PROXY_ARG := --build-arg http_proxy=$(http_proxy)
endif

help:
	@echo "make matlab          Use Matlab compiler to build Python package 'pycovis-postprocess'"
	@echo
	@echo "== Python tasks run in the Docker image =="
	@echo "make docker          Build the docker image"
	@echo "make force_docker    Build the docker image with --no-cache"
	@echo "make push            Push the docker image to DockerHub"
	@echo "make pytest          Run pytest in the docker image"
	@echo "make version         Demo which runs scripts/version.py in the Docker image"
	@echo
	@echo "make diffuse3.tar.gz_from_local"
	@echo "make diffuse3.7z_from_local"
	@echo "== Python tasks run locally (requires Matlab or Matlab runtime) =="
	@echo
	@echo "make local_pytest"
	@echo "make local_version"
	@echo "make local_test_process"


PP_SRC = $(shell find pycovis-postprocess -name "*.py" )

#== Tasks related to building and testing Docker scripts ==
docker:  ${PP_SRC}
	docker build ${HTTP_PROXY_ARG} --tag ${IMAGE_TAG} .

force_docker: ${PP_SRC}
	docker build ${HTTP_PROXY_ARG} --tag ${IMAGE_TAG} --no-cache .

push: build
	docker push ${IMAGE_TAG}

pytest: build
	mkdir -p tmp/pytest_output/
	docker run --rm \
						-v $(shell pwd)/../TestData/:/input \
						-v $(shell pwd)/tmp/pytest_output:/output \
						${IMAGE_TAG}  -m pytest

#						--workdir $(shell pwd)/pycovis-postprocess \

#== Some examples of running python scripts in the Docker image... ==
version: build
	docker run -i --rm ${IMAGE_TAG} pycovis-postprocess/scripts/version.py

TEST_DOCKER_MOUNTS = -v $(shell pwd)/../TestData/:/input \
											-v $(shell pwd)/output:/output \

diffuse3: diffuse3.tar.gz_from_local diffuse3.7z_from_local

diffuse3.tar.gz_from_local: build
	mkdir -p output/
	docker run --rm ${TEST_DOCKER_MOUNTS} \
 					${IMAGE_TAG} pycovis-postprocess/scripts/process.py /input/2019/10/24/COVIS-20191024T003346-diffuse3.tar.gz --output /output/

diffuse3.7z_from_local: build
	mkdir -p output/
	docker run --rm ${TEST_DOCKER_MOUNTS} \
 					${IMAGE_TAG} pycovis-postprocess/scripts/process.py /input/2019/10/24/COVIS-20191024T003346-diffuse3.7z --output /output/

imaging1: imaging1.tar.gz_from_local imaging1.7z_from_local

imaging1.tar.gz_from_local: build
	mkdir -p output/
	docker run --rm ${TEST_DOCKER_MOUNTS} \
 					${IMAGE_TAG} pycovis-postprocess/scripts/process.py /input/2019/10/24/COVIS-20191024T000002-imaging1.tar.gz --output /output/

imaging1.7z_from_local: build
	mkdir -p output/
	docker run --rm ${TEST_DOCKER_MOUNTS} \
 					${IMAGE_TAG} pycovis-postprocess/scripts/process.py /input/2019/10/24/COVIS-20191024T000002-imaging1.7z --output /output/


COVIS-20200421T133002-diffuse1.7z_from_local: build
	mkdir -p output/
	docker run --rm ${TEST_DOCKER_MOUNTS} \
 					${IMAGE_TAG} pycovis-postprocess/scripts/process.py --output /output/ /input/2020/04/21/COVIS-20200421T133002-diffuse1.7z



##== Tasks which run the Python library locally ==

LOCAL_MATLAB_INSTALL ?= /usr/local/MATLAB/R2020a
LOCAL_MATLAB_ENV = 	\
		PYTHONPATH=$(shell pwd)/pycovis-matlab:$(shell pwd)/pycovis-postprocess

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
	MATLAB_ARCH = glnx64a

	MATLAB_LD_LIBRARY_PATH=$(LOCAL_MATLAB_INSTALL)/runtime/$(MATLAB_ARCH):$(LOCAL_MATLAB_INSTALL)/bin/$(MATLAB_ARCH):$(LOCAL_MATLAB_INSTALL)/sys/os/$(MATLAB_ARCH):$(LOCAL_MATLAB_INSTALL)/sys/java/jre/$(MATLAB_ARCH)/jre/lib/amd64/native_threads:$(LOCAL_MATLAB_INSTALL)/sys/java/jre/$(MATLAB_ARCH)/jre/lib/amd64/server:$(LOCAL_MATLAB_INSTALL)/sys/java/jre/$(MATLAB_ARCH)/jre/lib/amd64
	LOCAL_MATLAB_ENV += LD_LIBRARY_PATH=$(MATLAB_LD_LIBRARY_PATH):$(LD_LIBRARY_PATH)
else ifeq ($(UNAME_S),Darwin)
## n.b.  Still have problems running locally on Mac ...  :-(
	MATLAB_ARCH = maci64

	MATLAB_DYLD_LIBRARY_PATH=$(LOCAL_MATLAB_INSTALL)/runtime/$(MATLAB_ARCH):$(LOCAL_MATLAB_INSTALL)/bin/$(MATLAB_ARCH):$(LOCAL_MATLAB_INSTALL)/sys/os/$(MATLAB_ARCH)
	LOCAL_MATLAB_ENV += DYLD_LIBRARY_PATH=$(MATLAB_DYLD_LIBRARY_PATH):$(DYLD_LIBRARY_PATH)
else
	$(error Unable to detect OS)
endif



local_pytest:
		${LOCAL_MATLAB_ENV} cd pycovis-postprocess && python3 -m pytest

local_version:
		${LOCAL_MATLAB_ENV} python3 pycovis-postprocess/scripts/version.py

local_diffuse3:
		${LOCAL_MATLAB_ENV} python3 pycovis-postprocess/scripts/process.py ../TestData/2019/10/04/COVIS-20191024T003346-diffuse3.tar.gz --output /output/COVIS-20191024T003346-diffuse3.mat



##== Matlab packaging rules. ==
## Must be run on a machine with Matlab Compiler SDK installed

## TODO.  How to ensure the version of Matlab / mcc used matches the
## runtime in the image?

COVIS_REPO=..
MCC=mcc

matlab: pycovis-matlab/setup.py ../master_program/ThirdParty/GSW/

GITREV=${shell git rev-parse HEAD }
GITTAG=${shell git describe --tags }


COVIS_CODE = $(shell find ${COVIS_REPO}/master_program/ -name *.m) \
						 $(shell find matlab/ -name *.m)


tmp/static_git_info.m:
	mkdir -p tmp/
	echo "function [gitinfo] = static_git_info()" > $@
	printf "   gitinfo = struct;\n" >> $@
	printf "   gitinfo.gitrev = '%s';\n" ${GITREV} >> $@
	printf "   gitinfo.gittags = '%s';\n" ${GITTAG} >> $@
	printf "   gitinfo.dockertag = '%s';\n" ${IMAGE} >> $@


pycovis-matlab/setup.py: $(COVIS_CODE) tmp/static_git_info.m
	${MCC} -v -d pycovis-matlab/ -W python:pycovis.matlab -T link:lib $^
	cp -r ${COVIS_REPO}/master_program/Inputs pycovis-matlab/


clean:
		rm -rf pycovis-matlab/ tmp/


.PHONY: build force_build pytest help matlab version imaging clean tmp/static_git_info.m
