

IMAGE = amarburg/covis-postprocess:R2018a

help:
	@echo "make matlab          Use Matlab compiler to build Python package 'pycovis-postprocess'"
	@echo "make build           Build the docker image"
	@echo "make force_build     Build the docker image with --no-cache"
	@echo "make push            Push the docker image to DockerHub"
	@echo "make test            Run pytest in the docker image"
	@echo "make version         Demo which runs scripts/version.py in the Docker image"


#== Tasks related to building and testing Docker scripts ==
build:       pycovis-matlab/setup.py pycovis-postprocess/pycovis/postprocess/*.py
	docker build --tag ${IMAGE} .

force_build: pycovis-matlab/setup.py
	docker build --tag ${IMAGE} --no-cache .

push: build
	docker push ${IMAGE}

test: build
	mkdir -p output/
	docker run --rm \
						-v $(shell pwd)/../covis-test-data/old-covis-nas1/raw/2011/10/01/:/input \
						-v $(shell pwd)/output:/output \
						${IMAGE} -m pytest --capture=no




#== Some examples of running python scripts in the Docker image... ==
version: build
	docker run -i --rm ${IMAGE} scripts/version.py

imaging: build
	mkdir -p output/
	docker run --rm \
					-v $(shell pwd)/../covis-test-data/old-covis-nas1/raw/2011/10/01/:/input \
					-v $(shell pwd)/output:/output \
 					${IMAGE} scripts/process_sweep.py


##== Matlab packaging rules. ==
## Must be run on a machine with Matlab Compiler SDK installed

## TODO.  How to ensure the version of Matlab / mcc used matches the
## runtime in the image?

COVIS_REPO=..
COVIS_INPUT=${COVIS_REPO}/Common/input
MCC=mcc

matlab: pycovis-matlab/setup.py

GITREV=${shell git rev-parse HEAD }
GITTAG=${shell git describe --tags --exact-match }

## How write a Makefile rule which will only rebuild if the Git commit has changed?
tmp/static_git_info.m: ../.git/index
	mkdir -p tmp/
	echo "function [gitinfo] = static_git_info()" > $@
	printf "   gitinfo = struct;\n" >> $@
	printf "   gitinfo.gitrev = '%s';\n" ${GITREV} >> $@
	printf "   gitinfo.gittags = '%s';\n" ${GITTAG} >> $@
	printf "   gitinfo.dockertag = '%s';\n" ${IMAGE} >> $@

	cat $@


MFILES = ${COVIS_REPO}/Common/*.m \
					${COVIS_REPO}/Imaging/*.m \
					${COVIS_REPO}/Diffuse/*.m \
					${COVIS_REPO}/Doppler/*.m

pycovis-matlab/setup.py: ${MFILES} tmp/static_git_info.m
	mkdir -p pycovis-matlab/
	${MCC} -v -w enable -d pycovis-matlab/ -W python:pycovis.matlab -T link:lib $^
	cp -r $(COVIS_INPUT) pycovis-matlab/


clean:
		rm -rf pycovis-matlab/ tmp/


.PHONY: build force_build test help matlab version imaging clean
