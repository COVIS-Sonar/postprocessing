FROM amarburg/matlab-runtime:latest

RUN apt-get update && apt-get install -y python3 python3-pip libgl1-mesa-dev libxt-dev

ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libexpat.so

RUN pip3 install pytest pytest-faulthandler

##
## Install pycovis.matlab
##
COPY pycovis-matlab /code/pycovis-matlab
WORKDIR /code/pycovis-matlab
RUN pip3 install -e .

##
## Install pycovis.postprocess
##
COPY pycovis-postprocess /code/pycovis-postprocess
WORKDIR /code/pycovis-postprocess
RUN pip3 install .

## Copy in the test functions
COPY test                /root/test
COPY scripts             /root/scripts

WORKDIR /code
RUN  ln -s /code/pycovis-matlab/input .

VOLUME  /input
VOLUME  /output

ENV LD_LIBRARY_PATH=$MATLAB_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
ENV COVIS_IN_DOCKER=1

ENTRYPOINT ["python3"]
