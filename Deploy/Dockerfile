FROM amarburg/matlab-runtime:2019b

RUN apt-get update && apt-get install -y python3 python3-pip libgl1-mesa-dev libxt-dev

ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libexpat.so

RUN pip3 install pytest

##
## Install pycovis.matlab
##
COPY pycovis-matlab /root/pycovis-matlab
WORKDIR /root/pycovis-matlab
RUN pip3 install -e .

##
## Install pycovis.postprocess
##
COPY pycovis-postprocess /root/pycovis-postprocess
WORKDIR /root/pycovis-postprocess
RUN pip3 install -e .

## Copy in the test functions
COPY test                /root/test
COPY scripts             /root/scripts

WORKDIR /root
RUN  ln -s /root/pycovis-matlab/Inputs .

VOLUME  /input
VOLUME  /output

ENV LD_LIBRARY_PATH=$MATLAB_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

ENTRYPOINT ["python3"]
