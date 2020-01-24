FROM amarburg/matlab-runtime:2019b

RUN apt-get update && apt-get install -y python3 python3-pip libgl1-mesa-dev patool \
                              libxt-dev p7zip-full netcat libarchive-dev  && \
          rm -rf /var/lib/apt/lists/*

ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libexpat.so

## Switch to non-root user "covis"
RUN groupadd -g 999 covis && \
    useradd -r -u 999 -g covis covis

WORKDIR /home/covis
RUN chown -R covis:covis /home/covis

USER covis

## Pre-stage some dependencies we know we need, so they get cached...
RUN pip3 install pytest python-decouple pyunpack patool minio

##
## Install pycovis.matlab
##
COPY --chown=covis:covis pycovis-matlab /home/covis/pycovis-matlab
WORKDIR /home/covis/pycovis-matlab
RUN pip3 install -e .

##
## Install pycovis.postprocess
##
COPY --chown=covis:covis pycovis-postprocess /home/covis/pycovis-postprocess
WORKDIR /home/covis/pycovis-postprocess
RUN pip3 install -e .

WORKDIR /home/covis

RUN  ln -s pycovis-matlab/Inputs .

VOLUME  /input
VOLUME  /output

ENV LD_LIBRARY_PATH=$MATLAB_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

ENTRYPOINT ["python3"]
