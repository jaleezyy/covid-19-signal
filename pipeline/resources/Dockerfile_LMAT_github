FROM ubuntu:18.04

LABEL base.image="ubunut"
LABEL version="1"
LABEL software="LMAT"
LABEL software.version="1.2.4/28fdb44"
LABEL description="Livermore Metagenome Analysis Toolkit"
LABEL website="https://github.com/LivGen/LMAT"
LABEL documentation="https://rawgit.com/LivGen/LMAT/master/doc/lmat-doc/index.html"

# install deps
RUN apt-get update && \
        apt-get install -y git python3 python3-dev \
            zlib1g-dev build-essential cmake mpi

# build and install LMAT
RUN git clone https://github.com/LivGen/LMAT && \
        cd LMAT && \
        git checkout 28fdb44 && \
        ./redoall.sh R gnu && \
        chmod +x bin/* && \
        cp bin/* /bin

WORKDIR /data/
ENV LMAT_DIR /data
