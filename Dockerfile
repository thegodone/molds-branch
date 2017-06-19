FROM ubuntu:17.04

RUN apt-get update && apt-get install -y \
 mpi-default-dev \
 libboost-mpi-dev \
 libopenblas-dev \
 liblapacke-dev \
 liblapacke-dev \
 libboost-mpi-dev \
 libboost-serialization-dev \
 libboost-thread-dev \
 libboost-system-dev \
 make \
 git 

RUN git clone https://github.com/thegodone/molds-branch.git
RUN cd molds-branch/src
RUN cp Makefile-GNU Makefile
RUN make clean
