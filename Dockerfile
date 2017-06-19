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


RUN git clone https://github.com/thegodone/OpenBLAS.git
RUN cd OpenBLAS
RUN make FC=gfortran
RUN make PREFIC=/usr/local/ install

RUN wget https://github.com/Reference-LAPACK/lapack/archive/v3.7.0.tar.gz
RUN wget http://www.netlib.org/blas/blast-forum/cblas.tgz

RUN tar -zxvf v3.7.0.tar.gz
RUN tar -C /cblas -zxvf cblas.gz

RUN git clone https://github.com/thegodone/molds-branch.git
RUN cd molds-branch/src
RUN cp Makefile-GNU Makefile
RUN make clean
