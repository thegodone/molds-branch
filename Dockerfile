FROM ubuntu:16.04

# install dependencies for the program (boost, mpi, g++, gfortran) and to download the archives and source code (make,)
RUN apt-get -q update
RUN apt-get install -qy \
 gfortran g++ libopenmpi-dev openmpi-bin make wget \
 libboost-mpi-dev  libboost-serialization-dev  libboost-thread-dev libboost-system-dev git vim

# this command work to build openblas with lapack header & include files!
RUN wget https://github.com/xianyi/OpenBLAS/archive/v0.2.19.tar.gz
RUN tar -zxvf v0.2.19.tar.gz
# compile OpenBLAS & LAPACK
RUN cd OpenBLAS-0.2.19 && make all BINARY=64 CC=/usr/bin/gcc FC=/usr/bin/gfortran USE_THREAD=0 INTERFACE64=1 1> make.log 2>make.err && make PREFIX=/usr/local/openblas install
# compile molds
RUN git clone https://github.com/thegodone/molds-branch.git
RUN cd molds-branch/src && cp Makefile_GNU Makefile && make clean && make -j8 && ./molds ../test/ch3Cl_am1_geo.in

