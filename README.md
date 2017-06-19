# molds-branch

This is a clone of the original projet in subversion:

https://osdn.jp/projects/molds/releases/

the goal is to add more atoms in the computation

Docker file is also added in order to simplify the compilation step:


1.install Docker on your computer from there https://www.docker.com/

2. copy the Dockerfile in this repository in a folder on your computer name "Dockers" for example

3. cd Dockers

4. docker build -t molds .

5. docker run -i -t molds /bin/bash

6. int the docker image you can do:
  cd molds-branch/src &&
  ./molds ../test/ch3I_am1_geo.in
  
7. close the docker image by doing:
  exit



