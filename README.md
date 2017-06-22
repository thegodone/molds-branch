# molds-branch

original copy v0.4 is located in SVN https://osdn.jp/projects/molds/releases/

MolDS ("Mol"ecular "D"ynamics simulation package with "S"emiempirical quantum chemistry) Under developement
Developers: Mikiya Fujii, Ph.D.(project lead), Katsuhiko Nishimra, and Michihiro Okuyama, Ph.D..
Other contributors: Michael Banck, Guillaume Godin

I add new atoms parameters : I, Br, P and we add RM1 support

To simplify compilation step all is available via a Docker file is also added in order to simplify the compilation step:


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


# generator of input file
You need rdkit to convert a smile to a input file using the following code
rdkit2molds.py

# advanced instructions  
doc/readme.md file
