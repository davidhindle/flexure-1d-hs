# flex-hs

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4571576.svg)](https://doi.org/10.5281/zenodo.4571576)

## Introduction

flex-hs is a fortran 90 code, consisting of a main code, and 2 modules, which implements the half station discretisation of the 1 dimensional, flexure equation with a variable coefficient. For a comprehensive explanation of the maths behind this model and why it is necessary, please see Hindle and Besson, 2021.
 

## Getting started

This is a "bare bones" repository to get the source code for the paper, Hindle and Besson (2021) published and available for download. These instructions are therefore not failsafe. If you are not experienced in
Fortran 90 usage and Linux shell environments, things may go wrong and you may not understand why. Moroever, it is important to note that this program has been developed, compiled and run using the GNU Fortran
compiler gfortran on a Linux system. The code employs some rather specific tricks, for instance 128 bit quadruple precision variables, and it cannot be guaranteed to work with other Fortran compilers and on 
other systems. I would thus welcome reports (up to the point of saturation) from any users on other platforms who find problems, and especially from users who find AND SOLVE those problems.
Steps to running the code
Download the main code and modules. 
flex-hs.f90 (main code)
module_pdma.f90 (pentadiagonal matrix solver)
module_plotflex.f90 (module containing plotting subroutines)

place all 3 files in a directory of your choice, best to create fresh one for the purpose.
Open a terminal window and navigate to the directory.
type the commands
gfortran -c module_pdma.f90
gfortran -c module_plotflex.f90
gfortran flex-hs.f90 module_pdma.f90 module_plotflex.f90 -o  hsflex (or any other name you want to give the executable).

you have now created the executable, hsflex. If you are an experienced Linux user, you will probably place the executable somewhere in your path, and make it globally available on your system. You need no elp
from me here.
For less experienced users however, the simple way to run the program is to type
./hsflex (or whatever else you called the program). 

### required input files 

Compiling the code is not everything. The code is structured to read very specific input files with very specific formats, with very specific meanings. Some idea of these can be got from reading the comments 
within the main source code, and some idea can be got from the paper, Hindle and Besson 2021. The main input file, param.txt also has a "header" line which is supposed to explain each variable on the line below.
Nevertheless, this is not ideal. 

The input files for the code are 
param.txt - a file that has grown and grown in size that now contains a bewildering set of inputs, some for controlling the code's default behaviours, and some simply giving the values of parameters for the problem. 
the variables in param.txt are, in order, on the third line of the file
EM - elastic modulus, Pa; (6.5e10 for continental lithosphere) 
h - "default" elastic thickness, metres;  all nodes are set to this value by default. Subsequently, the value will be changed only at those nodes specified within te.txt 
pm - mantle density, kgm^3; 
pload - load density, kgm^3; density of user defined load given as a load thickness per node, in load-flex.txt
pfill/pw - basin fill density, kgm^3; density of material filling basins created by flexure. Setting to 1.d3 is equivalent to water filled basins. Setting to 0, means empty basins. 
pcrust - crustal density, kgm^3; density of crust, effectively the density of material eroded from the plate if erosion is switched on.
g - gravitation, -9.81 ms^-2; THE CODE IS SET TO RUN WITH GRAVITY NEGATIVE!!! DO NOT CHANGE THIS VALUE.
dx - grid spacing, metres; space between nodes in the model. Effectively determines size of problem domain
v - poisson's ratio, generally taken to be 0.25;
pu-horiz - horizontal plate wide, background stress, Nm^-1; a messy parameter to use when elastic thickness varies. Often better to set to 0
iter tol - tolerance for ending iteration; the absolute difference between succeeding iterations of the value of u(x). When the iteration converges, usually, the change in value between iteration is very small
nodes - number of nodes - INTEGER; number of nodes must be an integer. Together with grid spacing, this determines physical dimension of problem domain.
iswitch - iteration on or off - LOGICAL; true = iteration is switched on, false = iteration is switched off, and solution is returned prior to iteration with all basins empty
pint - plot interval, INTEGER; the integer number of the increment of nodes between each recorded value of u. Setting to 5 means only every 5th node value of u is recorded.
bc3rd - boundary condition, LOGICAL; true = 3rd derivative and 2nd derivative boundary conditions are applied at left hand side of model. 
 
### or by running a minimal example:



## Installation

* Download the zip file or [clone](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository) the repository.
* Unzip the file, navigate to the directory and run the minimal example or one of the [Jupyter](https://jupyter.org/) notebooks.


## Required Python modules



## Authors
**David Hindle**, University of Göttingen, <dhindle-at-gwdg.de>


## Reference

This notebook has been published at Zenodo. Please cite the following reference if you publish any work that uses this notebook:

Hindle, D. (2021). flexure-1d-hs Zenodo. [https://doi.org/10.5281/zenodo.4571576](https://doi.org/10.5281/zenodo.4571576)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4571576.svg)](https://doi.org/10.5281/zenodo.4571576)


## License
This project is licensed under the GNU lesser general public license (LGPL v3). See the [LICENSE.txt](LICENSE.txt) file for details.


