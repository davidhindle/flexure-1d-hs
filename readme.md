## flex-hs

[![DOI](https://zenodo.org/badge/351581782.svg)](https://zenodo.org/badge/latestdoi/351581782)

## Introduction

flex-hs is a fortran 90 code, consisting of a main code, and 2 modules, which implements the half station discretisation of the 1 dimensional, flexure equation with a variable coefficient. For a comprehensive explanation of the maths behind this model and why it is necessary, please see Hindle and Besson, 2021.
 

## Important

This is a "bare bones" repository to get the source code for the paper, Hindle and Besson (2021) published and available for download. These instructions are therefore not failsafe. If you are not experienced in Fortran 90 usage and Linux shell environments, things may go wrong and you may not understand why. Moroever, it is important to note that this program has been developed, compiled and run using the GNU Fortran compiler gfortran on a Linux system. The code employs some rather specific tricks, for instance 128 bit quadruple precision variables, and it cannot be guaranteed to work with other Fortran compilers and on other systems. I would thus welcome reports (up to the point of saturation) from any users on other platforms who find problems, and especially from users who find AND SOLVE those problems.

## Requirements

Ideally <br/>
1) Linux or BSD or similar UNIX type environment <br/>
2) GNU fortran (gfortran) compiler <br/>
3) GMTv6.x installed (for visualisation) <br/><br/>
If you have these, you are pretty much okay. If you have a Windows machine, it may be possible to use the power shell and the g95 compiler to get things to work, but there are no guarantees and this has not been tested.

## Steps to running the code.

Download the main code and modules.<br/>
flex-hs.f90 (main code)<br/>
module_pdma.f90 (pentadiagonal matrix solver)<br/>
module_plotflex.f90 (module containing plotting subroutines)<br/>
Makefile<br/>
place all 4 files in a directory of your choice, best to create fresh one for the purpose.<br/>
Open a terminal window and navigate to the directory.<br/><br/>
Type make<br/>
You should now have an exectable file flex-hs.exe. If you are an experienced Linux user, you will probably place the executable somewhere in your path, and make it globally available on your system. You need no help from me here.<br/><br/>
For less experienced users however, the simple way to run the program is to type
`./flex-hs.exe`

## required input files 

Compiling the code is not everything. The code is structured to read very specific input files with very specific formats, with very specific meanings. Some idea of these can be got from reading the comments within the main source code, and some idea can be got from the paper, Hindle and Besson 2021. The main input file, param.txt also has a "header" line which is supposed to explain each variable on the line below.Nevertheless, this is not ideal.<br/><br/>
The input files for the code are<br/><br/>
`param.txt` - a file that has grown and grown in size that now contains a bewildering set of inputs, some for controlling the code's default behaviours, and some simply giving the values of parameters for the problem.<br/>
the variables in `param.txt` are, in order, on the third line of the file<br/><br/>
**EM** - elastic modulus, Pa; (6.5e10 for continental lithosphere) <br/>
**h **- "default" elastic thickness, metres;  all nodes are set to this value by default. Subsequently, the value will be changed only at those nodes specified within te.txt <br/>
**pm** - mantle density, kgm^3; <br/>
**pload **- load density, kgm^3; density of user defined load given as a load thickness per node, in load-flex.txt<br/>
**pfill/pw **- basin fill density, kgm^3; density of material filling basins created by flexure. Setting to 1.d3 is equivalent to water filled basins. Setting to 0, means empty basins. <br/>
**pcrust** - crustal density, kgm^3; density of crust, effectively the density of material eroded from the plate if erosion is switched on.<br/>
**g** - gravitation, -9.81 ms^-2; THE CODE IS SET TO RUN WITH GRAVITY NEGATIVE!!! DO NOT CHANGE THIS VALUE.<br/>
**dx** - grid spacing, metres; space between nodes in the model. Effectively determines size of problem domain<br/>
**v** - poisson's ratio, generally taken to be 0.25;<br/>
**pu-horiz** - horizontal plate wide, background stress, Nm^-1; a messy parameter to use when elastic thickness varies. Often better to set to 0<br/>
**iter tol **- tolerance for ending iteration; the absolute difference between succeeding iterations of the value of u(x). When the iteration converges, usually, the change in value between iteration is very small<br/>
**nodes** - number of nodes - INTEGER; number of nodes must be an integer. Together with grid spacing, this determines physical dimension of problem domain.<br/>
**iswitch** - iteration on or off - LOGICAL; true = iteration is switched on, false = iteration is switched off, and solution is returned prior to iteration with all basins empty<br/>
**pint** - plot interval, INTEGER; the integer number of the increment of nodes between each recorded value of u. Setting to 5 means only every 5th node value of u is recorded.<br/>
**bc3rd **- boundary condition, LOGICAL; true = 3rd derivative and 2nd derivative boundary conditions are applied at left hand side of model. <br/>
**comp**- erosion compensation for +ve u, LOGICAL; true = erosion of all areas of +ve u, generating upward force pcrust*g*u <br/><br/>
`te.txt` - set elastic thickness values at specific nodes. For instance, setting te values from nodes 24999 t0 25001 0.01m, and nodes 25600 25610 to 10m requires two lines in te.txt<br/><br/>
`24999 25001 0.01 `<br/>
`25600 25610 10 ` <br/><br/>
Every node is initialised in the code to have elastic thickness = **h**, taken from param.txt. Hence, if te.txt is empty, the model will have a constant elastic thickness everywhere, corresponding to the "background" value from param.txt<br/>
You can also use pwl_interp.f90 to automatically generate a piecewise, linear interpolated te.txt, where you give interpolation points and values of te at those points only, and the intervening nodes are linearly interpolated. The output from pwl_interp.f90 is called te.txt by default, and will be used by flex-hs automatically if it has been generated.<br/><br/>
`load-flex.txt` - set elastic thickness values at specific nodes. For instance, setting load values from nodes 25010 to 25110 at 3000m, requires one line in load-flex.txt<br/><br/>
`25010 25110 3000` <br/>
`<br/><br/>


## output and plots
Output is meant to be visualised directly by GMT v6, with the script plotflex6.gmt. The important output file is hs_udqp.txt which is the solution u(x) at all nodes, plotted with a kms x scale and u(x) in metres. 



## Installation

* Download the zip file or [clone](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository) the repository.
* Unzip the file, navigate to the directory and run the minimal example or one of the [Jupyter](https://jupyter.org/) notebooks.






## Authors
**David Hindle**, University of Göttingen, <dhindle-at-gwdg.de>


## Reference

This notebook has been published at Zenodo. Please cite the following reference if you publish any work that uses this notebook:

Hindle, D. (2021). flexure-1d-hs Zenodo. [https://doi.org/10.5281/zenodo.4642172](https://doi.org/10.5281/zenodo.4642172)

[![DOI](https://zenodo.org/badge/351581782.svg)](https://zenodo.org/badge/latestdoi/351581782)


## License
This project is licensed under the GNU lesser general public license (LGPL v3). See the [LICENSE.txt](LICENSE.txt) file for details.


