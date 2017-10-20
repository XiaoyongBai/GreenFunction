# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* This is a c++ repository for Green's function of elasto-statics and elasto-dyanmics. 

### How to compile the fortran library of multi-layer Green's function

* step 1: gfortran -c ./kinds.f90 ./haldgreen4.f90
* step 2: ar cr libhald.a haldgree4.o kinds.o
* The "libhald.a" is the library we want. To used functions defined in the fortran codes, it should be linked to your project. 
* In this project, it is linked to the GreenFunction.
* !!NOTE: The runtime library of fortran, namely, "libgfortran.a" must be linked as well to enalbe fortran function in cpp code.

### How to generate a Xcode project for MacOS? ###
* cmake -G Xcode ..


### Contribution guidelines ###

* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact