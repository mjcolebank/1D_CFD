Copyright (2019) M. J. Colebank, M. U. Qureshi, and M. S. Olufsen

Intro:

This folder contains a C++ code for a single bifurcation model of pulse wave propagation. 
The example code uses 3-element Windkessel model that is reducible to a pure resistance 
condition. A generic inflow waveform is used to drive the system. A Makefile is used to 
link all the modules in the right order. Current version of the Makefile is updated to 
work on an iMac (Make sure to modify the compiler versions and/or path links specific to 
your machine).



Files Description

1. Makefile: Creates all executable files for sor06, arteries, and tools files.
Makefile is originally written for Mac compiler, but can be adapted for Windows. 
We recommend using Visual Studio for help making on Windows platforms. 
There might be many other options available as well.

2. run_1D: This file can be run AS IS and will call the Makefile and subsequent
algorithms needed. Change geometric and hemodynamic parameters in this file
for each new example.

C++ Files

3. sor06.h:  Header file containing global parameters
4. sor06.c:  Defines the network and calls the solver and the prints the data in vessel specific files.
5. arteries.h: Header file declaring the class Tubeù and all its objects and functions, to be specified in arteries.C
6. arteries.c: Contains the main computational code, which uses Lax-Wenderoff scheme and a linear pressure-area relation. It also calls tools.C, tools.h (for root finding at bifurcations).
7. tools.c: Includes numerical algebraic tools for finding roots for the system of equations
8. tools.h: Used by tools.C and arteries.C

Input data

9. Qin_8192.dat: Generic inflow data created from InFlow.m

Other MATLAB files

10. run_1D.m: Plots the simulated data.
11. gnuplot.m: Rearrange the block matrix output
12. InFLow.m: Generates an inflow profile and saves it as Qin_8192.dat.
13. get_nominal_WK.m: Computes and non-dimensionalize vessel stiffness and the Windkessel parameters.

Output files

Code will generate a new data files with extension *.2d



