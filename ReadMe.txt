Intro:

This folder contains a C++ code for a single bifurcation model of pulse wave propagation. The example code uses 3-element Windkessel model that is reducible to a pure resistance condition. A generic inflow waveform is used to drive the system. A Makefile is used to link all the modules in the right order. Current version of the Makefile is updated to work on an iMac (Make sure to modify the compiler versions and/or path links specific to your machine).



Files Description

1. Makefile: Run it by typing "make" command in the command window to create an executable “sor06”. Google “Visual Studio For Windows” for help making on the Windows. There might be many other options available.
2. main.m: Run this file in the MATLAB After creating the executable. Change geometric and hemodynamic parameter in this film to create a new example

C++ Files

3. sor06.h:  Header file containing global parameters
4. sor06.C:  Defines the network and calls the solver and the prints the data in vessel specific files.
5. arteries.h: Header file declaring the class “Tube” and all its objects and functions, to be specified in arteries.C
6. arteries.C: Contains the main computational code, which uses Lax-Wenderoff scheme and a linear pressure-area relation. It also calls tools.C, tools.h (for root finding at bifurcations).
7. tools.C: Includes numerical algebraic tools for finding roots for the system of equations
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



