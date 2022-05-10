# Phys-Model
Global fit model using a physical description including a PDF. To run fit must use cycle.sh, to read the LHAPDF library.

Important files:

1) c_test
Does the actual fit. First calls all auxiliary functions to read fit data and fit parameters, then either fits and plots or just plots from existing fit file. More detailed info on auxiliary .txt files and the steps of the process in comments here

2) fitModel
Defines the class containing all necessary variables and methods for the fit workflow

3) aux_func
Contains auxiliary functions for fitting and plotting which don't depend on the full datasigma vector

4) cycle.sh
runs c_test(doFit), with doFit = 1 for fitting and doFit = 0 for just plotting (taken as input), then compiles latex file of plots and saves fit results in a folder

THINGS TO CHANGE
- shape part of the fit function still not very flexible. Devote some time to considering how to further improve this, might not be possible
- replace PTMNORM by a method that more accurately gives first prediction for L_QQ values
