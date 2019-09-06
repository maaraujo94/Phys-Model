# Phys-Model
Global fit model using a physical description including a PDF

Important files:
0) dataid_ref
consult to know which datasets are used in each fit, respective dataid and corrections to apply on cross-section measurements

1) c_test
Does the actual fit. First calls all auxiliary functions to read fit data and fit parameters, then either fits and plots or just plots from existing fit file

2) fitModel
Defines the class containing all necessary variables and methods for the fit workflow

3) aux_func_temp
Contains auxiliary functions for fitting and plotting which don't depend on the full datasigma vector

WORKFLOW: run c_test(doFit), with doFit = 1 for fitting and doFit = 0 for just plotting

THINGS TO CHANGE
- shape part of the fit function still not very flexible. Devote some time to considering how to further improve this, might not be possible