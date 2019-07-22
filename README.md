# Phys-Model
Global fit model using a physical description including a PDF

Important files:
0) dataid_ref
consult to know which datasets are used in each fit, respective dataid and corrections to apply on cross-section measurements

1) main_fit
Does the actual fit. Fit parameters initialized here

2) main_plot
Does the plotting of the results and pulls, from the fit results file fit.txt

3) aux_func
Contains all auxiliary functions for fitting and plotting

4) aux_input
Reads in all the data to the datasigma vector. Data to be read defined here

WORKFLOW: run main_fit, then main_plot

THINGS TO CHANGE
- Input done with a fixed stateid code that I input and in turn receive the corresponding states -> also determines which parts of the fit code I use
- Make list of stateid values to use an input parameter? 