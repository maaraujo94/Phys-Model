/* calls on the members of a fitting class (in fitModel) to
 - read options
 - get input
 - perform fit (optional)
 - save results and plot fits */

#import "fitModel.C"

// main function
// argument determines whether we fit or simply open fit.txt and plot from there
// called with > root -l 'c_test.C($tofit)'
void c_test(int tofit)
{
  // only input element given in main: lists of datasets and of parameters for the fit, fit and plot instructions
  const char* input_list = "state_list.txt"; //input files for each state+detector+y_bin, uncertainty method and correction factors. use # to remove a state from the fit
  const char* param_list = "param_list.txt"; //fit parameters. set last value to 1 to fix in the fit. norm parameters have relative uncertainty, shape have central value and uncertainty, nuis have uncertainty
  const char* method_list = "param_method.txt"; //sets which norm and nuisance parameters affect each dataset. Shouldn't need changes
  const char* plot_div = "plot_states.txt"; //sets which datasets are plotted together. Shouldn't need changes
  
  // generates a map btw masses and state names
  fitclass.init();

  // note that the FitModel class has been initialized in fitModel.C and only has to be called here
  fitclass.getDataList(input_list);  
  fitclass.getData();

  fitclass.getParams(param_list);
  fitclass.getParMethod(method_list);

  // either do the fit or just read fit results
  if(tofit == 1)  fitclass.doFit();
  else fitclass.readRes();

  // plot fit results
  fitclass.doPlot(plot_div);
  
  cout << "Fit option given was " << tofit << ", therefore";
  if(tofit == 1) cout << " fit was performed" << endl;
  else if(tofit == 0) cout << " previous fit results were read for the plotting" << endl;

  //delete pdf_ct;
}
