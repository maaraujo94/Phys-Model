//////////////////////////////////////////////////////////////////////////////////////////
// code using a class to save everything together and not have to double-input anything //
//////////////////////////////////////////////////////////////////////////////////////////

#import "fitModel.C"

// main function
// argument determines whether we fit or open fit.txt and plot from there
// called with > root -l 'c_test.C($tofit)'
void c_test(bool tofit)
{
  // only input element given in main: lists of datasets and of parameters for the fit, fit and plot instructions
  const char* input_list = "state_list.txt";
  const char* param_list = "param_list.txt";
  const char* method_list = "param_method.txt";
  const char* plot_div = "plot_states.txt";
  
  // generates a map btw masses and state names
  fitclass.init();

  // note that the FitModel class has been initialized in fitModel.C and only has to be called here
  fitclass.getDataList(input_list);  
  fitclass.getData();

  fitclass.getParams(param_list);
  fitclass.getParMethod(method_list);

  // either do the fit or just read fit results
  if(tofit)  fitclass.doFit();
  else fitclass.readRes();
  
  fitclass.doPlot(plot_div);

  cout << "Fit option given was " << tofit << ", therefore";
  if(tofit) cout << " fit was performed" << endl;
  else cout << " previous fit results were read for the plotting" << endl;

  // testing normalization attribution
  int nstates = fitclass.states;
  double norm;
  vector <double> fitpar = fitclass.fitpar;
  int counter = 0;
  for(int id = 0; id < nstates; id++) {
    norm = 1;
    for(int i = 0; i < fitclass.norm_att[id].size(); i++)
      norm *= fitpar[fitclass.norm_att[id][i]];
    cout << id << " " << fitclass.auxnames[1][id] << " " << fitclass.auxnames[0][id] << ": ybin " << fitclass.datasigma[1][counter] << "-" << fitclass.datasigma[2][counter] << ", norm = " << norm << endl;
    counter+=fitclass.ns[id];
  }
  
  
}
