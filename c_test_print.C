//////////////////////////////////////////////////////////////////////////////////////////
// code using a class to save everything together and not have to double-input anything //
//////////////////////////////////////////////////////////////////////////////////////////

#import "fitModel_print.C"

// main function
// argument determines whether we fit or open fit.txt and plot from there
// called with > root -l 'c_test.C($tofit)'
void c_test_print(int tofit)
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
  if(tofit == 1)  fitclass.doFit();
  else fitclass.readRes();

  if(tofit!=2) fitclass.doPlot(plot_div);
  else {
    int npartot = fitclass.nparam[0]+fitclass.nparam[1]+fitclass.nparam[2];
    double chipar[npartot];
    for(int i = 0; i < npartot; i++)
      chipar[i] = fitclass.fitpar[i];
    fitclass.myFunction(chipar);
  }
  
  cout << "Fit option given was " << tofit << ", therefore";
  if(tofit) cout << " fit was performed" << endl;
  else if(tofit == 0) cout << " previous fit results were read for the plotting" << endl;
  else cout << "just called minimum functon" << endl;

  /*vector <vector <double> > datasigma = fitclass.datasigma;
  int l = datasigma.size();
  int lin = datasigma[0].size();
  cout << l << " " << lin << endl;
  cout << datasigma[0][0] << endl;
  for(int i = 0; i < l; i++)
    {
      for(int j = 0; j < lin; j++)
	cout << datasigma[i][j] << " ";
      cout << endl << endl;
      }*/
}
