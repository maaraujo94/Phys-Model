/* calls on the members of a fitting class (in fitModel) to
 - read options
 - get input
 - perform fit (optional)
 - plot and save results
 "print" means it's the version used to store the cosalpha distribution */

#import "fitModel_print.C"

// main function
// reads fit result, calls minimum() and stores cosa values scanned
// called with > root -l c_test_print.C
void c_test_print()
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

  // read fit results
  fitclass.readRes();

  // call minimum once to store cosalpha
  int npartot = fitclass.nparam[0]+fitclass.nparam[1]+fitclass.nparam[2];
  double chipar[npartot];
  for(int i = 0; i < npartot; i++)
    chipar[i] = fitclass.fitpar[i];
  fitclass.myFunction(chipar);
  
  cout << "c_test_print: called minimum function" << endl;
}
