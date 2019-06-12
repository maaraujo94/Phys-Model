////////////////////////////////////////////////////////////////////
// function that initializes fit values and performs the fitting  //
////////////////////////////////////////////////////////////////////

#include "TFitter.h"
#include "TMath.h"
#include "TROOT.h"
#include <fstream>
#include <iostream>
using namespace std;

#include "aux_func.C"
#include "aux_input.C"

vector <vector <double> > datasigma(10);

const int nstates = 18;
int ns[nstates] = {0};

//psi', chi_c2, chi_c1, jpsi, ups(3S), ups(2S), ups(1S)
float mass[7] = {3.686, 3.556, 3.511, 3.097, 10.355, 10.023, 9.460};

//////////////////////////////////////
//part: the log-likelihood function
/////////////////////////////////////

//defining the log-likelihood
double myFunction(double Lpsip, double Lchic2, double Lchic1, double Ljpsi, double Lups3, double Lups2, double Lchib2, double Lchib1, double Lups1, double A, double beta, double tau, double rho, double delta, double b, double c, double d, double e, double brpsipdm, double brpsipdp, double brc2jpsi, double brc1jpsi, double brjpsidm, double brb2ups1, double brb1ups1, double brups3dm, double brups2dm, double brups1dm, double lcms, double latlas, double latlasups, double lcms13, double ptmmin)
{
  Lchib2 = Lchic2;
  Lchib1 = Lchic1;
  double xinorm = 4, ssnorm = 7000, ynorm = 0, Mnorm = mass[0], loglike = 0.;
  int pos = 0, stateid;

  //defining the arguments of sig(). Must redefine pars[0] when fitting data at different sqrt(s) (13 TeV, etc)
  double pars[14] = {ssnorm/Mnorm, xinorm, ynorm, A, beta, tau, rho, delta, 1., Mnorm, b, c, d, e};
  //normalize all states at xinorm, ynorm
  double signc = sig(pars);

  //store all state-specific parameters here
  double Mstate[nstates] = {mass[0], mass[0], mass[1], mass[2], mass[3], mass[4], mass[4], mass[5], mass[5], mass[6], mass[6], mass[0], mass[3], mass[4], mass[5], mass[6], mass[3], mass[6]};
  double Lstate[nstates-2] = {Lpsip, Lpsip, Lchic2, Lchic1, Ljpsi, Lups3, Lups3, Lups2, Lups2, Lups1, Lups1, Lpsip, Ljpsi, Lups3, Lups2, Lups1};
  double lumistate[nstates-2] = {lcms, latlas, latlas, latlas, lcms, lcms, latlasups, lcms, latlasups, lcms, latlasups, lcms13, lcms13, lcms13, lcms13, lcms13};
  double BRstate[nstates] = {brpsipdm, brpsipdp*brjpsidm, brc2jpsi*brjpsidm, brc1jpsi*brjpsidm, brjpsidm, brups3dm, brups3dm, brups2dm, brups2dm, brups1dm, brups1dm, brpsipdm, brjpsidm, brups3dm, brups2dm, brups1dm, brc2jpsi/brc1jpsi, brc2jpsi/brc1jpsi};

  //cycle for cross section measurements
  for(int stateid = 0; stateid < nstates-2; stateid++)
    if(ns[stateid] != 0) {
      loglike += statecont(datasigma, pars, signc, ptmmin, Mstate[stateid], Lstate[stateid], lumistate[stateid], BRstate[stateid], pos, ns[stateid]);
      pos += ns[stateid];
    }
  //cycle for cross section ratio measurements
  for(int stateid = nstates-2; stateid < nstates; stateid++)
    if(ns[stateid] != 0) {
      loglike += ratiocont(datasigma, pars, signc, ptmmin, Mstate[stateid], Lchic2, Lchic1, BRstate[stateid], pos, ns[stateid]);
      pos += ns[stateid];
    } 
  
  //constraints for the BR nuisance parameters
  loglike -= 0.5*((brpsipdm-1)/(0.9/7.9))*((brpsipdm-1)/(0.9/7.9));
  loglike -= 0.5*((brpsipdp-1)/(0.3/34.46))*((brpsipdp-1)/(0.3/34.46));
  loglike -= 0.5*((brc2jpsi-1)/(0.7/19.2))*((brc2jpsi-1)/(0.7/19.2));
  loglike -= 0.5*((brc1jpsi-1)/(1.2/33.9))*((brc1jpsi-1)/(1.2/33.9));
  loglike -= 0.5*((brjpsidm-1)/(0.033/5.961))*((brjpsidm-1)/(0.033/5.961));
  loglike -= 0.5*((brb2ups1-1)/(1.2/19.1))*((brb2ups1-1)/(1.2/19.1));
  loglike -= 0.5*((brb1ups1-1)/(2.2/33.9))*((brb1ups1-1)/(2.2/33.9));
  loglike -= 0.5*((brups3dm-1)/(0.21/2.18))*((brups3dm-1)/(0.21/2.18));
  loglike -= 0.5*((brups2dm-1)/(0.17/1.93))*((brups2dm-1)/(0.17/1.93));
  loglike -= 0.5*((brups1dm-1)/(0.05/2.48))*((brups1dm-1)/(0.05/2.48));
  
  //constraints for the luminosity nuisance parameters
  loglike -= 0.5*((lcms-1)/2.2e-2)*((lcms-1)/2.2e-2);
  loglike -= 0.5*((latlas-1)/1.8e-2)*((latlas-1)/1.8e-2);
  loglike -= 0.5*((latlasups-1)/3.9e-2)*((latlasups-1)/3.9e-2);
  loglike -= 0.5*((lcms13-1)/2.2e-2)*((lcms13-1)/2.2e-2);
  
  //I am returning the chisquare because I can't figure out how to change the fit method (must check this)
  return -2*loglike;
}

//Minuit-specific "wrapping"
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg)
{
  result = myFunction(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], par[10], par[11], par[12], par[13], par[14], par[15], par[16], par[17], par[18], par[19], par[20], par[21], par[22], par[23], par[24], par[25], par[26], par[27], par[28], par[29], par[30], par[31], par[32]);
}

//main function: collects data, does the fit, plots the results
void fit()
{
  //define variables
  
  ofstream outf;
  int ndf, aux_int = 0;
  
  const int nparam = 32;
  double param[nparam], ptmmin = 2., chisquare, L_est[7];
  float PTMNORM = 4.;
  
  //////////////////////////////////
  //part: getting data from files
  /////////////////////////////////

  input(datasigma, ns);
  
  ////////////////////////////
  //part: fitting the data
  ///////////////////////////

  //define function to minimize, set initial parameters
  TFitter* fit = new TFitter(nparam+1);
  fit->SetFCN(minuitFunction);

  //calculate estimates for L_QQ
  int len[7] = {ns[0]+ns[1], ns[2], ns[3], ns[4], ns[5]+ns[6], ns[7]+ns[8], ns[9]+ns[10]};
  for(int i = 0; i < 7; i++) {
    L_est[i] = lest(datasigma, PTMNORM, aux_int, len[i]);
    aux_int += len[i];
  }

  fit->SetParameter(0, "L_psip", L_est[0], L_est[0]/10, 0, 0);
  fit->SetParameter(1, "L_chic2", L_est[1], L_est[1]/10, 0, 0);
  fit->SetParameter(2, "L_chic1", L_est[2], L_est[2]/10, 0, 0);
  fit->SetParameter(3, "L_jpsi", L_est[3], L_est[3]/10, 0, 0);
  fit->SetParameter(4, "L_ups3", L_est[4]*200, L_est[4]*20, 0, 0);
  fit->SetParameter(5, "L_ups2", L_est[5]*200, L_est[5]*20, 0, 0);
  fit->SetParameter(6, "L_chib2", L_est[5]*200, L_est[5]*20, 0, 0);
  fit->FixParameter(6);
  fit->SetParameter(7, "L_chib1", L_est[5]*200, L_est[5]*20, 0, 0);
  fit->FixParameter(7);
  fit->SetParameter(8, "L_ups1", L_est[6]*200, L_est[6]*20, 0, 0);

  fit->SetParameter(9, "A", 0., 0.1, 0, 0);
  fit->FixParameter(9);
  fit->SetParameter(10, "beta", 1.5, 0.1, 0, 0);
  fit->SetParameter(11, "tau", 1., 0.1, 0, 0);
  fit->SetParameter(12, "rho", 2.4, 0.1, 0, 0);
  fit->FixParameter(12);
  fit->SetParameter(13, "delta", 0.1, 0.1, 0, 0);
  fit->FixParameter(13);

  fit->SetParameter(14, "b", 1.8, 0.1, 0, 0);
  fit->SetParameter(15, "c", 6, 0.1, 0, 0);
  fit->SetParameter(16, "d", 0., 0.1, 0, 0);
  fit->FixParameter(16);
  fit->SetParameter(17, "e", 0., 0.1, 0, 0);
  fit->FixParameter(17);

  fit->SetParameter(18, "br_pp_dm", 1., 0.1, 0, 0);
  fit->SetParameter(19, "br_pp_jdp", 1., 0.01, 0, 0);
  fit->SetParameter(20, "br_c2_jpsi", 1., 0.01, 0, 0);
  fit->SetParameter(21, "br_c1_jpsi", 1., 0.01, 0, 0);
  fit->SetParameter(22, "br_jpsi_dm", 1., 0.01, 0, 0);
  fit->SetParameter(23, "br_b2_ups1", 1, 0.01, 0, 0);
  fit->SetParameter(24, "br_b1_ups1", 1, 0.01, 0, 0);
  fit->SetParameter(25, "br_ups3_dm", 1, 0.1, 0, 0);
  fit->SetParameter(26, "br_ups2_dm", 1, 0.1, 0, 0);
  fit->SetParameter(27, "br_ups1_dm", 1, 0.01, 0, 0);
  fit->SetParameter(28, "L_CMS", 1, 0.1, 0, 0);
  fit->SetParameter(29, "L_ATLAS", 1, 0.1, 0, 0);
  fit->SetParameter(30, "L_ATL_ups", 1, 0.1, 0, 0);
  fit->SetParameter(31, "L_CMS_13", 1, 0.1, 0, 0);
  aux_int = 14;
  
  fit->SetParameter(32, "ptm_min", ptmmin, 1, 0, 0);
  fit->FixParameter(32);

  fit->ExecuteCommand("MIGRAD",0,0);

  //////////////////////////////
  // part: storing fit results
  //////////////////////////////
  
  //register the fit parameters to an array
  for(int i = 0; i < nparam; i++)
    param[i] = fit->GetParameter(i);
  param[6] = param[1];
  param[7] = param[2];
  
  chisquare = myFunction(param[0], param[1], param[2], param[3], param[4], param[5], param[6], param[7], param[8], param[9], param[10], param[11], param[12], param[13], param[14], param[15], param[16], param[17], param[18], param[19], param[20], param[21], param[22], param[23], param[24], param[25], param[26], param[27], param[28], param[29], param[30], param[31], ptmmin);

  ndf = 0;
  //calculating the number of data points used in the fit
  for(int i = 0; i < (int)datasigma[0].size(); i++)
    {
      if(datasigma[0][i] >= ptmmin)
	ndf += 1;
      else if(datasigma[0][i] < ptmmin && (datasigma[3][i] == 16 || datasigma[3][i] == 17))
      ndf += 1;
    }
  cout << ndf << " data points" << endl;
  //I am adding one degree of freedom for each of the nuisance parameters
  ndf += aux_int;
  //the final ndf is obtained by subtracting the number of free parameters
  ndf -= fit->GetNumberFreeParameters();
  cout << "ndf " << ndf << endl;
     
  cout << "chi^2 norm : " << chisquare/ndf  << endl;
  cout << "chi^2 prob : " << TMath::Prob(chisquare,ndf)  << endl;

  //////////////////////////////////////
  //part: saving results to an outfile
  //////////////////////////////////////
  
  outf.open("fit.txt");
  for(int i = 0; i < nparam; i++)
    outf << param[i] << " ";
  for(int i = 0; i < nparam; i++)
    outf << fit->GetParError(i) << " ";
  outf << ptmmin << " ";
  outf << chisquare << " " << ndf << endl;
  outf.close();

}
