#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TFitter.h"
#include "TMath.h"
#include <fstream>
#include <iostream>
using namespace std;

double sqsmNorm = 9.;
double cosaNorm = 0.;

// function to parse a string into components separated by "deli"
vector< string > parseString( string line, string deli) {
  vector< string > out;
  string aux = line;
  size_t pos = 0;
  
  while( (pos = aux.find(deli)) != string::npos)
    {
      out.push_back( aux.substr( 0, pos ) );
      aux.erase( 0, pos + deli.length() );	    
    }
  out.push_back( aux );
  
  return out;
}

//auxiliary classes for the three types of parameters (normalization, shape and nuisance)
class ShapePar{
public:
  string name;
  double val;
  double unc;
  int fix;
  void display() {
    cout << name << " " << val << " " << unc << " " << fix << endl;
  }
};

class NormPar{
public:
  string name;
  double mult;
  int fix;
  void display() {
    cout << name << " " << mult << " " << fix << endl;
  }
};

class NuisPar{
public:
  string name;
  double normunc;
  int fix;
  void display() {
    cout << name << " " << normunc << " " << fix << endl;
  }
};

//the parametrization of the pdf function
double pdf(double x, double b, double c, double d, double e)
{
  if (x>1)
    return 0;
  return pow(x, -b) * pow(1-x, c) * (1 + d*sqrt(x) + e*x);
}

//the function to be integrated
//par-> [0]: sqsfm, [1]: pTM, [2]: y, [3]: y4, [4]: A, [5]: beta, [6]: tau, [7]: rho, [8]: delta, [9]: b, [10]: c, [11]: d, [12]: e
double integf(double *par)
{
  double eps = 1e-5;
  
  double p = sqrt( par[1] * par[1] * pow( cosh(0.5*(par[2]-par[3])), 2) + pow( sinh(0.5*(par[2]-par[3])), 2) );
  double sqsm = p + sqrt( 1 + p*p );
  double s = sqsm*sqsm;
  double en = (s+1) / (2*sqsm);
  double x1 = sqsm/par[0] * exp( 0.5*(par[2]+par[3]) );
  double x2 = sqsm/par[0] * exp( -0.5*(par[2]+par[3]) );
  double cosa = (sqrt( 1 + par[1]*par[1] ) * sinh(0.5*(par[2]-par[3]))) / p;

  double fx1 = pdf(x1,par[9],par[10],par[11], par[12]);
  double fx2 = pdf(x2,par[9],par[10],par[11], par[12]);

  //the dsigma/dt^star function already multiplied by sqrt(s^star)*p^star
  double fcos = (pow(en + p*cosa, -par[7]) + pow(en - p*cosa, -par[7])) * pow(1+eps-cosa*cosa, -par[8]) * pow(en, par[7]) / 2;
  double fs;
  if(abs(par[4]) < 1e-5)
    fs = pow(p, par[6])*pow(s,-par[5]-par[6]/2);//pow(1+s/par[4], -par[5]-par[6]/2);
  else
    fs = pow(p, par[6])*pow(1+s/par[4], -par[5]-par[6]/2);
  double dsdt = fcos * fs;

  return fx1 * fx2 * dsdt;
}

//the function that performs the integration, corresponds to dsigma/dxidy
//par-> [0]: sqsfm, [1]: pT/M, [2]: y, [3]: L, [4]: M, [5]: A, [6]: beta, [7]: tau, [8]: rho, [9]: delta, [10]: b, [11]: c, [12]: d, [13]: e
double sig(vector <double> par)
{
  double sum = 0;
  
  double ylim = acosh( 1/par[1] * ( par[0] - sqrt(1+par[1]*par[1]) * cosh(par[2])));
  double fin = 26;
  double dy = 2*ylim/(fin-1.);
  double y4 = -ylim;
  double pars[13]={par[0], par[1], par[2], y4, par[5], par[6], par[7], par[8], par[9], par[10], par[11], par[12], par[13]};

  for(int i=0; i < fin; i++)
    {
      pars[3] = y4;
      sum += ( par[1] / pow(par[0], 2) ) * integf(pars) * dy * ( par[3] / pow(par[4],5));
      y4 += dy;
    }
 
  //get the normalization value
  double sNorm = sqsmNorm*sqsmNorm;
  
  double pNorm = (sNorm-1)/(2*sqsmNorm);
  double eNorm = (sNorm+1)/(2*sqsmNorm);

  double fcosNorm = (pow(eNorm + pNorm*cosaNorm, -par[8]) + pow(eNorm - pNorm*cosaNorm, -par[8])) * pow(1+1e-5-cosaNorm*cosaNorm, -par[9]) * pow(eNorm, par[8]) / 2;
  double fsNorm = pow(pNorm, par[7])*pow(sNorm,-par[6]-par[7]/2);
  double dsdtNorm = fcosNorm*fsNorm;
  
  return sum/dsdtNorm;
}

//sigma function for plotting
//takes params individually rather than through pointer 
double sigplot(double sqsfm, double pTM, double y, double L, double M, double A, double beta, double tau, double rho, double delta, double b, double c, double d, double e)
{
  double sum = 0;
  
  double ylim = acosh( 1/pTM * ( sqsfm - sqrt(1+pTM*pTM) * cosh(y)));
  double fin = 26;
  double dy = 2*ylim/(fin-1.);
  double y4 = -ylim;
  double pars[13]={sqsfm, pTM, y, y4, A, beta, tau, rho, delta, b, c, d, e};
  
  for(int i=0; i < fin; i++)
    {
      pars[3] = y4;
      sum += ( pTM / pow(sqsfm, 2) ) * integf(pars) * dy * ( L / pow(M,5));
      y4 += dy;
    }

  //get the normalization value
  double sNorm = sqsmNorm*sqsmNorm;
  
  double pNorm = (sNorm-1)/(2*sqsmNorm);
  double eNorm = (sNorm+1)/(2*sqsmNorm);

  double fcosNorm = (pow(eNorm + pNorm*cosaNorm, -rho) + pow(eNorm - pNorm*cosaNorm, -rho)) * pow(1+1e-5-cosaNorm*cosaNorm, -delta) * pow(eNorm, rho) / 2;
  double fsNorm = pow(pNorm, tau)*pow(sNorm,-beta-tau/2);
  double dsdtNorm = fcosNorm*fsNorm;
  
  return sum/dsdtNorm;
 
  return sum;
}

//calculates the average pT/M of a bin by weighing it with the model cs
//uint and lint are the upper (sigma*pT/M) and lower (sigma) integrals
double avgptm(double lbound, double ubound, double lybound, double uybound, vector <double> par)
{
  double uint = 0, lint = 0, dpt = ubound-lbound, dy = uybound-lybound;
  int npt = 8, ny = 4;
  
  for(int xpt = 0; xpt < npt; xpt++) {
    par[1] = lbound + xpt*dpt/(npt-1.)+1e-10;
    for(int xy = 0; xy < ny; xy++) {
      par[2] = lybound + xy*dy/(ny-1.);
      uint += sig(par)*par[1];
      lint += sig(par);
    }
  }
  return uint / lint;
}

int getStyle(string det) {
  if(det == "CMS") return 20;
  if(det == "ATLAS") return 25;
  if(det == "LHCb") return 22;
  else return 1;
}

void aRange(string det, string state, double *apos) {
  apos[0] = -1;
  apos[3] = 9.99e3;
  if(det == "CMS" && state == "jpsi") {
    apos[1] = 1.01e-5;
    apos[2] = 39.9;
  }
  else {
    apos[1] = 1.01e-4;
    apos[2] = 10.9;
  }
}

// get relative position on an axis (pi, pf)
double getPos(double pi, double pf, double mult, bool isLog) {
  if(isLog) return pow(10, log10(pi)+mult*(log10(pf)-log10(pi)));
  else return pi + mult*(pf-pi);
}
