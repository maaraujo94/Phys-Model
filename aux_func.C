/* auxiliary macro with all the functions that are independent of the fit class
 - aux classes for param types, string parser, style functions
 - all functions that make up the cross-section model used */

#include "LHAPDF/LHAPDF.h"
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

using namespace LHAPDF;
using namespace std;

// normalization point (sqrt(s), xi, y)*
double sqsNorm = 7000.;
double xiNorm = 5.;
double yNorm = 1.;

PDF *pdf_ct = mkPDF("CT14nnlo", 0);

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

//the parametrization of the pdf function with LHAPDF NNLO_CTEQ14 
double pdf_new(double x, double q2, double s)
{
  double x_min = q2 / s;
  if ( x > 1 || x < x_min)
    return 0;
  return pdf_ct->xfxQ2(21, x, q2);
}

// the alpha function (currently unused)
double alpha(double q2)
{
  return pdf_ct->alphasQ2(q2);
}

//the function to be integrated
//par-> [0]: sqrt(s)/M, [1]: xi, [2]: y, [3]: cosalpha, [4]: M, [5]: beta, [6]: rho, [7]: delta
double integf(double *par)
{
  double eps = 1e-5;

  double pstar = par[1] / sqrt(1.+eps-par[3]*par[3]);
  double sqstar = pstar + sqrt(1.+pstar*pstar);
  double sstar = sqstar*sqstar;
  double Estar = (sstar+1)/(2*sqrt(sstar));
  
  double pLstar = pstar * par[3];
  double y_hat = 0.5*log((Estar + pLstar)/(Estar - pLstar));
  double y_gg = par[2] - y_hat;
  double x1 = sqstar/par[0]*exp(y_gg);
  double x2 = sqstar/par[0]*exp(-y_gg);

  double s_hat = sstar*par[4]*par[4];
  double s = pow(par[0]*par[4], 2);
  double fx1 = pdf_new(x1, s_hat, s);
  double fx2 = pdf_new(x2, s_hat, s);

  //the dsigma/dt^star function multiplied by mass-normalized jacobian
  double gfunc = ( pow( Estar+pstar*par[3], -par[6] ) + pow( Estar-pstar*par[3], -par[6] ) ) * ( pow( Estar, par[6] ) / 2.) * pow( 1 + eps - par[3]*par[3], -par[7]);
  double hfunc = pow(sstar, -par[5]);
 
  double dsdt = gfunc * hfunc;
  double jac = pstar/Estar * (sstar-1.);
  
  return fx1 * fx2 * jac * dsdt;
}

//the function that performs the integration, corresponds to dsigma/dxidy
//par-> [0]: sqrt(s)/M, [1]: xi, [2]: y, [3]: M, [4]: beta, [5]: rho, [6]: delta
double sig(double *par)
{
  double sum = 0;

  double cosa_lim = 1;
  double d_cosa = 0.1;
  int n_cosa = (2.*cosa_lim)/d_cosa;

  double cosa = -cosa_lim;
  double pars[8]={par[0], par[1], par[2], cosa, par[3], par[4], par[5], par[6]};
  
  for(int i=0; i <= n_cosa; i++)
    {
      pars[3] = cosa;
      sum += 1. / ( pow(par[3],5) * par[1]) * d_cosa * integf(pars);
      cosa += d_cosa;
    }
  
  return sum;
}

// the cross section function with normalization applied
//par-> [0]: sqrt(s)/M, [1]: xi, [2]: y, [3]: L, [4]: M, [5]: beta, [6]: rho, [7]: delta
double sigN(vector <double> par)
{
  double pars[7] = {par[0], par[1], par[2], par[4], par[5], par[6], par[7]};
  double sigFull = sig(pars);

  pars[0] = sqsNorm/par[4];
  pars[1] = xiNorm;
  pars[2] = yNorm;
  double sigNorm = sig(pars);

  return par[3] * sigFull/sigNorm;
}

//sigma function for plotting
//takes params individually rather than through pointer
double sigplot(double sqs_M, double xi, double y, double L, double M, double beta, double rho, double delta)
{
  double pars[7] = {sqs_M, xi, y, M, beta, rho, delta};
  double sigFull = sig(pars);

  pars[0] = sqsNorm/M;
  pars[1] = xiNorm;
  pars[2] = yNorm;
  double sigNorm = sig(pars);

  return L * sigFull/sigNorm;
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
      uint += sigN(par)*par[1];
      lint += sigN(par);
    }
  }
  return uint / lint;
}

// define marker style from data "experiment" label
int getStyle(string det) {
  if(det == "CMS") return 20;
  if(det == "ATLAS") return 25;
  if(det == "LHCb") return 22;
  else return 1;
}

// define plot range depending on data label
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
