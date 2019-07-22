///////////////////////////////////////////////////////////////////////
// macro to read files outputed by main_fit_cycle.C and plot results //
///////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include <fstream>
using namespace std;

#include "aux_func.C"
#include "aux_input.C"

void plot_pars()
{
  ifstream file;
  //these n values must match those of main_fit_cycle.C
  const int n_tau = 4, n_s = 3, n_L = 3;

  double param[n_L][n_s][n_s][n_s][n_tau][nparam], eparam[n_L][n_s][n_s][n_s][n_tau][nparam], minimum[n_L][n_s][n_s][n_s][n_tau], aux, chiprob[n_L][n_s][n_s][n_s][n_tau];
  int ndf[n_L][n_s][n_s][n_s][n_tau], par;
  float taupts[n_tau], tauerr[n_tau];
  
  //parameters to plot. 0/1 depending on which J/psi or Ups(1S) are being fitted
  int toplot[4] = {1, 3, 5, 7};
  const int ntotal = n_L * n_s * n_s * n_s;
  
  char names[nparam][100] = {"L(Jpsi)", "L(Upsilon(1S))", "A", "beta", "tau", "rho", "delta", "b", "c", "d", "e", "BR(jpsidm)", "BR(ups1dm)", "L(CMS)", "L(ATLAS)", "L(LHCb,psi)", "L(LHCb,Upsilon)"};

  //cycle to read results for all fits in cycle
  for(int iL = 0; iL < n_L; iL++)
    for(int ibeta= 0; ibeta < n_s; ibeta++)
      for(int irho = 0; irho < n_s; irho++)
	for(int ib = 0; ib < n_s; ib++) 
	  for(int itau = 0; itau < n_tau; itau++) {
	    string name = Form("fit_cycle/fit_%d%d%d%d%d.txt", iL, ibeta, irho, ib, itau); 
	    file.open(name);
	    for(int i = 0; i < nparam; i++)
	      file >> param[iL][ibeta][irho][ib][itau][i];
	    for(int i = 0; i < nparam; i++)
	      file >> eparam[iL][ibeta][irho][ib][itau][i];
	    file >> aux;
	    file >> minimum[iL][ibeta][irho][ib][itau];
	    file >> ndf[iL][ibeta][irho][ib][itau];
	    chiprob[iL][ibeta][irho][ib][itau] = TMath::Prob(minimum[iL][ibeta][irho][ib][itau], ndf[iL][ibeta][irho][ib][itau]);
	    file.close();
	  }
	  
  for(int i = 0; i < n_tau; i++) {
    taupts[i] = i;
    tauerr[i] = 0;
  }
  
  //part: plotting

  TCanvas *c = new TCanvas("param", "param", 700, 700);
  c->SetLogy();
  
  int xi = -1, xf = n_tau;
  //float yi[5] = {10, 0, 0, 0, 7};
  float yi[4] = {1e-1, 1e-1, 1e-1, 1e-2};
  //float yf[5] = {30, 5, 4, 4, 20};
  float yf[4] = {1000, 10, 10, 10};
  
  float parpts[n_L][n_s][n_s][n_s][n_tau], parerr[n_L][n_s][n_s][n_s][n_tau];

  //cycle plots the four parameters we want
  for(int i = 0; i< 4; i++){
    par = toplot[i];
    TH1F *fc = c->DrawFrame(xi, yi[i], xf, yf[i]);
    fc->SetXTitle("tau");
    string name = Form("%s vs tau", names[par]);
    fc->SetYTitle(names[par]);
    fc->GetYaxis()->SetTitleOffset(1.2);
    fc->SetTitle(name.c_str());
    c->Modified();
    
    for(int iL = 0; iL < n_L; iL++)
      for(int ibeta= 0; ibeta < n_s; ibeta++)
	for(int irho = 0; irho < n_s; irho++)
	  for(int ib = 0; ib < n_s; ib++) 
	    for(int itau = 0; itau < n_tau; itau++) {
	      parpts[iL][ibeta][irho][ib][itau] = param[iL][ibeta][irho][ib][itau][par];
	      parerr[iL][ibeta][irho][ib][itau] = eparam[iL][ibeta][irho][ib][itau][par];
	    }
    //cycle to define TGraphAsymmErrors
    TGraphErrors** gc = new TGraphErrors*[ntotal];
    
    int itot = 0;
    for(int iL = 0; iL < n_L; iL++)
      for(int ibeta= 0; ibeta < n_s; ibeta++)
	for(int irho = 0; irho < n_s; irho++)
	  for(int ib = 0; ib < n_s; ib++) {
	    gc[itot] = new TGraphErrors(n_tau, taupts, parpts[iL][ibeta][irho][ib], tauerr, parerr[iL][ibeta][irho][ib]);
	    gc[itot]->SetMarkerStyle(20);
	    gc[itot]->SetMarkerSize(.75);
	    gc[itot]->SetMarkerColor((itot%9)+1);
	    gc[itot]->SetLineColor((itot%9)+1);
	    gc[itot]->Draw("PL");
	    itot++;
	  }
    
    //save plot
    string savename = Form("plots_cycle/%svsTau.pdf", names[par]);
    const char* save = savename.c_str();
    c->SaveAs(save);
    c->Clear();
    
  }

  //repeat what's done above but to plot chi square
  //c->SetLogy(0);
  TH1F *fc = c->DrawFrame(xi, 100, xf, 1e10);
  fc->SetXTitle("tau");
  fc->SetYTitle("chisquare");
  fc->GetYaxis()->SetTitleOffset(1.2);
  fc->SetTitle("chisquare vs tau");
  c->Modified();
  
  for(int iL = 0; iL < n_L; iL++)
    for(int ibeta= 0; ibeta < n_s; ibeta++)
      for(int irho = 0; irho < n_s; irho++)
	for(int ib = 0; ib < n_s; ib++) 
	  for(int itau = 0; itau < n_tau; itau++) {
	    parpts[iL][ibeta][irho][ib][itau] = minimum[iL][ibeta][irho][ib][itau];
	    parerr[iL][ibeta][irho][ib][itau] =0;
	  }

  //cycle to define TGraphAsymmErrors
  TGraphErrors** gc = new TGraphErrors*[ntotal];
    
  int itot = 0;
  for(int iL = 0; iL < n_L; iL++)
    for(int ibeta= 0; ibeta < n_s; ibeta++)
      for(int irho = 0; irho < n_s; irho++)
	for(int ib = 0; ib < n_s; ib++) {
	  gc[itot] = new TGraphErrors(n_tau, taupts, parpts[iL][ibeta][irho][ib], tauerr, parerr[iL][ibeta][irho][ib]);
	  gc[itot]->SetMarkerStyle(20);
	  gc[itot]->SetMarkerSize(.75);
	  gc[itot]->SetMarkerColor((itot%9)+1);
	  gc[itot]->SetLineColor((itot%9)+1);
	  gc[itot]->Draw("PL");
	  itot++;
	}
  
  //save plot
  string savename = "plots_cycle/chisquarevsTau.pdf";
  const char* save = savename.c_str();
  c->SaveAs(save);
  //c->Clear();
  
  c->Destructor();
}
