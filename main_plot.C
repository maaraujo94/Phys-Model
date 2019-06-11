#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TAxis.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TStyle.h"
#include <fstream>
using namespace std;

#include "aux_func.C"
#include "aux_input.C"

vector <vector <double> > datasigma(10);
//psi', chi_c2, chi_c1, jpsi, ups(3S), ups(2S), ups(1S)
float mass[7] = {3.686, 3.556, 3.511, 3.097, 10.355, 10.023, 9.460};

void plot()
{
  gROOT->SetBatch(1);
  
  ifstream file;
  ofstream tex;
  const int nstates = 18;
  int ns[nstates] = {0}, counter = 0, ndf;
  const int nparam = 32;
  double param[nparam], eparam[nparam], ptmmin = 0, minimum, chiprob;

  //////////////////////////////////
  //part: getting data from files
  /////////////////////////////////

  input(datasigma, ns);

  /////////////////////////////
  //part: reading fit results
  /////////////////////////////
  
  file.open("fit.txt");
  for(int i = 0; i < nparam; i++)
    file >> param[i];
  for(int i = 0; i < nparam; i++)
    file >> eparam[i];
  file >> ptmmin;
  file >> minimum;
  file >> ndf;
  chiprob = TMath::Prob(minimum, ndf);

  char names[nparam][100] = {"L_{\\psi(2S)}","L_{\\chi_{c2}}","L_{\\chi_{c1}}","L_{J/\\psi}","L_{\\Upsilon(3S)}","L_{\\Upsilon(2S)}","L_{\\chi_{b2}}","L_{\\chi_{b1}}","L_{\\Upsilon(1S)}","A","\\beta","\\tau","\\rho","\\delta","b","c","d","e","BR_{ppdm}","BR_{ppjdp}","BR_{c2jpsi}","BR_{c1jpsi}","BR_{jpsidm}","BR_{b2ups1}","BR_{b1ups1}","BR_{ups3dm}","BR_{ups2dm}","BR_{ups1dm}","L_{CMS}","L_{ATLAS}","L_{ATLAS}(\\Upsilon)","L_{CMS}(13)"};

  //writing the fit parameters to a pretty tex file
  tex.open("plots/fitp.tex");

  tex << "\\begin{table}" << endl;
  tex << "\\centering" << endl;
  tex << "\\begin{tabular}{c|c|c}" << endl;
  tex << "Parameter & Value & Uncertainty \\\\" << endl;
  tex << "\\hline" << endl;
  for(int i = 0; i < nparam; i++)
    {
      tex << "$" << names[i] << "$ & " << param[i] << " & ";
      if(eparam[i] == 0) tex << "fixed \\\\" << endl;
      else tex << eparam[i] << " \\\\" << endl;
    }
  tex << "min $p_T/M$ & " << ptmmin << " & fixed \\\\" << endl;
  tex << "\\end{tabular}" << endl;
  tex << "\\caption{Fit parameters ($\\chi^2$ / ndf = " << minimum << " / " << ndf << " = " << minimum/ndf << ")}" << endl;
  tex << "\\end{table}" << endl;
  tex.close();

  /////////////////////////
  //part: plotting data
  /////////////////////////
  
  //auxiliary functions storing plotting variables for all 14 plots
  //psi', chic2, chic1, jpsi, ups3, ups2, ups1, psi'13, jpsi13, ups313, ups213, ups113, chicr(state->Lpos2), chibr(state->Lpos2)
  const int nplots = 14;
  int ndata[nplots] = {2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1};
  int Lpos[nplots] = {0, 1, 2, 3, 4, 5, 8, 0, 3, 4, 5, 8, 1, 6};
  int state[nplots] = {0, 1, 2, 3, 4, 5, 6, 0, 3, 4, 5, 6, 2, 7};
  double mqq[nplots];
  string savename[nplots] = {"psiprime", "chic2", "chic1", "jpsi", "ups3S", "ups2S", "ups1S", "psiprime_13", "jpsi_13", "ups3S_13", "ups2S_13", "ups1S_13", "chicr", "chibr"};
  for(int i = 0; i < nplots; i++) {
    mqq[i] = mass[state[i]];
    //savename[i] = "plots/"+savename[i]+"_cs.pdf";
    savename[i] = "plots/"+savename[i];
  }
  mqq[nplots-2] = mass[3];
  mqq[nplots-1] = mass[6];

  //xplot - plot number, dst - dataset number
  int xplot = 0, dst = 0;
  double lumibr[2];
  string legtitles[2];
  int mkrStyle[2], len[2];
  
  //psiprime 7 TeV
  //customizable options
  lumibr[0] = param[28]*param[18];
  lumibr[1] = param[29]*param[19]*param[22];
  legtitles[0] = "CMS #psi(2S)";
  legtitles[1] = "ATLAS #psi(2S)";
  mkrStyle[0] = 20;
  mkrStyle[1] = 25;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;
  
  //chic2 7 TeV
  //customizable options
  lumibr[0] = param[29]*param[20]*param[22];
  legtitles[0] = "ATLAS #chi_{c2}";
  mkrStyle[0] = 25;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //chic1 7 TeV
  //customizable options
  lumibr[0] = param[21]*param[22];
  legtitles[0] = "ATLAS #chi_{c1}";
  mkrStyle[0] = 25;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;
  
  //jpsi 7 TeV
  //customizable options
  lumibr[0] = param[28]*param[22];
  legtitles[0] = "CMS J/#psi";
  mkrStyle[0] = 20;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //ups(3S) 7 TeV
  //customizable options
  lumibr[0] = param[28]*param[25];
  lumibr[1] = param[25]*param[30];
  legtitles[0] = "CMS #Upsilon(3S)";
  legtitles[1] = "ATLAS #Upsilon(3S)";
  mkrStyle[0] = 20;
  mkrStyle[1] = 25;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //ups(2S) 7 TeV
  //customizable options
  lumibr[0] = param[28]*param[26];
  lumibr[1] = param[26]*param[30];
  legtitles[0] = "CMS #Upsilon(2S)";
  legtitles[1] = "ATLAS #Upsilon(2S)";
  mkrStyle[0] = 20;
  mkrStyle[1] = 25;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //ups(1S) 7 TeV
  //customizable options
  lumibr[0] = param[28]*param[27];
  lumibr[1] = param[27]*param[30];
  legtitles[0] = "CMS #Upsilon(1S)";
  legtitles[1] = "ATLAS #Upsilon(1S)";
  mkrStyle[0] = 20;
  mkrStyle[1] = 25;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //psiprime 13 TeV
  //customizable options
  lumibr[0] = param[31]*param[18];
  legtitles[0] = "CMS #psi(2S)";
  mkrStyle[0] = 20;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //jpsi 13 TeV
  //customizable options
  lumibr[0] = param[31]*param[22];
  legtitles[0] = "CMS J/#psi";
  mkrStyle[0] = 20;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //ups(3S) 13 TeV
  //customizable options
  lumibr[0] = param[31]*param[25];
  legtitles[0] = "CMS #Upsilon(3S)";
  mkrStyle[0] = 20;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //ups(2S) 13 TeV
  //customizable options
  lumibr[0] = param[31]*param[26];
  legtitles[0] = "CMS #Upsilon(2S)";
  mkrStyle[0] = 20;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //ups(1S) 13 TeV
  //customizable options
  lumibr[0] = param[31]*param[27];
  legtitles[0] = "CMS #Upsilon(1S)";
  mkrStyle[0] = 20;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //chic2 / chic1 ratio 7 TeV
  //customizable options
  lumibr[0] = param[20]/param[21];
  legtitles[0] = "CMS #chi_{c2} / #chi_{c1}";
  mkrStyle[0] = 20;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = rplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //chib2 / chib1 ratio 7 TeV
  //customizable options
  lumibr[0] = param[23]/param[24];
  legtitles[0] = "CMS #chi_{b2} / #chi_{b1}";
  mkrStyle[0] = 20;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = rplot(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  ////////////////////////////////////////////////////////////////////
  //part: plotting the "pulls" - maintaining the uncertainty, though
  ////////////////////////////////////////////////////////////////////

  counter = 0;
  xplot = 0;
  dst = 0;
  
  //psiprime 7 TeV
  //customizable options
  lumibr[0] = param[28]*param[18];
  lumibr[1] = param[29]*param[19]*param[22];
  legtitles[0] = "CMS #psi(2S)";
  legtitles[1] = "ATLAS #psi(2S)";
  mkrStyle[0] = 20;
  mkrStyle[1] = 25;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = cspull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;
  
  
  /*double cs, y, pt, dpt, dy, npt;
  int nmin = 3, ny=4, nsteps;
  double dn = (8.-3.)/39., signc = sigplot(datasigma[9][0]/(mass[0]), 4, 0, param[9], param[10], param[11], param[12], param[13], 1., mass[0], param[14], param[15], param[16], param[17]);
  
  //7 TeV psiprime data
  // CMS psiprime
  counter=0;
  float pintpsipccs[2][npsipccs];
  for(int i=0; i<npsipccs; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[0]), pt, y, param[9], param[10], param[11], param[12], param[13], param[0], mass[0], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintpsipccs[0][i] = (datapsipccs[1][i] - cs ) / cs;
      pintpsipccs[1][i] = datapsipccs[2][i] / cs;
    }
 
  counter+=npsipccs;
  //ATLAS psiprime
  float pintpsipacs[2][npsipacs];
  for(int i=0; i<npsipacs; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[0]), pt, y, param[9], param[10], param[11], param[12], param[13], param[0], mass[0], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintpsipacs[0][i] = (datapsipacs[1][i] - cs ) / cs;
      pintpsipacs[1][i] = datapsipacs[2][i] / cs;
    }
  
  TCanvas *ppsip = new TCanvas("pull psip", "pull psip", 700, 700);
  TH1F *fppsip = ppsip->DrawFrame(xi, -2, xf, 2);
  fppsip->SetXTitle("p_{T}/M");
  fppsip->SetYTitle("pulls");
  fppsip->GetYaxis()->SetTitleOffset(1);
  ppsip->Modified();
  ppsip->SetTitle("");

  TF1 *zero = new TF1("zero", "0", xi, xf);
  zero->SetLineColor(kBlue);
  zero->SetLineStyle(7);
  zero->Draw("lsame");
  
  //plot CMS data
  TGraphAsymmErrors *pcpull = new TGraphAsymmErrors(npsipccs, datapsipccs[0], pintpsipccs[0], datapsipccs[3], datapsipccs[4], pintpsipccs[1], pintpsipccs[1]);
  pcpull->SetMarkerStyle(20);
  pcpull->SetLineColor(kBlack);
  pcpull->SetMarkerColor(kBlack);
  pcpull->SetMarkerSize(.75);
  pcpull->Draw("P");

  TGraphAsymmErrors *papull = new TGraphAsymmErrors(npsipacs, datapsipacs[0], pintpsipacs[0], datapsipacs[3], datapsipacs[4], pintpsipacs[1], pintpsipacs[1]);
  papull->SetMarkerStyle(25);
  papull->SetLineColor(kBlack);
  papull->SetMarkerColor(kBlack);
  papull->SetMarkerSize(.75);
  papull->Draw("P");

  //text on the plot
  TLatex lppsip;
  lppsip.SetTextSize(0.03);
  lppsip.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lppsip.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lppsip.DrawLatex(0, -1.7, "pp 7 TeV");
  
  //draw legend
  TLegend *legppsip = new TLegend(0.6, 0.75, 0.9, 0.9);
  legppsip->SetTextSize(0.03);
  legppsip->AddEntry(pcpull, "CMS #psi(2S)", "p");
  legppsip->AddEntry(papull, "ATLAS #psi(2S)", "p");
  legppsip->Draw();
  //save psiprime cross section plot
  ppsip->SaveAs("plots/psiprime_pull.pdf");
  
  //7 TeV chic2 data
  counter+=npsipacs;
  float pintchic2cs[2][nchic2cs];
  for(int i=0; i<nchic2cs; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[1]), pt, y, param[9], param[10], param[11], param[12], param[13], param[1], mass[1], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);

      pintchic2cs[0][i] = (datachic2cs[1][i] - cs ) / cs;
      pintchic2cs[1][i] = datachic2cs[2][i] / cs;
    }

  TCanvas *pchic2 = new TCanvas("pull chic2", "pull chic2", 700, 700);
  TH1F *fpchic2 = pchic2->DrawFrame(xi, -2, xf, 2);
  fpchic2->SetXTitle("p_{T}/M");
  fpchic2->SetYTitle("pulls");
  fpchic2->GetYaxis()->SetTitleOffset(1);
  pchic2->Modified();
  pchic2->SetTitle("");

  zero -> Draw("lsame");

  //plot CMS data
  TGraphAsymmErrors *c2pull = new TGraphAsymmErrors(nchic2cs, datachic2cs[0], pintchic2cs[0], datachic2cs[3], datachic2cs[4], pintchic2cs[1], pintchic2cs[1]);
  c2pull->SetMarkerStyle(20);
  c2pull->SetLineColor(kBlack);
  c2pull->SetMarkerColor(kBlack);
  c2pull->SetMarkerSize(.75);
  c2pull->Draw("P");
  
  //text on the plot
  TLatex lpchic2;
  lpchic2.SetTextSize(0.03);
  lpchic2.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpchic2.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpchic2.DrawLatex(0, -1.7, "pp 7 TeV");
  
  //draw legend
  TLegend *legpchic2 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legpchic2->SetTextSize(0.03);
  legpchic2->AddEntry(c2pull, "ATLAS #chi_{c2}", "p");
  legpchic2->Draw();
  //save psiprime cross section plot
  pchic2->SaveAs("plots/chic2_pull.pdf");

  //7 TeV chic1 data
  counter+=nchic2cs;
  float pintchic1cs[2][nchic1cs];
  for(int i=0; i<nchic1cs; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[2]), pt, y, param[9], param[10], param[11], param[12], param[13], param[2], mass[2], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);

      pintchic1cs[0][i] = (datachic1cs[1][i] - cs ) / cs;
      pintchic1cs[1][i] = datachic1cs[2][i] / cs;
    }

  TCanvas *pchic1 = new TCanvas("pull chic1", "pull chic1", 700, 700);
  TH1F *fpchic1 = pchic1->DrawFrame(xi, -2, xf, 2);
  fpchic1->SetXTitle("p_{T}/M");
  fpchic1->SetYTitle("pulls");
  fpchic1->GetYaxis()->SetTitleOffset(1);
  pchic1->Modified();
  pchic1->SetTitle("");

  zero -> Draw("lsame");

  //plot CMS data
  TGraphAsymmErrors *c1pull = new TGraphAsymmErrors(nchic1cs, datachic1cs[0], pintchic1cs[0], datachic1cs[3], datachic1cs[4], pintchic1cs[1], pintchic1cs[1]);
  c1pull->SetMarkerStyle(20);
  c1pull->SetLineColor(kBlack);
  c1pull->SetMarkerColor(kBlack);
  c1pull->SetMarkerSize(.75);
  c1pull->Draw("P");
  
  //text on the plot
  TLatex lpchic1;
  lpchic1.SetTextSize(0.03);
  lpchic1.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpchic1.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpchic1.DrawLatex(0, -1.7, "pp 7 TeV");
  
  //draw legend
  TLegend *legpchic1 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legpchic1->SetTextSize(0.03);
  legpchic1->AddEntry(c1pull, "ATLAS #chi_{c1}", "p");
  legpchic1->Draw();
  //save psiprime cross section plot
  pchic1->SaveAs("plots/chic1_pull.pdf");
  
  //7 TeV jpsi data
  counter+=nchic1cs;
  float pintjpsics[2][njpsics];
  for(int i=0; i<njpsics; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[3]), pt, y, param[9], param[10], param[11], param[12], param[13], param[3], mass[3], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);

      pintjpsics[0][i] = (datajpsics[1][i] - cs ) / cs;
      pintjpsics[1][i] = datajpsics[2][i] / cs;
    }

  TCanvas *pjpsi = new TCanvas("pull jpsi", "pull jpsi", 700, 700);
  TH1F *fpjpsi = pjpsi->DrawFrame(xi, -2, xf, 2);
  fpjpsi->SetXTitle("p_{T}/M");
  fpjpsi->SetYTitle("pulls");
  fpjpsi->GetYaxis()->SetTitleOffset(1);
  pjpsi->Modified();
  pjpsi->SetTitle("");

  zero -> Draw("lsame");

  //plot CMS data
  TGraphAsymmErrors *jppull = new TGraphAsymmErrors(njpsics, datajpsics[0], pintjpsics[0], datajpsics[3], datajpsics[4], pintjpsics[1], pintjpsics[1]);
  jppull->SetMarkerStyle(20);
  jppull->SetLineColor(kBlack);
  jppull->SetMarkerColor(kBlack);
  jppull->SetMarkerSize(.75);
  jppull->Draw("P");
  
  //text on the plot
  TLatex lpjpsi;
  lpjpsi.SetTextSize(0.03);
  lpjpsi.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpjpsi.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpjpsi.DrawLatex(0, -1.7, "pp 7 TeV");
  
  //draw legend
  TLegend *legpjpsi = new TLegend(0.6, 0.8, 0.9, 0.9);
  legpjpsi->SetTextSize(0.03);
  legpjpsi->AddEntry(c1pull, "CMS J/#psi", "p");
  legpjpsi->Draw();
  //save psiprime cross section plot
  pjpsi->SaveAs("plots/jpsi_pull.pdf");

  //7 TeV ups(3S) data
  // CMS ups(3S)
  counter+=njpsics+nchicr+nchibr;
  float pintups3ccs[2][nups3ccs];
  for(int i=0; i<nups3ccs; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[4]), pt, y, param[9], param[10], param[11], param[12], param[13], param[4], mass[4], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintups3ccs[0][i] = (dataups3ccs[1][i] - cs ) / cs;
      pintups3ccs[1][i] = dataups3ccs[2][i] / cs;
    }
 
  counter+=nups3ccs;
  //ATLAS ups(3S)
  float pintups3acs[2][nups3acs];
  for(int i=0; i<nups3acs; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[4]), pt, y, param[9], param[10], param[11], param[12], param[13], param[4], mass[4], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintups3acs[0][i] = (dataups3acs[1][i] - cs ) / cs;
      pintups3acs[1][i] = dataups3acs[2][i] / cs;
    }    
  
  TCanvas *pups3 = new TCanvas("pull ups3", "pull ups3", 700, 700);
  TH1F *fpups3 = pups3->DrawFrame(xi, -2, xf, 2);
  fpups3->SetXTitle("p_{T}/M");
  fpups3->SetYTitle("pulls");
  fpups3->GetYaxis()->SetTitleOffset(1);
  pups3->Modified();
  pups3->SetTitle("");

  zero->Draw("lsame");
  
  //plot CMS data
  TGraphAsymmErrors *u3cpull = new TGraphAsymmErrors(nups3ccs, dataups3ccs[0], pintups3ccs[0], dataups3ccs[3], dataups3ccs[4], pintups3ccs[1], pintups3ccs[1]);
  u3cpull->SetMarkerStyle(20);
  u3cpull->SetLineColor(kBlack);
  u3cpull->SetMarkerColor(kBlack);
  u3cpull->SetMarkerSize(.75);
  u3cpull->Draw("P");

  for(int i=0; i<nups3ccs; i++)
    if(dataups3ccs[0][i] < ptmmin)
      u3cpull->RemovePoint(0);
  u3cpull->Draw("P");
  
  TGraphAsymmErrors *u3apull = new TGraphAsymmErrors(nups3acs, dataups3acs[0], pintups3acs[0], dataups3acs[3], dataups3acs[4], pintups3acs[1], pintups3acs[1]);
  u3apull->SetMarkerStyle(25);
  u3apull->SetLineColor(kBlack);
  u3apull->SetMarkerColor(kBlack);
  u3apull->SetMarkerSize(.75);
  u3apull->Draw("P");

  for(int i=0; i<nups3acs; i++)
    if(dataups3acs[0][i] < ptmmin)
      u3apull->RemovePoint(0);
  u3apull->Draw("P");
  
  //text on the plot
  TLatex lpups3;
  lpups3.SetTextSize(0.03);
  lpups3.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpups3.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpups3.DrawLatex(0, -1.7, "pp 7 TeV");
  
  //draw legend
  TLegend *legpups3 = new TLegend(0.6, 0.75, 0.9, 0.9);
  legpups3->SetTextSize(0.03);
  legpups3->AddEntry(u3cpull, "CMS #Upsilon(3S)", "p");
  legpups3->AddEntry(u3apull, "ATLAS #Upsilon(3S)", "p");
  legpups3->Draw();
  //save ups3 cross section plot
  pups3->SaveAs("plots/ups3S_pull.pdf");

  //7 TeV ups(2S) data
  // CMS ups(2S)
  counter+=nups3acs;
  float pintups2ccs[2][nups2ccs];
  for(int i=0; i<nups2ccs; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[5]), pt, y, param[9], param[10], param[11], param[12], param[13], param[5], mass[5], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintups2ccs[0][i] = (dataups2ccs[1][i] - cs ) / cs;
      pintups2ccs[1][i] = dataups2ccs[2][i] / cs;
    }
 
  counter+=nups2ccs;
  //ATLAS ups(2S)
  float pintups2acs[2][nups2acs];
  for(int i=0; i<nups2acs; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[5]), pt, y, param[9], param[10], param[11], param[12], param[13], param[5], mass[5], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintups2acs[0][i] = (dataups2acs[1][i] - cs ) / cs;
      pintups2acs[1][i] = dataups2acs[2][i] / cs;
    }    
  
  TCanvas *pups2 = new TCanvas("pull ups2", "pull ups2", 700, 700);
  TH1F *fpups2 = pups2->DrawFrame(xi, -2, xf, 2);
  fpups2->SetXTitle("p_{T}/M");
  fpups2->SetYTitle("pulls");
  fpups2->GetYaxis()->SetTitleOffset(1);
  pups2->Modified();
  pups2->SetTitle("");

  zero->Draw("lsame");
  
  //plot CMS data
  TGraphAsymmErrors *u2cpull = new TGraphAsymmErrors(nups2ccs, dataups2ccs[0], pintups2ccs[0], dataups2ccs[3], dataups2ccs[4], pintups2ccs[1], pintups2ccs[1]);
  u2cpull->SetMarkerStyle(20);
  u2cpull->SetLineColor(kBlack);
  u2cpull->SetMarkerColor(kBlack);
  u2cpull->SetMarkerSize(.75);
  u2cpull->Draw("P");

  for(int i=0; i<nups2ccs; i++)
    if(dataups2ccs[0][i] < ptmmin)
      u2cpull->RemovePoint(0);
  u2cpull->Draw("P");
  
  TGraphAsymmErrors *u2apull = new TGraphAsymmErrors(nups2acs, dataups2acs[0], pintups2acs[0], dataups2acs[3], dataups2acs[4], pintups2acs[1], pintups2acs[1]);
  u2apull->SetMarkerStyle(25);
  u2apull->SetLineColor(kBlack);
  u2apull->SetMarkerColor(kBlack);
  u2apull->SetMarkerSize(.75);
  u2apull->Draw("P");

  for(int i=0; i<nups2acs; i++)
    if(dataups2acs[0][i] < ptmmin)
      u2apull->RemovePoint(0);
  u2apull->Draw("P");
  
  //text on the plot
  TLatex lpups2;
  lpups2.SetTextSize(0.03);
  lpups2.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpups2.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpups2.DrawLatex(0, -1.7, "pp 7 TeV");
  
  //draw legend
  TLegend *legpups2 = new TLegend(0.6, 0.75, 0.9, 0.9);
  legpups2->SetTextSize(0.03);
  legpups2->AddEntry(u2cpull, "CMS #Upsilon(2S)", "p");
  legpups2->AddEntry(u2apull, "ATLAS #Upsilon(2S)", "p");
  legpups2->Draw();
  //save ups2 cross section plot
  pups2->SaveAs("plots/ups2S_pull.pdf");

  //7 TeV ups(1S) data
  // CMS ups(1S)
  counter+=nups2acs;
  float pintups1ccs[2][nups1ccs];
  for(int i=0; i<nups1ccs; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[6]), pt, y, param[9], param[10], param[11], param[12], param[13], param[8], mass[6], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintups1ccs[0][i] = (dataups1ccs[1][i] - cs ) / cs;
      pintups1ccs[1][i] = dataups1ccs[2][i] / cs;
    }
 
  counter+=nups1ccs;
  //ATLAS ups(1S)
  float pintups1acs[2][nups1acs];
  for(int i=0; i<nups1acs; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[6]), pt, y, param[9], param[10], param[11], param[12], param[13], param[8], mass[6], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintups1acs[0][i] = (dataups1acs[1][i] - cs ) / cs;
      pintups1acs[1][i] = dataups1acs[2][i] / cs;
    }    
  
  TCanvas *pups1 = new TCanvas("pull ups1", "pull ups1", 700, 700);
  TH1F *fpups1 = pups1->DrawFrame(xi, -2, xf, 2);
  fpups1->SetXTitle("p_{T}/M");
  fpups1->SetYTitle("pulls");
  fpups1->GetYaxis()->SetTitleOffset(1);
  pups1->Modified();
  pups1->SetTitle("");

  zero->Draw("lsame");
  
  //plot CMS data
  TGraphAsymmErrors *u1cpull = new TGraphAsymmErrors(nups1ccs, dataups1ccs[0], pintups1ccs[0], dataups1ccs[3], dataups1ccs[4], pintups1ccs[1], pintups1ccs[1]);
  u1cpull->SetMarkerStyle(20);
  u1cpull->SetLineColor(kBlack);
  u1cpull->SetMarkerColor(kBlack);
  u1cpull->SetMarkerSize(.75);
  u1cpull->Draw("P");

  for(int i=0; i<nups1ccs; i++)
    if(dataups1ccs[0][i] < ptmmin)
      u1cpull->RemovePoint(0);
  u1cpull->Draw("P");
  
  TGraphAsymmErrors *u1apull = new TGraphAsymmErrors(nups1acs, dataups1acs[0], pintups1acs[0], dataups1acs[3], dataups1acs[4], pintups1acs[1], pintups1acs[1]);
  u1apull->SetMarkerStyle(25);
  u1apull->SetLineColor(kBlack);
  u1apull->SetMarkerColor(kBlack);
  u1apull->SetMarkerSize(.75);
  u1apull->Draw("P");

  for(int i=0; i<nups1acs; i++)
    if(dataups1acs[0][i] < ptmmin)
      u1apull->RemovePoint(0);
  u1apull->Draw("P");
  
  //text on the plot
  TLatex lpups1;
  lpups1.SetTextSize(0.03);
  lpups1.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpups1.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpups1.DrawLatex(0, -1.7, "pp 7 TeV");
  
  //draw legend
  TLegend *legpups1 = new TLegend(0.6, 0.75, 0.9, 0.9);
  legpups1->SetTextSize(0.03);
  legpups1->AddEntry(u1cpull, "CMS #Upsilon(1S)", "p");
  legpups1->AddEntry(u1apull, "ATLAS #Upsilon(1S)", "p");
  legpups1->Draw();
  //save ups2 cross section plot
  pups1->SaveAs("plots/ups1S_pull.pdf");
  
  //13 TeV psiprime data
  counter += nups1acs;
  float pintpsipcs13[2][npsipccs13];
  for(int i=0; i<npsipccs13; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[0]), pt, y, param[9], param[10], param[11], param[12], param[13], param[0], mass[0], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintpsipcs13[0][i] = (datapsipccs13[1][i] - cs ) / cs;
      pintpsipcs13[1][i] = datapsipccs13[2][i] / cs;
    }

  TCanvas *ppsip13 = new TCanvas("pull psip13", "pull psip13", 700, 700);
  TH1F *fppsip13 = ppsip13->DrawFrame(xi, -2, xf, 2);
  fppsip13->SetXTitle("p_{T}/M");
  fppsip13->SetYTitle("pulls");
  fppsip13->GetYaxis()->SetTitleOffset(1);
  ppsip13->Modified();
  ppsip13->SetTitle("");

  zero -> Draw("lsame");

  //plot CMS data
  TGraphAsymmErrors *pc13pull = new TGraphAsymmErrors(npsipccs13, datapsipccs13[0], pintpsipcs13[0], datapsipccs13[3], datapsipccs13[4], pintpsipcs13[1], pintpsipcs13[1]);
  pc13pull->SetMarkerStyle(20);
  pc13pull->SetLineColor(kBlack);
  pc13pull->SetMarkerColor(kBlack);
  pc13pull->SetMarkerSize(.75);
  pc13pull->Draw("P");

  //text on the plot
  TLatex lppsip13;
  lppsip13.SetTextSize(0.03);
  lppsip13.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lppsip13.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lppsip13.DrawLatex(0, -1.7, "pp 13 TeV");
  
  //draw legend
  TLegend *legppsip13 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legppsip13->SetTextSize(0.03);
  legppsip13->AddEntry(pc13pull, "CMS #psi(2S)", "p");
  legppsip13->Draw();
  //save psiprime cross section plot
  ppsip13->SaveAs("plots/psiprime_pull_13.pdf");
  
  //13 TeV jpsi data
  counter+=npsipccs13;
  float pintjpsics13[2][njpsics13];
  for(int i=0; i<njpsics13; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[3]), pt, y, param[9], param[10], param[11], param[12], param[13], param[3], mass[3], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);

      pintjpsics13[0][i] = (datajpsics13[1][i] - cs ) / cs;
      pintjpsics13[1][i] = datajpsics13[2][i] / cs;
    }
  
  TCanvas *pjpsi13 = new TCanvas("pull jpsi13", "pull jpsi13", 700, 700);
  TH1F *fpjpsi13 = pjpsi13->DrawFrame(xi, -2, xf, 2);
  fpjpsi13->SetXTitle("p_{T}/M");
  fpjpsi13->SetYTitle("pulls");
  fpjpsi13->GetYaxis()->SetTitleOffset(1);
  pjpsi13->Modified();
  pjpsi13->SetTitle("");

  zero -> Draw("lsame");

  //plot CMS data
  TGraphAsymmErrors *jp13pull = new TGraphAsymmErrors(njpsics13, datajpsics13[0], pintjpsics13[0], datajpsics13[3], datajpsics13[4], pintjpsics13[1], pintjpsics13[1]);
  jp13pull->SetMarkerStyle(20);
  jp13pull->SetLineColor(kBlack);
  jp13pull->SetMarkerColor(kBlack);
  jp13pull->SetMarkerSize(.75);
  jp13pull->Draw("P");

  //text on the plot
  TLatex lpjpsi13;
  lpjpsi13.SetTextSize(0.03);
  lpjpsi13.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpjpsi13.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpjpsi13.DrawLatex(0, -1.7, "pp 13 TeV");
  
  //draw legend
  TLegend *legpjpsi13 = new TLegend(0.6, 0.8, 0.9, 0.9);
  legpjpsi13->SetTextSize(0.03);
  legpjpsi13->AddEntry(jp13pull, "CMS J/#psi", "p");
  legpjpsi13->Draw();
  //save psiprime cross section plot
  pjpsi13->SaveAs("plots/jpsi_pull_13.pdf");

  //13 TeV ups(3S) data
  // CMS ups(3S)
  counter+=njpsics13;
  float pintups3ccs13[2][nups3ccs13];
  for(int i=0; i<nups3ccs13; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[4]), pt, y, param[9], param[10], param[11], param[12], param[13], param[4], mass[4], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintups3ccs13[0][i] = (dataups3ccs13[1][i] - cs ) / cs;
      pintups3ccs13[1][i] = dataups3ccs13[2][i] / cs;
    }
 
  TCanvas *pups313 = new TCanvas("pull ups3 13", "pull ups3 13", 700, 700);
  TH1F *fpups313 = pups313->DrawFrame(xi, -2, xf, 2);
  fpups313->SetXTitle("p_{T}/M");
  fpups313->SetYTitle("pulls");
  fpups313->GetYaxis()->SetTitleOffset(1);
  pups313->Modified();
  pups313->SetTitle("");

  zero->Draw("lsame");
  
  //plot CMS data
  TGraphAsymmErrors *u3cpull13 = new TGraphAsymmErrors(nups3ccs13, dataups3ccs13[0], pintups3ccs13[0], dataups3ccs13[3], dataups3ccs13[4], pintups3ccs13[1], pintups3ccs13[1]);
  u3cpull13->SetMarkerStyle(20);
  u3cpull13->SetLineColor(kBlack);
  u3cpull13->SetMarkerColor(kBlack);
  u3cpull13->SetMarkerSize(.75);
  u3cpull13->Draw("P");

  for(int i=0; i<nups3ccs13; i++)
    if(dataups3ccs13[0][i] < ptmmin)
      u3cpull13->RemovePoint(0);
  u3cpull13->Draw("P");
  
  //text on the plot
  TLatex lpups313;
  lpups313.SetTextSize(0.03);
  lpups313.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpups313.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpups313.DrawLatex(0, -1.7, "pp 13 TeV");
  
  //draw legend
  TLegend *legpups313 = new TLegend(0.6, 0.75, 0.9, 0.9);
  legpups313->SetTextSize(0.03);
  legpups313->AddEntry(u3cpull13, "CMS #Upsilon(3S)", "p");
  legpups313->Draw();
  //save ups3 cross section plot
  pups313->SaveAs("plots/ups3S_pull_13.pdf");

  //13 TeV ups(2S) data
  // CMS ups(2S)
  counter+=nups3ccs13;
  float pintups2ccs13[2][nups2ccs13];
  for(int i=0; i<nups2ccs13; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[5]), pt, y, param[9], param[10], param[11], param[12], param[13], param[5], mass[5], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintups2ccs13[0][i] = (dataups2ccs13[1][i] - cs ) / cs;
      pintups2ccs13[1][i] = dataups2ccs13[2][i] / cs;
    }
 
  TCanvas *pups213 = new TCanvas("pull ups2 13", "pull ups2 13", 700, 700);
  TH1F *fpups213 = pups213->DrawFrame(xi, -2, xf, 2);
  fpups213->SetXTitle("p_{T}/M");
  fpups213->SetYTitle("pulls");
  fpups213->GetYaxis()->SetTitleOffset(1);
  pups213->Modified();
  pups213->SetTitle("");

  zero->Draw("lsame");
  
  //plot CMS data
  TGraphAsymmErrors *u2cpull13 = new TGraphAsymmErrors(nups2ccs13, dataups2ccs13[0], pintups2ccs13[0], dataups2ccs13[3], dataups2ccs13[4], pintups2ccs13[1], pintups2ccs13[1]);
  u2cpull13->SetMarkerStyle(20);
  u2cpull13->SetLineColor(kBlack);
  u2cpull13->SetMarkerColor(kBlack);
  u2cpull13->SetMarkerSize(.75);
  u2cpull13->Draw("P");

  for(int i=0; i<nups2ccs13; i++)
    if(dataups2ccs13[0][i] < ptmmin)
      u2cpull13->RemovePoint(0);
  u2cpull13->Draw("P");
  
  //text on the plot
  TLatex lpups213;
  lpups213.SetTextSize(0.03);
  lpups213.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpups213.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpups213.DrawLatex(0, -1.7, "pp 13 TeV");
  
  //draw legend
  TLegend *legpups213 = new TLegend(0.6, 0.75, 0.9, 0.9);
  legpups213->SetTextSize(0.03);
  legpups213->AddEntry(u2cpull13, "CMS #Upsilon(2S)", "p");
  legpups213->Draw();
  //save ups2 cross section plot
  pups213->SaveAs("plots/ups2S_pull_13.pdf");

  //7 TeV ups(1S) data
  // CMS ups(1S)
  counter+=nups2ccs13;
  float pintups1ccs13[2][nups1ccs13];
  for(int i=0; i<nups1ccs13; i++)
    {
      cs=0.;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs+=sigplot(datasigma[9][i+counter]/(mass[6]), pt, y, param[9], param[10], param[11], param[12], param[13], param[8], mass[6], param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(signc*nsteps);
      
      pintups1ccs13[0][i] = (dataups1ccs13[1][i] - cs ) / cs;
      pintups1ccs13[1][i] = dataups1ccs13[2][i] / cs;
    }
 
  TCanvas *pups113 = new TCanvas("pull ups1 13", "pull ups1 13", 700, 700);
  TH1F *fpups113 = pups113->DrawFrame(xi, -2, xf, 2);
  fpups113->SetXTitle("p_{T}/M");
  fpups113->SetYTitle("pulls");
  fpups113->GetYaxis()->SetTitleOffset(1);
  pups113->Modified();
  pups113->SetTitle("");

  zero->Draw("lsame");
  
  //plot CMS data
  TGraphAsymmErrors *u1cpull13 = new TGraphAsymmErrors(nups1ccs13, dataups1ccs13[0], pintups1ccs13[0], dataups1ccs13[3], dataups1ccs13[4], pintups1ccs13[1], pintups1ccs13[1]);
  u1cpull13->SetMarkerStyle(20);
  u1cpull13->SetLineColor(kBlack);
  u1cpull13->SetMarkerColor(kBlack);
  u1cpull13->SetMarkerSize(.75);
  u1cpull13->Draw("P");

  for(int i=0; i<nups1ccs13; i++)
    if(dataups1ccs13[0][i] < ptmmin)
      u1cpull13->RemovePoint(0);
  u1cpull13->Draw("P");
  
  //text on the plot
  TLatex lpups113;
  lpups113.SetTextSize(0.03);
  lpups113.DrawLatex(8, 1.7, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpups113.DrawLatex(8, 1.3, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpups113.DrawLatex(0, -1.7, "pp 13 TeV");
  
  //draw legend
  TLegend *legpups113 = new TLegend(0.6, 0.75, 0.9, 0.9);
  legpups113->SetTextSize(0.03);
  legpups113->AddEntry(u1cpull13, "CMS #Upsilon(1S)", "p");
  legpups113->Draw();
  //save ups2 cross section plot
  pups113->SaveAs("plots/ups1S_pull_13.pdf");*/
}
