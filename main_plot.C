/////////////////////////////////////////////////////////
// function that reads the fit results and plots them  //
/////////////////////////////////////////////////////////

#include "TMath.h"
#include "TROOT.h"
#include <fstream>
using namespace std;

#include "aux_func.C"
#include "aux_input.C"

vector <vector <double> > datasigma(10);

void plot()
{
  gROOT->SetBatch(1);
  
  ifstream file;
  ofstream tex;
  int ns[nstates] = {0}, counter = 0, ndf;
  double param[nparam], eparam[nparam], ptmmin = 0, minimum, chiprob;

  //////////////////////////////////
  //part: getting data from files
  //////////////////////////////////

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

  char names[nparam][100] = {"L_{J/\\psi}", "L_{\\Upsilon(1S)}", "A", "\\beta", "\\tau", "\\rho", "\\delta", "b", "c", "d", "e", "BR_{jpsidm}", "BR_{ups1dm}", "L_{CMS}", "L_{ATLAS}", "L_{LHCb}(\\psi)", "L_{LHCb}(\\Upsilon)"};

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

  ///////////////////////////////
  //part: plotting data + pulls
  ///////////////////////////////

  //auxiliary functions storing plotting variables for all plots
  //jpsi,  ups1

  const int nplots = 4;
  int ndata[nplots] = {1, 5, 3, 5};
  int Lpos[nplots] = {0, 0, 1, 1};
  int state[nplots] = {0, 0, 4, 4};
  double mqq[nplots];
  string savename[nplots] = {"jpsi_CA", "jpsi_LHCb", "ups1S_CA", "ups1S_LHCb"};
  for(int i = 0; i < nplots; i++) {
    mqq[i] = mass[state[i]];
    savename[i] = "plots/"+savename[i];
  }
  
  //xplot - plot number, dst - dataset number
  int xplot = 0, dst = 0;
  double lumibr[5];
  string legtitles[5];
  int mkrStyle[5], mkrCol[5], len[5], fplot[5] = {1, 1, 1, 1, 1};
  
  //CMS+ATLAS jpsi 7 TeV
  //customizable options
  lumibr[0] = param[11]*param[13];
  legtitles[0] = "CMS J/#psi";
  mkrStyle[0] = 20;
  mkrCol[0] = kBlack;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, mkrCol, fplot, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //LHCb jpsi 7 TeV
  for(int i = 0; i < 5; i++) {
    lumibr[i] = param[15];
    legtitles[i] = Form("LHCb J/#psi y%d", i+1);
    mkrStyle[i] = 22;
    mkrCol[i] = i+1;
  }
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, mkrCol, fplot, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //ups(1S) 7 TeV
  //customizable options
  lumibr[0] = param[12]*param[13];
  legtitles[0] = "CMS #Upsilon(1S)";
  mkrStyle[0] = 20;
  mkrCol[0] = kBlack;
  for(int i = 0; i < 2; i++) {
    lumibr[i+1] = param[12]*param[14];
    legtitles[i+1] = Form("ATLAS #Upsilon(1S) y%d", i+1);
    mkrStyle[i+1] = 25;
    mkrCol[i+1] = i+1;
  }
  //fplot[2] = 1;
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, mkrCol, fplot, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;
  
  for(int i = 0; i < 5; i++) {
    lumibr[i] = param[12]*param[16];
    legtitles[i] = Form("LHCb #Upsilon(1S) y%d", i+1);
    mkrStyle[i] = 22;
    mkrCol[i] = i+1;
  }
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, mkrCol, fplot, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  
  /*const int nplots = 2;
  int ndata[nplots] = {6, 8};
  int Lpos[nplots] = {0, 1};
  int state[nplots] = {0, 4};
  double mqq[nplots];
  string savename[nplots] = {"jpsi", "ups1S"};
  for(int i = 0; i < nplots; i++) {
    mqq[i] = mass[state[i]];
    savename[i] = "plots/"+savename[i];
  }
  
  //xplot - plot number, dst - dataset number
  int xplot = 0, dst = 0;
  double lumibr[8];
  string legtitles[8];
  int mkrStyle[8], mkrCol[8], len[8];
  
  //jpsi 7 TeV
  //customizable options
  lumibr[0] = param[11]*param[13];
  legtitles[0] = "CMS J/#psi";
  mkrStyle[0] = 20;
  mkrCol[0] = kBlack;
  for(int i = 0; i < 5; i++) {
    lumibr[i+1] = param[15];
    legtitles[i+1] = Form("LHCb J/#psi y%d", i+1);
    mkrStyle[i+1] = 22;
    mkrCol[i+1] = i+1;
  }
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, mkrCol, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;

  //ups(1S) 7 TeV
  //customizable options
  lumibr[0] = param[12]*param[13];
  legtitles[0] = "CMS #Upsilon(1S)";
  mkrStyle[0] = 20;
  mkrCol[0] = kBlack;
  for(int i = 0; i < 2; i++) {
    lumibr[i+1] = param[12]*param[14];
    legtitles[i+1] = Form("ATLAS #Upsilon(1S) y%d", i+1);
    mkrStyle[i+1] = 25;
    mkrCol[i+1] = i+1;
  }
  for(int i = 0; i < 5; i++) {
    lumibr[i+3] = param[12]*param[16];
    legtitles[i+3] = Form("LHCb #Upsilon(1S) y%d", i+1);
    mkrStyle[i+3] = 22;
    mkrCol[i+3] = i+1;
  }
  //this structure is constant for all plots
  for(int j = 0; j < ndata[xplot]; j++)
    len[j] = ns[dst+j];
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, mkrCol, legtitles, savename[xplot]);
  dst+=ndata[xplot];
  xplot++;*/
}
