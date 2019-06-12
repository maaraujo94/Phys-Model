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

  ///////////////////////////////
  //part: plotting data + pulls
  ///////////////////////////////
  
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
  counter = csplotpull(datasigma, counter, ndata[xplot], len, mqq[xplot], param, Lpos[xplot], state[xplot], lumibr, ptmmin, minimum, ndf, chiprob, mkrStyle, legtitles, savename[xplot]);
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
}
