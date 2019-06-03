#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
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
float mass[7]={3.686, 3.556, 3.511, 3.097, 10.355, 10.023, 9.460};

////////////////////////////////////////////////////////////////////////////////////////////
//part: definition of auxiliary functions to obtain sig(pT, y, sqrt(s)/M, fitting params)
///////////////////////////////////////////////////////////////////////////////////////////

void plot()
{
  ifstream file;
  ofstream tex;
  std::string data;
  int ns[18]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  char nums[12][100];
  double param[36];
  double eparam[36];
  double ptmmin=0;
  int stateid=0, counter=0;
  int ndf;
  
  int pos_psip=0, pos_chic2=0, pos_chic1=0, pos_jpsi=0, pos_ups3=0, pos_ups2=0, pos_ups1=0;
  double pt_psip=100, pt_chic2=100, pt_chic1=100, pt_jpsi=100, pt_ups3=100, pt_ups2=100, pt_ups1=100;
  double minimum, chinorm, chiprob, cs_psip, cs_chic2, cs_chic1, cs_jpsi, cs_ups3, cs_ups2, cs_chib2, cs_chib1, cs_ups1, lambda_chic2, lambda_chic1, lambda_jpsi, lambda_ups3, lambda_ups2, lambda_chib2, lambda_chib1, lambda_ups1, ratio, aux;

  double dataout[75];

  //////////////////////////////////
  //part: getting data from files
  /////////////////////////////////

  input(datasigma, ns);
  
  //part: reading fit results

  file.open("fit.txt");
  for(int i=0; i<75; i++)
    file >> dataout[i];
  for(int i=0; i<36; i++)
    param[i] = dataout[i];
  for(int i=0; i<36; i++)
    eparam[i] = dataout[i+36];
  //for(int i=18; i<31; i++)
  //  param[i] = 1;
  ptmmin = dataout[72];
  //ptmmin = 1e-5;
  minimum = dataout[73];
  ndf = dataout[74];
  chiprob = TMath::Prob(minimum, ndf);

  char names[36][100] = {"L_{\\psi(2S)}","L_{\\chi_{c2}}","L_{\\chi_{c1}}","L_{J/\\psi}","L_{\\Upsilon(3S)}","L_{\\Upsilon(2S)}","L_{\\chi_{b2}}","L_{\\chi_{b1}}","L_{\\Upsilon(1S)}","A","\\beta","\\tau","\\rho","\\delta","b","c","d","e","BR_{ppdm}","BR_{ppjdp}","BR_{c2jpsi}","BR_{c1jpsi}","BR_{jpsidm}","BR_{b2ups1}","BR_{b1ups1}","BR_{ups3dm}","BR_{ups2dm}","BR_{ups1dm}","L_{CMS}","L_{ATLAS}","L_{ATLAS}(\\Upsilon)","L(13/7)_{c\\overline c}","L(13/7)_{b\\overline b}","L_{CMS}(13)","\\Delta M_{c\\overline c}","\\Delta M_{b\\overline b}"};
  
  tex.open("plots/fitp.tex");

  tex << "\\begin{table}" << endl;
  tex << "\\centering" << endl;
  tex << "\\begin{tabular}{c|c|c}" << endl;
  tex << "Parameter & Value & Uncertainty \\\\" << endl;
  tex << "\\hline" << endl;
  for(int i=0; i<36; i++)
    {
      tex << "$" << names[i] << "$ & " << param[i] << " & ";
      if(eparam[i]==0) tex << "fixed \\\\" << endl;
      else tex << eparam[i] << " \\\\" << endl;
    }
  tex << "min $p_T/M$ & " << ptmmin << " & fixed \\\\" << endl;
  tex << "\\end{tabular}" << endl;
  tex << "\\caption{Fit parameters ($\\chi^2$ / ndf = " << minimum << " / " << ndf << " = " << minimum/ndf << ")}" << endl;
  tex << "\\end{table}" << endl;
  tex.close();
  
  //part: plotting data
  float xi = -2, xf = 49.9, yi = 1.01e-5, yf = 9.99e2;
  
  //plot of psiprime results
  //cross section plot
  //CMS data
  const int npsipccs=ns[0];
  float datapsipccs[5][npsipccs];
  for(int i=0; i<npsipccs; i++)
    {
      datapsipccs[0][i]=avgptm(datasigma[5][i], datasigma[6][i], datasigma[7][i], datasigma[8][i], datasigma[9][i]/(param[34]+mass[0]), param[9], param[10], param[11], param[12], param[13], param[0], param[34]+mass[0], param[14], param[15], param[16], param[17], 0);
      datapsipccs[1][i]=datasigma[1][i]*param[28]*param[18];
      datapsipccs[2][i]=datasigma[2][i]*param[28]*param[18];
      datapsipccs[3][i]=datapsipccs[0][i]-datasigma[5][i];
      datapsipccs[4][i]=datasigma[6][i]-datapsipccs[0][i];
    }

  counter+=npsipccs;
  //ATLAS data
  const int npsipacs=ns[1];
  float datapsipacs[5][npsipacs];
  for(int i=0; i<npsipacs; i++)
    {
      datapsipacs[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter], datasigma[8][i+counter], datasigma[9][i+counter]/(param[34]+mass[0]), param[9], param[10], param[11], param[12], param[13], param[0], param[34]+mass[0], param[14], param[15], param[16], param[17], 0);
      datapsipacs[1][i]=datasigma[1][i+counter]*param[29]*param[19]*param[22];
      datapsipacs[2][i]=datasigma[2][i+counter]*param[29]*param[19]*param[22];
      datapsipacs[3][i]=datapsipacs[0][i]-datasigma[5][i+counter];
      datapsipacs[4][i]=datasigma[6][i+counter]-datapsipacs[0][i];
    }

  TCanvas *psipcs = new TCanvas("psip cross section", "psip cross section", 700, 700);
  psipcs->SetLogy();
  
  TH1F *fpsipcs = psipcs->DrawFrame(xi, yi, xf, yf);
  fpsipcs->SetXTitle("p_{T}/M");
  fpsipcs->SetYTitle("d#sigma / d#xidy (nb/GeV)");
  fpsipcs->GetYaxis()->SetTitleOffset(1);
  psipcs->Modified();
  psipcs->SetTitle("");

  //plot CMS data
  TGraphAsymmErrors *pcsection = new TGraphAsymmErrors(npsipccs, datapsipccs[0], datapsipccs[1], datapsipccs[3], datapsipccs[4], datapsipccs[2], datapsipccs[2]);
  pcsection->SetMarkerStyle(20);
  pcsection->SetLineColor(kBlack);
  pcsection->SetMarkerColor(kBlack);
  pcsection->SetMarkerSize(.75);
  pcsection->Draw("P");

  //plot ATLAS data
  TGraphAsymmErrors *pasection = new TGraphAsymmErrors(npsipacs, datapsipacs[0], datapsipacs[1], datapsipacs[3], datapsipacs[4], datapsipacs[2], datapsipacs[2]);
  pasection->SetMarkerStyle(25);
  pasection->SetMarkerColor(kBlack);
  pasection->SetLineColor(kBlack);
  pasection->SetMarkerSize(.75);
  pasection->Draw("P");

  //plot the fitted function and each contribution
  TF1 *fitpsipcs = new TF1("psip cs fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitpsipcs->SetParameter(0, datasigma[9][0]/(param[34]+mass[0]));
  fitpsipcs->SetParameter(1, datasigma[8][0]/2);
  for(int i=2; i<7; i++)
    fitpsipcs->SetParameter(i, param[i+7]);
  fitpsipcs->SetParameter(7, param[0]);
  fitpsipcs->SetParameter(8, param[34]+mass[0]);
  for(int i=9; i<13; i++)
    fitpsipcs->SetParameter(i, param[i+5]);
  fitpsipcs->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitpsipcs->SetParameter(14, param[34]+mass[0]);
  fitpsipcs->SetLineColor(kBlue);
  fitpsipcs->Draw("lsame");
  
  //text on the plot
  TLatex lpsipcs;
  lpsipcs.SetTextSize(0.03);
  lpsipcs.DrawLatex(8, 20, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpsipcs.DrawLatex(8, 7, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpsipcs.DrawLatex(0, 2.e-5, "pp 7 TeV");

  //draw legend
  TLegend *legpsipcs = new TLegend(0.6, 0.7, 0.9, 0.9);
  legpsipcs->SetTextSize(0.03);
  legpsipcs->AddEntry(pcsection, "CMS #psi(2S)", "p");
  legpsipcs->AddEntry(pasection, "ATLAS #psi(2S)", "p");
  legpsipcs->AddEntry(fitpsipcs, "model", "l");
  legpsipcs->Draw();
  //save psiprime cross section plot
  psipcs->SaveAs("plots/psiprime_cs.pdf");
   
  //plot of chic2 results
  //cross section plot
  //ATLAS data
  counter+=npsipacs;
  const int nchic2cs=ns[2];
  float datachic2cs[5][nchic2cs];  
  for(int i=0; i<nchic2cs; i++)
    {
      datachic2cs[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter], datasigma[8][i+counter], datasigma[9][i+counter]/(param[34]+mass[1]), param[9], param[10], param[11], param[12], param[13], param[1], param[34]+mass[1], param[14], param[15], param[16], param[17], 1);
      datachic2cs[1][i]=datasigma[1][i+counter]*param[29]*param[20]*param[22];
      datachic2cs[2][i]=datasigma[2][i+counter]*param[29]*param[20]*param[22];
      datachic2cs[3][i]=datachic2cs[0][i]-datasigma[5][i+counter];
      datachic2cs[4][i]=datasigma[6][i+counter]-datachic2cs[0][i];
    }

  TCanvas *chic2cs = new TCanvas("chic2 cross section", "chic2 cross section", 700, 700);
  chic2cs->SetLogy();
  TH1F *fchic2cs = chic2cs->DrawFrame(xi, yi, xf, yf); 
  fchic2cs->SetXTitle("p_{T}/M");
  fchic2cs->SetYTitle("d#sigma / d#xidy (nb/GeV)");
  chic2cs->Modified();
  chic2cs->SetTitle("");

  //plot data
  TGraphAsymmErrors *c2section = new TGraphAsymmErrors(nchic2cs, datachic2cs[0], datachic2cs[1], datachic2cs[3], datachic2cs[4], datachic2cs[2], datachic2cs[2]);
  c2section->SetMarkerStyle(20);
  c2section->SetLineColor(kBlack);
  c2section->SetMarkerColor(kBlack);
  c2section->SetMarkerSize(.75);
  c2section->Draw("P");
  
  //plot the fitted function and the direct contributions
  TF1 *fitchic2cs = new TF1("psip cs fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitchic2cs->SetParameter(0, datasigma[9][counter]/(param[34]+mass[1]));
  fitchic2cs->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitchic2cs->SetParameter(i, param[i+7]);
  fitchic2cs->SetParameter(7, param[1]);
  fitchic2cs->SetParameter(8, param[34]+mass[1]);
  for(int i=9; i<13; i++)
    fitchic2cs->SetParameter(i, param[i+5]);
  fitchic2cs->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitchic2cs->SetParameter(14, param[34]+mass[0]);
  fitchic2cs->SetLineColor(kBlue);
  fitchic2cs->Draw("lsame");
 
  //text on the plot
  TLatex lchic2cs;
  lchic2cs.SetTextSize(0.03);
  lchic2cs.DrawLatex(8, 50, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lchic2cs.DrawLatex(8, 20, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lchic2cs.DrawLatex(0, 2.e-5, "pp 7 TeV");
  
  //draw legend
  TLegend *legchic2cs = new TLegend(0.6, 0.75, 0.9, 0.9);
  legchic2cs->SetTextSize(0.03);
  legchic2cs->AddEntry(c2section, "ATLAS #chi_{c2}", "p");
  legchic2cs->AddEntry(fitchic2cs, "model", "l");
  legchic2cs->Draw();
  //save chic2 cross section plot 
  chic2cs->SaveAs("plots/chic2_cs.pdf");
  
  //plot of chic1 results
  //cross section plot
  //ATLAS data
  counter+=nchic2cs;
  const int nchic1cs=ns[3];
  float datachic1cs[5][nchic1cs];  
  for(int i=0; i<nchic1cs; i++)
    {
      datachic1cs[0][i] = avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter], datasigma[8][i+counter], datasigma[9][i+counter]/(param[34]+mass[2]), param[9], param[10], param[11], param[12], param[13], param[2], param[34]+mass[2], param[14], param[15], param[16], param[17], 2);
      datachic1cs[1][i]=datasigma[1][i+counter]*param[29]*param[21]*param[22];
      datachic1cs[2][i]=datasigma[2][i+counter]*param[29]*param[21]*param[22];
      datachic1cs[3][i]=datachic1cs[0][i]-datasigma[5][i+counter];
      datachic1cs[4][i]=datasigma[6][i+counter]-datachic1cs[0][i];
    }
  
  TCanvas *chic1cs = new TCanvas("chic1 cross section", "chic1 cross section", 700, 700);
  chic1cs->SetLogy();
  TH1F *fchic1cs = chic1cs->DrawFrame(xi, yi, xf, yf);
  fchic1cs->SetXTitle("p_{T}/M");
  fchic1cs->SetYTitle("d#sigma / d#xidy (nb/GeV)");
  chic1cs->Modified();
  chic1cs->SetTitle("");

  //plot data
  TGraphAsymmErrors *c1section = new TGraphAsymmErrors(nchic1cs, datachic1cs[0], datachic1cs[1], datachic1cs[3], datachic1cs[4], datachic1cs[2], datachic1cs[2]);
  c1section->SetMarkerStyle(20);
  c1section->SetLineColor(kBlack);
  c1section->SetMarkerColor(kBlack);
  c1section->SetMarkerSize(.75);
  c1section->Draw("P");
  
  //plot the fitted function and the direct contributions
  TF1 *fitchic1cs = new TF1("psip cs fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitchic1cs->SetParameter(0, datasigma[9][counter]/(param[34]+mass[2]));
  fitchic1cs->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitchic1cs->SetParameter(i, param[i+7]);
  fitchic1cs->SetParameter(7, param[2]);
  fitchic1cs->SetParameter(8, param[34]+mass[2]);
  for(int i=9; i<13; i++)
    fitchic1cs->SetParameter(i, param[i+5]);
  fitchic1cs->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitchic1cs->SetParameter(14, param[34]+mass[0]);
  fitchic1cs->SetLineColor(kBlue);
  fitchic1cs->Draw("lsame");

  //text on the plot
  TLatex lchic1cs;
  lchic1cs.SetTextSize(0.03);
  lchic1cs.DrawLatex(8, 50, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lchic1cs.DrawLatex(8, 20, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lchic1cs.DrawLatex(0, 2.e-5, "pp 7 TeV");

  //draw legend
  TLegend *legchic1cs = new TLegend(0.6, 0.75, 0.9, 0.9);
  legchic1cs->SetTextSize(0.03);
  legchic1cs->AddEntry(c1section, "ATLAS #chi_{c1}", "p");
  legchic1cs->AddEntry(fitchic1cs, "model", "l");
  legchic1cs->Draw();
  //save chic1 cross section plot
  chic1cs->SaveAs("plots/chic1_cs.pdf");
  
  //plot of j/psi results
  //cross section plot
  //CMS data
  counter+=nchic1cs;
  const int njpsics=ns[4];
  float datajpsics[5][njpsics];
  for(int i=0; i<njpsics; i++)
    {
      datajpsics[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter], datasigma[8][i+counter], datasigma[9][i+counter]/(param[34]+mass[3]), param[9], param[10], param[11], param[12], param[13], param[3], param[34]+mass[3], param[14], param[15], param[16], param[17], 3);
      datajpsics[1][i]=datasigma[1][i+counter]*param[28]*param[22];
      datajpsics[2][i]=datasigma[2][i+counter]*param[28]*param[22];
      datajpsics[3][i]=datajpsics[0][i]-datasigma[5][i+counter];
      datajpsics[4][i]=datasigma[6][i+counter]-datajpsics[0][i];
    }
  
  TCanvas *jpsics = new TCanvas("jpsi cross section", "jpsi cross section", 700, 700);
  jpsics->SetLogy();
  TH1F *fjpsics = jpsics->DrawFrame(xi, yi, xf, yf);
  fjpsics->SetXTitle("p_{T}/M");
  fjpsics->SetYTitle("d#sigma / d#xidy (nb/GeV)");
  jpsics->Modified();
  jpsics->SetTitle("");
  
  //plot data
  TGraphAsymmErrors *jsection = new TGraphAsymmErrors(njpsics, datajpsics[0], datajpsics[1], datajpsics[3], datajpsics[4], datajpsics[2], datajpsics[2]);
  jsection->SetMarkerStyle(20);
  jsection->SetLineColor(kBlack);
  jsection->SetMarkerColor(kBlack);
  jsection->SetMarkerSize(.75);
  jsection->Draw("P");

  //plot the fitted function
  TF1 *fitjpsics = new TF1("jpsi cs fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitjpsics->SetParameter(0, datasigma[9][counter]/(param[34]+mass[3]));
  fitjpsics->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitjpsics->SetParameter(i, param[i+7]);
  fitjpsics->SetParameter(7, param[3]);
  fitjpsics->SetParameter(8, param[34]+mass[3]);
  for(int i=9; i<13; i++)
    fitjpsics->SetParameter(i, param[i+5]);
  fitjpsics->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitjpsics->SetParameter(14, param[34]+mass[0]);
  fitjpsics->SetLineColor(kBlue);
  fitjpsics->Draw("lsame");
  
  //text on the plot
  TLatex ljpsics;
  ljpsics.SetTextSize(0.03);
  ljpsics.DrawLatex(8, 20, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  ljpsics.DrawLatex(8, 7, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  ljpsics.DrawLatex(0, 2.e-5, "pp 7 TeV");

  //draw legend
  TLegend *legjpsics = new TLegend(0.6, 0.75, 0.9, 0.9);
  legjpsics->SetTextSize(0.03);
  legjpsics->AddEntry(jsection, "CMS J/#psi", "p");
  legjpsics->AddEntry(fitjpsics, "model", "l");
  legjpsics->Draw();
  //save jpsi cross section plot
  jpsics->SaveAs("plots/jpsi_cs.pdf");

  //plot of chic2/chic1 ratio
  counter+=njpsics;
  const int nchicr=ns[5];
  float datachicr[5][nchicr];  
  for(int i=0; i<nchicr; i++)
    {
      datachicr[0][i]=datasigma[0][i+counter];
      datachicr[1][i]=datasigma[1][i+counter]*param[20]/param[21];
      datachicr[2][i]=datasigma[2][i+counter]*param[20]/param[21];
      datachicr[3][i]=datachicr[0][i]-datasigma[5][i+counter];
      datachicr[4][i]=datasigma[6][i+counter]-datachicr[0][i];
      }

  TCanvas *chicr = new TCanvas("chic2/chic1 cross section ratio", "chic2/chic1 cross section ratio", 700, 700);
  TH1F *fchicr=chicr->DrawFrame(0.01, 0, 9, 1.39);
  fchicr->SetXTitle("p_{T}/M");
  fchicr->SetYTitle("ratio");
  chicr->Modified();
  chicr->SetTitle("");

  //plot data
  TGraphAsymmErrors *cratio = new TGraphAsymmErrors(nchicr, datachicr[0], datachicr[1], datachicr[3], datachicr[4], datachicr[2], datachicr[2]);
  cratio->SetMarkerStyle(20);
  cratio->SetLineColor(kBlack);
  cratio->SetMarkerColor(kBlack);
  cratio->SetMarkerSize(.75);
  cratio->Draw("P");

  //plot the fitted function
  TF1 *fitchicr = new TF1("psip cs fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [13], [9], [10], [11], [12])/sigplot([0], x, [1], [2], [3], [4], [5], [6], [8], [13], [9], [10], [11], [12])", 0, 10);
  fitchicr->SetParameter(0, datasigma[9][counter]/(param[34]+mass[3]));
  fitchicr->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitchicr->SetParameter(i, param[i+7]);
  fitchicr->SetParameter(7, param[1]);
  fitchicr->SetParameter(8, param[2]);
  for(int i=9; i<13; i++)
    fitchicr->SetParameter(i, param[i+5]);
  fitchicr->SetParameter(13, param[34]+mass[3]);
  fitchicr->SetLineColor(kBlue);
  fitchicr->Draw("lsame");
  
  //text on the plot
  TLatex lchicr;
  lchicr.SetTextSize(0.03);
  lchicr.DrawLatex(0.5, 1.3, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lchicr.DrawLatex(0.5, 1.2, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lchicr.DrawLatex(0.25, 0.05, "pp 7 TeV");

  //draw legend
  TLegend *legchicr = new TLegend(0.6, 0.75, 0.9, 0.9);
  legchicr->SetTextSize(0.03);
  legchicr->AddEntry(cratio, "CMS #chi_{c2} / #chi_{c1}", "p");
  legchicr->AddEntry(fitchicr, "model", "l");
  legchicr->Draw();
  //save chic ratio plot
  chicr->SaveAs("plots/chic_ratio.pdf");

  //plot of chib2/chib1 ratio
  counter+=nchicr;
  const int nchibr=ns[6];
  float datachibr[5][nchibr];  
  for(int i=0; i<nchibr; i++)
    {
      datachibr[0][i]=datasigma[0][i+counter];
      datachibr[1][i]=datasigma[1][i+counter]*param[23]/param[24];
      datachibr[2][i]=datasigma[2][i+counter]*param[23]/param[24];
      datachibr[3][i]=datachibr[0][i]-datasigma[5][i+counter];
      datachibr[4][i]=datasigma[6][i+counter]-datachibr[0][i];
      }

  TCanvas *chibr = new TCanvas("chib2/chib1 cross section ratio", "chib2/chib1 cross section ratio", 700, 700);
  TH1F *fchibr=chibr->DrawFrame(0.01, 0, 9, 1.39);
  fchibr->SetXTitle("p_{T}/M");
  fchibr->SetYTitle("ratio");
  chibr->Modified();
  chibr->SetTitle("");

  //plot data
  TGraphAsymmErrors *bratio = new TGraphAsymmErrors(nchibr, datachibr[0], datachibr[1], datachibr[3], datachibr[4], datachibr[2], datachibr[2]);
  bratio->SetMarkerStyle(20);
  bratio->SetLineColor(kBlack);
  bratio->SetMarkerColor(kBlack);
  bratio->SetMarkerSize(.75);
  bratio->Draw("P");

  //plot the fitted function
  TF1 *fitchibr = new TF1("chib2/chib1 ratio fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [13], [9], [10], [11], [12])/sigplot([0], x, [1], [2], [3], [4], [5], [6], [8], [13], [9], [10], [11], [12])", 0, 10);
  fitchibr->SetParameter(0, datasigma[9][counter]/(param[35]+mass[6]));
  fitchibr->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitchibr->SetParameter(i, param[i+7]);
  fitchibr->SetParameter(7, param[6]);
  fitchibr->SetParameter(8, param[7]);
  for(int i=9; i<13; i++)
    fitchibr->SetParameter(i, param[i+5]);
  fitchibr->SetParameter(13, param[35]+mass[6]);
  fitchibr->SetLineColor(kBlue);
  fitchibr->Draw("lsame");
  
  //text on the plot
  TLatex lchibr;
  lchibr.SetTextSize(0.03);
  lchibr.DrawLatex(0.5, 1.3, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lchibr.DrawLatex(0.5, 1.2, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lchibr.DrawLatex(0.25, 0.05, "pp 7 TeV");

  //draw legend
  TLegend *legchibr = new TLegend(0.6, 0.75, 0.9, 0.9);
  legchibr->SetTextSize(0.03);
  legchibr->AddEntry(bratio, "CMS #chi_{b2} / #chi_{b1}", "p");
  legchibr->AddEntry(fitchibr, "model", "l");
  legchibr->Draw();
  //save chib ratio plot
  chibr->SaveAs("plots/chib_ratio.pdf");

  //CMS ups(3S)
  const int nups3ccs=ns[7];
  counter+=nchibr;
  float dataups3ccs[5][nups3ccs];
  for(int i=0; i<nups3ccs; i++)
    {
      dataups3ccs[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter], datasigma[8][i+counter], datasigma[9][i+counter]/(param[35]+mass[4]), param[9], param[10], param[11], param[12], param[13], param[4], param[35]+mass[4], param[14], param[15], param[16], param[17], 4);
      dataups3ccs[1][i]=datasigma[1][i+counter]*param[28]*param[25];
      dataups3ccs[2][i]=datasigma[2][i+counter]*param[28]*param[25];
      dataups3ccs[3][i]=dataups3ccs[0][i]-datasigma[5][i+counter];
      dataups3ccs[4][i]=datasigma[6][i+counter]-dataups3ccs[0][i];
    }
  
  counter+=nups3ccs;
  //ATLAS ups(3S)
  const int nups3acs=ns[8];
  float dataups3acs[5][nups3acs];
  for(int i=0; i<nups3acs; i++)
    {
      dataups3acs[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter], datasigma[8][i+counter], datasigma[9][i+counter]/(param[35]+mass[4]), param[9], param[10], param[11], param[12], param[13], param[4], param[35]+mass[4], param[14], param[15], param[16], param[17], 4);
      dataups3acs[1][i]=datasigma[1][i+counter]*param[25]*param[30];
      dataups3acs[2][i]=datasigma[2][i+counter]*param[25]*param[30];
      dataups3acs[3][i]=dataups3acs[0][i]-datasigma[5][i+counter];
      dataups3acs[4][i]=datasigma[6][i+counter]-dataups3acs[0][i];
    }

  TCanvas *ups3cs = new TCanvas("ups3 cross section", "ups3 cross section", 700, 700);
  ups3cs->SetLogy();
  
  TH1F *fups3cs = ups3cs->DrawFrame(xi, yi, xf, yf);
  fups3cs->SetXTitle("p_{T}/M");
  fups3cs->SetYTitle("d#sigma / dp_{T} (nb/GeV)");
  ups3cs->Modified();
  ups3cs->SetTitle("");

  //plot CMS data
  TGraphAsymmErrors *u3csection = new TGraphAsymmErrors(nups3ccs, dataups3ccs[0], dataups3ccs[1], dataups3ccs[3], dataups3ccs[4], dataups3ccs[2], dataups3ccs[2]);
  u3csection->SetMarkerStyle(20);
  u3csection->SetLineColor(kBlack);
  u3csection->SetMarkerColor(kBlack);
  u3csection->SetMarkerSize(.75);
  u3csection->Draw("P");

  //plot ATLAS data
  TGraphAsymmErrors *u3asection = new TGraphAsymmErrors(nups3acs, dataups3acs[0], dataups3acs[1], dataups3acs[3], dataups3acs[4], dataups3acs[2], dataups3acs[2]);
  u3asection->SetMarkerStyle(25);
  u3asection->SetMarkerColor(kBlack);
  u3asection->SetLineColor(kBlack);
  u3asection->SetMarkerSize(.75);
  u3asection->Draw("P");
  
  //plot the fitted function and each contribution
  TF1 *fitups3cs = new TF1("ups3 cs fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitups3cs->SetParameter(0, datasigma[9][counter]/(param[35]+mass[4]));
  fitups3cs->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitups3cs->SetParameter(i, param[i+7]);
  fitups3cs->SetParameter(7, param[4]);
  fitups3cs->SetParameter(8, param[35]+mass[4]);
  for(int i=9; i<13; i++)
    fitups3cs->SetParameter(i, param[i+5]);
  fitups3cs->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitups3cs->SetParameter(14, param[34]+mass[0]);
  fitups3cs->SetLineColor(kBlue);
  fitups3cs->Draw("lsame");
  
  //text on the plot
  TLatex lups3cs;
  lups3cs.SetTextSize(0.03);
  lups3cs.DrawLatex(8, 50, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lups3cs.DrawLatex(8, 20, Form("P(#chi^{2},ndf) = %.2g%%", 100*chiprob));
  lups3cs.DrawLatex(0, 2.e-5, "pp 7 TeV");

  //draw legend
  TLegend *legups3cs = new TLegend(0.6, 0.7, 0.9, 0.9);
  legups3cs->SetTextSize(0.03);
  legups3cs->AddEntry(u3csection, "CMS #Upsilon(3S)", "p");
  legups3cs->AddEntry(u3asection, "ATLAS #Upsilon(3S)", "p");
  legups3cs->AddEntry(fitups3cs, "model", "l");
  legups3cs->Draw();
  //save ups(3S) cross section plot
  ups3cs->SaveAs("plots/ups3S_cs.pdf");

  //CMS ups(2S)
  const int nups2ccs=ns[9];
  counter+=nups3acs;
  float dataups2ccs[5][nups2ccs];
  for(int i=0; i<nups2ccs; i++)
    {
      dataups2ccs[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter],datasigma[8][i+counter], datasigma[9][i+counter]/(param[35]+mass[5]), param[9], param[10], param[11], param[12], param[13], param[5], param[35]+mass[5], param[14], param[15], param[16], param[17], 5);
      dataups2ccs[1][i]=datasigma[1][i+counter]*param[28]*param[26];
      dataups2ccs[2][i]=datasigma[2][i+counter]*param[28]*param[26];
      dataups2ccs[3][i]=dataups2ccs[0][i]-datasigma[5][i+counter];
      dataups2ccs[4][i]=datasigma[6][i+counter]-dataups2ccs[0][i];
    }
  
  counter+=nups2ccs;
  //ATLAS ups(2S)
  const int nups2acs=ns[10];
  float dataups2acs[5][nups2acs];
  for(int i=0; i<nups2acs; i++)
    {
      dataups2acs[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter],datasigma[8][i+counter], datasigma[9][i+counter]/(param[35]+mass[5]), param[9], param[10], param[11], param[12], param[13], param[5], param[35]+mass[5], param[14], param[15], param[16], param[17], 5);
      dataups2acs[1][i]=datasigma[1][i+counter]*param[30]*param[26];
      dataups2acs[2][i]=datasigma[2][i+counter]*param[30]*param[26];
      dataups2acs[3][i]=dataups2acs[0][i]-datasigma[5][i+counter];
      dataups2acs[4][i]=datasigma[6][i+counter]-dataups2acs[0][i];
    }

  TCanvas *ups2cs = new TCanvas("ups2 cross section", "ups2 cross section", 700, 700);
  ups2cs->SetLogy();
  
  TH1F *fups2cs = ups2cs->DrawFrame(xi, yi, xf, yf);
  fups2cs->SetXTitle("p_{T}/M");
  fups2cs->SetYTitle("d#sigma / dp_{T} (nb/GeV)");
  ups2cs->Modified();
  ups2cs->SetTitle("");

  //plot CMS data
  TGraphAsymmErrors *u2csection = new TGraphAsymmErrors(nups2ccs, dataups2ccs[0], dataups2ccs[1], dataups2ccs[3], dataups2ccs[4], dataups2ccs[2], dataups2ccs[2]);
  u2csection->SetMarkerStyle(20);
  u2csection->SetLineColor(kBlack);
  u2csection->SetMarkerColor(kBlack);
  u2csection->SetMarkerSize(.75);
  u2csection->Draw("P");

  //plot ATLAS data
  TGraphAsymmErrors *u2asection = new TGraphAsymmErrors(nups2acs, dataups2acs[0], dataups2acs[1], dataups2acs[3], dataups2acs[4], dataups2acs[2], dataups2acs[2]);
  u2asection->SetMarkerStyle(25);
  u2asection->SetMarkerColor(kBlack);
  u2asection->SetLineColor(kBlack);
  u2asection->SetMarkerSize(.75);
  u2asection->Draw("P");

  //plot the fitted function and each contribution
  TF1 *fitups2cs = new TF1("ups2 cs fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitups2cs->SetParameter(0, datasigma[9][counter]/(param[35]+mass[5]));
  fitups2cs->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitups2cs->SetParameter(i, param[i+7]);
  fitups2cs->SetParameter(7, param[5]);
  fitups2cs->SetParameter(8, param[35]+mass[5]);
  for(int i=9; i<13; i++)
    fitups2cs->SetParameter(i, param[i+5]);
  fitups2cs->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitups2cs->SetParameter(14, param[34]+mass[0]);
  fitups2cs->SetLineColor(kBlue);
  fitups2cs->Draw("lsame");

  //text on the plot
  TLatex lups2cs;
  lups2cs.SetTextSize(0.03);
  lups2cs.DrawLatex(8, 50, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lups2cs.DrawLatex(8, 20, Form("P(#chi^{2},ndf) = %.2g%%", 100*chiprob));
  lups2cs.DrawLatex(0, 2.e-5, "pp 7 TeV");

  //draw legend
  TLegend *legups2cs = new TLegend(0.6, 0.7, 0.9, 0.9);
  legups2cs->SetTextSize(0.03);
  legups2cs->AddEntry(u2csection, "CMS #Upsilon(2S)", "p");
  legups2cs->AddEntry(u2asection, "ATLAS #Upsilon(2S)", "p");
  legups2cs->AddEntry(fitups2cs, "model", "l");
  legups2cs->Draw();
  //save ups(2S) cross section plot
  ups2cs->SaveAs("plots/ups2S_cs.pdf");

  //CMS ups(1S)
  const int nups1ccs=ns[11];
  counter+=nups2acs;
  float dataups1ccs[5][nups1ccs];
  for(int i=0; i<nups1ccs; i++)
    {
      dataups1ccs[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter],datasigma[8][i+counter], datasigma[9][i+counter]/(param[35]+mass[6]), param[9], param[10], param[11], param[12], param[13], param[8], param[35]+mass[6], param[14], param[15], param[16], param[17], 6);
      dataups1ccs[1][i]=datasigma[1][i+counter]*param[28]*param[27];
      dataups1ccs[2][i]=datasigma[2][i+counter]*param[28]*param[27];
      dataups1ccs[3][i]=dataups1ccs[0][i]-datasigma[5][i+counter];
      dataups1ccs[4][i]=datasigma[6][i+counter]-dataups1ccs[0][i];
    }
  
  counter+=nups1ccs;
  //ATLAS ups(1S)
  const int nups1acs=ns[12];
  float dataups1acs[5][nups1acs];
  for(int i=0; i<nups1acs; i++)
    {
      dataups1acs[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter],datasigma[8][i+counter], datasigma[9][i+counter]/(param[35]+mass[6]), param[9], param[10], param[11], param[12], param[13], param[8], param[35]+mass[6], param[14], param[15], param[16], param[17], 6);
      dataups1acs[1][i]=datasigma[1][i+counter]*param[30]*param[27];
      dataups1acs[2][i]=datasigma[2][i+counter]*param[30]*param[27];
      dataups1acs[3][i]=dataups1acs[0][i]-datasigma[5][i+counter];
      dataups1acs[4][i]=datasigma[6][i+counter]-dataups1acs[0][i];
    }
  
  TCanvas *ups1cs = new TCanvas("ups1 cross section", "ups1 cross section", 700, 700);
  ups1cs->SetLogy();
  
  TH1F *fups1cs = ups1cs->DrawFrame(xi, yi, xf, yf);
  fups1cs->SetXTitle("p_{T}/M");
  fups1cs->SetYTitle("d#sigma / dp_{T} (nb/GeV)");
  ups1cs->Modified();
  ups1cs->SetTitle("");

  //plot CMS data
  TGraphAsymmErrors *u1csection = new TGraphAsymmErrors(nups1ccs, dataups1ccs[0], dataups1ccs[1], dataups1ccs[3], dataups1ccs[4], dataups1ccs[2], dataups1ccs[2]);
  u1csection->SetMarkerStyle(20);
  u1csection->SetLineColor(kBlack);
  u1csection->SetMarkerColor(kBlack);
  u1csection->SetMarkerSize(.75);
  u1csection->Draw("P");

  //plot ATLAS data
  TGraphAsymmErrors *u1asection = new TGraphAsymmErrors(nups1acs, dataups1acs[0], dataups1acs[1], dataups1acs[3], dataups1acs[4], dataups1acs[2], dataups1acs[2]);
  u1asection->SetMarkerStyle(25);
  u1asection->SetMarkerColor(kBlack);
  u1asection->SetLineColor(kBlack);
  u1asection->SetMarkerSize(.75);
  u1asection->Draw("P");
  
  //plot the fitted function and each contribution
  TF1 *fitups1cs = new TF1("ups1 cs fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitups1cs->SetParameter(0, datasigma[9][counter]/(param[35]+mass[6]));
  fitups1cs->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitups1cs->SetParameter(i, param[i+7]);
  fitups1cs->SetParameter(7, param[8]);
  fitups1cs->SetParameter(8, param[35]+mass[6]);
  for(int i=9; i<13; i++)
    fitups1cs->SetParameter(i, param[i+5]);
  fitups1cs->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitups1cs->SetParameter(14, param[34]+mass[0]);
  fitups1cs->SetLineColor(kBlue);
  fitups1cs->Draw("lsame");

  //text on the plot
  TLatex lups1cs;
  lups1cs.SetTextSize(0.03);
  lups1cs.DrawLatex(8, 50, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lups1cs.DrawLatex(8, 20, Form("P(#chi^{2},ndf) = %.2g%%", 100*chiprob));
  lups1cs.DrawLatex(0, 2.e-5, "pp 7 TeV");

  //draw legend
  TLegend *legups1cs = new TLegend(0.6, 0.7, 0.9, 0.9);
  legups1cs->SetTextSize(0.03);
  legups1cs->AddEntry(u1csection, "CMS #Upsilon(1S)", "p");
  legups1cs->AddEntry(u1asection, "ATLAS #Upsilon(1S)", "p");
  legups1cs->AddEntry(fitups1cs, "model", "l");
  legups1cs->Draw();
  //save ups(1S) cross section plot
  ups1cs->SaveAs("plots/ups1S_cs.pdf");

  //13 TeV psiprime
  counter+=nups1acs;
  const int npsipccs13=ns[13];
  float datapsipccs13[5][npsipccs13];
  for(int i=0; i<npsipccs13; i++)
    {
      datapsipccs13[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter], datasigma[8][i+counter], datasigma[9][i+counter]/(param[34]+mass[0]), param[9], param[10], param[11], param[12], param[13], param[0], param[34]+mass[0], param[14], param[15], param[16], param[17], 0);
      datapsipccs13[1][i]=datasigma[1][i+counter]*param[33]*param[18];
      datapsipccs13[2][i]=datasigma[2][i+counter]*param[33]*param[18];
      datapsipccs13[3][i]=datapsipccs13[0][i]-datasigma[5][i+counter];
      datapsipccs13[4][i]=datasigma[6][i+counter]-datapsipccs13[0][i];
    }
  
  TCanvas *psipcs13 = new TCanvas("psip13 cross section", "psip13 cross section", 700, 700);
  psipcs13->SetLogy();
  
  TH1F *fpsipcs13 = psipcs13->DrawFrame(xi, yi, xf, yf);
  fpsipcs13->SetXTitle("p_{T}/M");
  fpsipcs13->SetYTitle("d#sigma / d#xidy (nb/GeV)");
  fpsipcs13->GetYaxis()->SetTitleOffset(1);
  psipcs13->Modified();
  psipcs13->SetTitle("");

  //plot CMS data
  TGraphAsymmErrors *pcsection13 = new TGraphAsymmErrors(npsipccs13, datapsipccs13[0], datapsipccs13[1], datapsipccs13[3], datapsipccs13[4], datapsipccs13[2], datapsipccs13[2]);
  pcsection13->SetMarkerStyle(20);
  pcsection13->SetLineColor(kBlack);
  pcsection13->SetMarkerColor(kBlack);
  pcsection13->SetMarkerSize(.75);
  pcsection13->Draw("P");

  //plot the fitted function and each contribution
  TF1 *fitpsipcs13 = new TF1("psip13 cs fit", "[15]*sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitpsipcs13->SetParameter(0, datasigma[9][counter]/(param[34]+mass[0]));
  fitpsipcs13->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitpsipcs13->SetParameter(i, param[i+7]);
  fitpsipcs13->SetParameter(7, param[0]);
  fitpsipcs13->SetParameter(8, param[34]+mass[0]);
  for(int i=9; i<13; i++)
    fitpsipcs13->SetParameter(i, param[i+5]);
  fitpsipcs13->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitpsipcs13->SetParameter(14, param[34]+mass[0]);
  fitpsipcs13->SetParameter(15, param[31]);
  fitpsipcs13->SetLineColor(kBlue);
  fitpsipcs13->Draw("lsame");

  //text on the plot
  TLatex lpsipcs13;
  lpsipcs13.SetTextSize(0.03);
  lpsipcs13.DrawLatex(8, 20, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lpsipcs13.DrawLatex(8, 7, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  lpsipcs13.DrawLatex(0, 2.e-5, "pp 13 TeV");

  //draw legend
  TLegend *legpsipcs13 = new TLegend(0.6, 0.75, 0.9, 0.9);
  legpsipcs13->SetTextSize(0.03);
  legpsipcs13->AddEntry(pcsection13, "CMS #psi(2S)", "p");
  legpsipcs13->AddEntry(fitpsipcs13, "model", "l");
  legpsipcs13->Draw();
  //save psiprime cross section plot
  psipcs13->SaveAs("plots/psiprime_cs_13.pdf");
  
  //plot of j/psi results
  //cross section plot
  //CMS data
  counter+=npsipccs13;
  const int njpsics13=ns[14];
  float datajpsics13[5][njpsics13];
  for(int i=0; i<njpsics13; i++)
    {
      datajpsics13[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter], datasigma[8][i+counter], datasigma[9][i+counter]/(param[34]+mass[3]), param[9], param[10], param[11], param[12], param[13], param[3], param[34]+mass[3], param[14], param[15], param[16], param[17], 3);
      datajpsics13[1][i]=datasigma[1][i+counter]*param[33]*param[22];
      datajpsics13[2][i]=datasigma[2][i+counter]*param[33]*param[22];
      datajpsics13[3][i]=datajpsics13[0][i]-datasigma[5][i+counter];
      datajpsics13[4][i]=datasigma[6][i+counter]-datajpsics13[0][i];
    }
  
  TCanvas *jpsics13 = new TCanvas("jpsi13 cross section", "jpsi13 cross section", 700, 700);
  jpsics13->SetLogy();
  
  TH1F *fjpsics13 = jpsics13->DrawFrame(xi, yi, xf, yf);
  fjpsics13->SetXTitle("p_{T}/M");
  fjpsics13->SetYTitle("d#sigma / d#xidy (nb/GeV)");
  jpsics13->Modified();
  jpsics13->SetTitle("");
  
  //plot data
  TGraphAsymmErrors *jsection13 = new TGraphAsymmErrors(njpsics13, datajpsics13[0], datajpsics13[1], datajpsics13[3], datajpsics13[4], datajpsics13[2], datajpsics13[2]);
  jsection13->SetMarkerStyle(20);
  jsection13->SetLineColor(kBlack);
  jsection13->SetMarkerColor(kBlack);
  jsection13->SetMarkerSize(.75);
  jsection13->Draw("P");

  //plot the fitted function
  TF1 *fitjpsics13 = new TF1("jpsi cs fit", "[15]*sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitjpsics13->SetParameter(0, datasigma[9][counter]/(param[34]+mass[3]));
  fitjpsics13->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitjpsics13->SetParameter(i, param[i+7]);
  fitjpsics13->SetParameter(7, param[3]);
  fitjpsics13->SetParameter(8, param[34]+mass[3]);
  for(int i=9; i<13; i++)
    fitjpsics13->SetParameter(i, param[i+5]);
  fitjpsics13->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitjpsics13->SetParameter(14, param[34]+mass[0]);
  fitjpsics13->SetParameter(15, param[31]);
  fitjpsics13->SetLineColor(kBlue);
  fitjpsics13->Draw("lsame");

  //text on the plot
  TLatex ljpsics13;
  ljpsics13.SetTextSize(0.03);
  ljpsics13.DrawLatex(8, 20, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  ljpsics13.DrawLatex(8, 7, Form("P(#chi^{2},ndf) = %.0f%%", 100*chiprob));
  ljpsics13.DrawLatex(0, 2.e-5, "pp 13 TeV");

  //draw legend
  TLegend *legjpsics13 = new TLegend(0.6, 0.75, 0.9, 0.9);
  legjpsics13->SetTextSize(0.03);
  legjpsics13->AddEntry(jsection13, "CMS J/#psi", "p");
  legjpsics13->AddEntry(fitjpsics13, "model", "l");
  legjpsics13->Draw();
  //save jpsi cross section plot
  jpsics13->SaveAs("plots/jpsi_cs_13.pdf");

  //CMS ups(3S)
  const int nups3ccs13=ns[15];
  counter+=njpsics13;
  float dataups3ccs13[5][nups3ccs13];
  for(int i=0; i<nups3ccs13; i++)
    {
      dataups3ccs13[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter], datasigma[8][i+counter], datasigma[9][i+counter]/(param[35]+mass[4]), param[9], param[10], param[11], param[12], param[13], param[4], param[35]+mass[4], param[14], param[15], param[16], param[17], 4);
      dataups3ccs13[1][i]=datasigma[1][i+counter]*param[33]*param[25];
      dataups3ccs13[2][i]=datasigma[2][i+counter]*param[33]*param[25];
      dataups3ccs13[3][i]=dataups3ccs13[0][i]-datasigma[5][i+counter];
      dataups3ccs13[4][i]=datasigma[6][i+counter]-dataups3ccs13[0][i];
    }

  TCanvas *ups3cs13 = new TCanvas("ups313 cross section", "ups313 cross section", 700, 700);
  ups3cs13->SetLogy();
  
  TH1F *fups3cs13 = ups3cs13->DrawFrame(xi, yi, xf, yf);
  fups3cs13->SetXTitle("p_{T}/M");
  fups3cs13->SetYTitle("d#sigma / dp_{T} (nb/GeV)");
  ups3cs13->Modified();
  ups3cs13->SetTitle("");

  //plot CMS data
  TGraphAsymmErrors *u3csection13 = new TGraphAsymmErrors(nups3ccs13, dataups3ccs13[0], dataups3ccs13[1], dataups3ccs13[3], dataups3ccs13[4], dataups3ccs13[2], dataups3ccs13[2]);
  u3csection13->SetMarkerStyle(20);
  u3csection13->SetLineColor(kBlack);
  u3csection13->SetMarkerColor(kBlack);
  u3csection13->SetMarkerSize(.75);
  u3csection13->Draw("P");

  //plot the fitted function and each contribution
  TF1 *fitups3cs13 = new TF1("ups3 cs fit", "[15]*sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitups3cs13->SetParameter(0, datasigma[9][counter]/(param[35]+mass[4]));
  fitups3cs13->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitups3cs13->SetParameter(i, param[i+7]);
  fitups3cs13->SetParameter(7, param[4]);
  fitups3cs13->SetParameter(8, param[35]+mass[4]);
  for(int i=9; i<13; i++)
    fitups3cs13->SetParameter(i, param[i+5]);
  fitups3cs13->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitups3cs13->SetParameter(14, param[34]+mass[0]);
  fitups3cs13->SetParameter(15, param[32]);
  fitups3cs13->SetLineColor(kBlue);
  fitups3cs13->Draw("lsame");
  
  //text on the plot
  TLatex lups3cs13;
  lups3cs13.SetTextSize(0.03);
  lups3cs13.DrawLatex(8, 50, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lups3cs13.DrawLatex(8, 20, Form("P(#chi^{2},ndf) = %.2g%%", 100*chiprob));
  lups3cs13.DrawLatex(0, 2.e-5, "pp 13 TeV");

  //draw legend
  TLegend *legups3cs13 = new TLegend(0.6, 0.7, 0.9, 0.9);
  legups3cs13->SetTextSize(0.03);
  legups3cs13->AddEntry(u3csection13, "CMS #Upsilon(3S)", "p");
  legups3cs13->AddEntry(fitups3cs13, "model", "l");
  legups3cs13->Draw();
  //save ups(3S) cross section plot
  ups3cs13->SaveAs("plots/ups3S_cs_13.pdf");

  //CMS ups(2S)
  const int nups2ccs13=ns[16];
  counter+=nups3ccs13;
  float dataups2ccs13[5][nups2ccs13];
  for(int i=0; i<nups2ccs13; i++)
    {
      dataups2ccs13[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter],datasigma[8][i+counter], datasigma[9][i+counter]/(param[35]+mass[5]), param[9], param[10], param[11], param[12], param[13], param[5], param[35]+mass[5], param[14], param[15], param[16], param[17], 5);
      dataups2ccs13[1][i]=datasigma[1][i+counter]*param[33]*param[26];
      dataups2ccs13[2][i]=datasigma[2][i+counter]*param[33]*param[26];
      dataups2ccs13[3][i]=dataups2ccs13[0][i]-datasigma[5][i+counter];
      dataups2ccs13[4][i]=datasigma[6][i+counter]-dataups2ccs13[0][i];
    }

  TCanvas *ups2cs13 = new TCanvas("ups213 cross section", "ups213 cross section", 700, 700);
  ups2cs13->SetLogy();
  
  TH1F *fups2cs13 = ups2cs13->DrawFrame(xi, yi, xf, yf);
  fups2cs13->SetXTitle("p_{T}/M");
  fups2cs13->SetYTitle("d#sigma / dp_{T} (nb/GeV)");
  ups2cs13->Modified();
  ups2cs13->SetTitle("");

  //plot CMS data
  TGraphAsymmErrors *u2csection13 = new TGraphAsymmErrors(nups2ccs13, dataups2ccs13[0], dataups2ccs13[1], dataups2ccs13[3], dataups2ccs13[4], dataups2ccs13[2], dataups2ccs13[2]);
  u2csection13->SetMarkerStyle(20);
  u2csection13->SetLineColor(kBlack);
  u2csection13->SetMarkerColor(kBlack);
  u2csection13->SetMarkerSize(.75);
  u2csection13->Draw("P");

  //plot the fitted function and each contribution
  TF1 *fitups2cs13 = new TF1("ups2 cs fit", "[15]*sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitups2cs13->SetParameter(0, datasigma[9][counter]/(param[35]+mass[5]));
  fitups2cs13->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitups2cs13->SetParameter(i, param[i+7]);
  fitups2cs13->SetParameter(7, param[5]);
  fitups2cs13->SetParameter(8, param[35]+mass[5]);
  for(int i=9; i<13; i++)
    fitups2cs13->SetParameter(i, param[i+5]);
  fitups2cs13->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitups2cs13->SetParameter(14, param[34]+mass[0]);
  fitups2cs13->SetParameter(15, param[32]);
  fitups2cs13->SetLineColor(kBlue);
  fitups2cs13->Draw("lsame");

  //text on the plot
  TLatex lups2cs13;
  lups2cs13.SetTextSize(0.03);
  lups2cs13.DrawLatex(8, 50, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lups2cs13.DrawLatex(8, 20, Form("P(#chi^{2},ndf) = %.2g%%", 100*chiprob));
  lups2cs13.DrawLatex(0, 2.e-5, "pp 13 TeV");

  //draw legend
  TLegend *legups2cs13 = new TLegend(0.6, 0.7, 0.9, 0.9);
  legups2cs13->SetTextSize(0.03);
  legups2cs13->AddEntry(u2csection13, "CMS #Upsilon(2S)", "p");
  legups2cs13->AddEntry(fitups2cs13, "model", "l");
  legups2cs13->Draw();
  //save ups(2S) cross section plot
  ups2cs13->SaveAs("plots/ups2S_cs_13.pdf");

  //CMS ups(1S)
  const int nups1ccs13=ns[17];
  counter+=nups2ccs13;
  float dataups1ccs13[5][nups1ccs13];
  for(int i=0; i<nups1ccs13; i++)
    {
      dataups1ccs13[0][i]=avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter],datasigma[8][i+counter], datasigma[9][i+counter]/(param[35]+mass[6]), param[9], param[10], param[11], param[12], param[13], param[8], param[35]+mass[6], param[14], param[15], param[16], param[17], 6);
      dataups1ccs13[1][i]=datasigma[1][i+counter]*param[33]*param[27];
      dataups1ccs13[2][i]=datasigma[2][i+counter]*param[33]*param[27];
      dataups1ccs13[3][i]=dataups1ccs13[0][i]-datasigma[5][i+counter];
      dataups1ccs13[4][i]=datasigma[6][i+counter]-dataups1ccs13[0][i];
    }
  
  TCanvas *ups1cs13 = new TCanvas("ups113 cross section", "ups113 cross section", 700, 700);
  ups1cs13->SetLogy();
  
  TH1F *fups1cs13 = ups1cs13->DrawFrame(xi, yi, xf, yf);
  fups1cs13->SetXTitle("p_{T}/M");
  fups1cs13->SetYTitle("d#sigma / dp_{T} (nb/GeV)");
  ups1cs13->Modified();
  ups1cs13->SetTitle("");

  //plot CMS data
  TGraphAsymmErrors *u1csection13 = new TGraphAsymmErrors(nups1ccs13, dataups1ccs13[0], dataups1ccs13[1], dataups1ccs13[3], dataups1ccs13[4], dataups1ccs13[2], dataups1ccs13[2]);
  u1csection13->SetMarkerStyle(20);
  u1csection13->SetLineColor(kBlack);
  u1csection13->SetMarkerColor(kBlack);
  u1csection13->SetMarkerSize(.75);
  u1csection13->Draw("P");
  
  //plot the fitted function and each contribution
  TF1 *fitups1cs13 = new TF1("ups1 cs fit", "[15]*sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, xf);
  fitups1cs13->SetParameter(0, datasigma[9][counter]/(param[35]+mass[6]));
  fitups1cs13->SetParameter(1, datasigma[8][counter]/2);
  for(int i=2; i<7; i++)
    fitups1cs13->SetParameter(i, param[i+7]);
  fitups1cs13->SetParameter(7, param[8]);
  fitups1cs13->SetParameter(8, param[35]+mass[6]);
  for(int i=9; i<13; i++)
    fitups1cs13->SetParameter(i, param[i+5]);
  fitups1cs13->SetParameter(13, datasigma[9][0]/(param[34]+mass[0]));
  fitups1cs13->SetParameter(14, param[34]+mass[0]);
  fitups1cs13->SetParameter(15, param[32]);
  fitups1cs13->SetLineColor(kBlue);
  fitups1cs13->Draw("lsame");

  //text on the plot
  TLatex lups1cs13;
  lups1cs13.SetTextSize(0.03);
  lups1cs13.DrawLatex(8, 50, Form("#chi^{2}/ndf = %.0f/%d", minimum, ndf));
  lups1cs13.DrawLatex(8, 20, Form("P(#chi^{2},ndf) = %.2g%%", 100*chiprob));
  lups1cs13.DrawLatex(0, 2.e-5, "pp 13 TeV");

  //draw legend
  TLegend *legups1cs13 = new TLegend(0.6, 0.7, 0.9, 0.9);
  legups1cs13->SetTextSize(0.03);
  legups1cs13->AddEntry(u1csection13, "CMS #Upsilon(1S)", "p");
  legups1cs13->AddEntry(fitups1cs13, "model", "l");
  legups1cs13->Draw();
  //save ups(1S) cross section plot
  ups1cs13->SaveAs("plots/ups1S_cs_13.pdf");
  
  //part: plotting the "pulls" - maintaining the uncertainty, though
  double cs, y, pt, dpt, dy, npt;
  int nmin = 3, ny=4, nsteps;
  double dn = (8.-3.)/39., signc = sigplot(datasigma[9][0]/(param[34]+mass[0]), 4, 0, param[9], param[10], param[11], param[12], param[13], 1., param[34]+mass[0], param[14], param[15], param[16], param[17]);
  
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[34]+mass[0]), pt, y, param[9], param[10], param[11], param[12], param[13], param[0], param[34]+mass[0], param[14], param[15], param[16], param[17]);
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[34]+mass[0]), pt, y, param[9], param[10], param[11], param[12], param[13], param[0], param[34]+mass[0], param[14], param[15], param[16], param[17]);
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[34]+mass[1]), pt, y, param[9], param[10], param[11], param[12], param[13], param[1], param[34]+mass[1], param[14], param[15], param[16], param[17]);
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[34]+mass[2]), pt, y, param[9], param[10], param[11], param[12], param[13], param[2], param[34]+mass[2], param[14], param[15], param[16], param[17]);
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[34]+mass[3]), pt, y, param[9], param[10], param[11], param[12], param[13], param[3], param[34]+mass[3], param[14], param[15], param[16], param[17]);
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[35]+mass[4]), pt, y, param[9], param[10], param[11], param[12], param[13], param[4], param[35]+mass[4], param[14], param[15], param[16], param[17]);
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[35]+mass[4]), pt, y, param[9], param[10], param[11], param[12], param[13], param[4], param[35]+mass[4], param[14], param[15], param[16], param[17]);
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[35]+mass[5]), pt, y, param[9], param[10], param[11], param[12], param[13], param[5], param[35]+mass[5], param[14], param[15], param[16], param[17]);
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[35]+mass[5]), pt, y, param[9], param[10], param[11], param[12], param[13], param[5], param[35]+mass[5], param[14], param[15], param[16], param[17]);
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[35]+mass[6]), pt, y, param[9], param[10], param[11], param[12], param[13], param[8], param[35]+mass[6], param[14], param[15], param[16], param[17]);
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
	      cs+=sigplot(datasigma[9][i+counter]/(param[35]+mass[6]), pt, y, param[9], param[10], param[11], param[12], param[13], param[8], param[35]+mass[6], param[14], param[15], param[16], param[17]);
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
	      cs+=param[31]*sigplot(datasigma[9][i+counter]/(param[34]+mass[0]), pt, y, param[9], param[10], param[11], param[12], param[13], param[0], param[34]+mass[0], param[14], param[15], param[16], param[17]);
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
	      cs+=param[31]*sigplot(datasigma[9][i+counter]/(param[34]+mass[3]), pt, y, param[9], param[10], param[11], param[12], param[13], param[3], param[34]+mass[3], param[14], param[15], param[16], param[17]);
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
	      cs+=param[32]*sigplot(datasigma[9][i+counter]/(param[35]+mass[4]), pt, y, param[9], param[10], param[11], param[12], param[13], param[4], param[35]+mass[4], param[14], param[15], param[16], param[17]);
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
	      cs+=param[32]*sigplot(datasigma[9][i+counter]/(param[35]+mass[5]), pt, y, param[9], param[10], param[11], param[12], param[13], param[5], param[35]+mass[5], param[14], param[15], param[16], param[17]);
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
	      cs+=param[32]*sigplot(datasigma[9][i+counter]/(param[35]+mass[6]), pt, y, param[9], param[10], param[11], param[12], param[13], param[8], param[35]+mass[6], param[14], param[15], param[16], param[17]);
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
  pups113->SaveAs("plots/ups1S_pull_13.pdf");
}
