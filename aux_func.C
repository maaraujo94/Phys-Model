////////////////////////////////////////////////////////////////////////
// auxiliary file that defines the functions used in main macros      //
////////////////////////////////////////////////////////////////////////

#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TStyle.h"

//////////////////////////////////////////
// part: auxiliary functions for fitting
//////////////////////////////////////////

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
//par-> [0]: sqsfm, [1]: pT/M, [2]: y, [3]: A, [4]: beta, [5]: tau, [6]: rho, [7]: delta, [8]: L, [9]: M, [10]: b, [11]: c, [12]: d, [13]: e
double sig(double *par)
{
  double sum = 0;
  
  double ylim = acosh( 1/par[1] * ( par[0] - sqrt(1+par[1]*par[1]) * cosh(par[2])));
  double final = 26;
  double dy = 2*ylim/(final-1.);
  double y4 = -ylim;
  double pars[13]={par[0], par[1], par[2], y4, par[3], par[4], par[5], par[6], par[7], par[10], par[11], par[12], par[13]};
  
  for(int i=0; i < final; i++)
    {
      pars[3] = y4;
      sum += ( par[1] / pow(par[0], 2) ) * integf(pars) * dy * ( par[8] / pow(par[9],5));
      y4 += dy;
    }
 
  return sum;
}

// general contribution to chisquare for any state
// mqq can be defined as the user prefers: M(state) or fixed Mc, Mb
double statecont(vector <vector <double> > &datasigma, double *pars, double scale, double ptmmin,
		 double mqq, double Lstate, double lumi, double br,
		 int pos_i, int length)
{
  double m_psip = 3.686, dpt, dy, ratio, factor = 0, cs;
  int nmin = 3, nmax = 8, ny = 4, npt, nsteps;
  double dn = (double)(nmax-nmin) / 39.;

  pars[0] = datasigma[9][pos_i]/mqq;
  pars[9] = mqq;
  pars[8] = Lstate;

  //cycle over every position in datasigma for that state
  //for each data pt, averages model over pt/M, y
  for(int i = pos_i; i < pos_i+length; i++)
    if(datasigma[0][i] >= ptmmin) {

      dpt = datasigma[6][i]-datasigma[5][i];
      dy = datasigma[8][i]-datasigma[7][i];
      npt = nmin + int(0.5 + dn*(dpt*m_psip-1.));
      nsteps = npt*ny;

      cs = 0.;
      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pars[1] = datasigma[5][i]+xpt*dpt/(npt-1.)+1e-10;
	  for(int xy = 0; xy < ny; xy++)
	    {
	      pars[2] = datasigma[7][i]+xy*dy/(ny-1.);
	      cs+=sig(pars);
	    }
	}
      cs/=(scale*nsteps);

      //lambda=0.;
      //ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);
      ratio = 1.;

      factor-=0.5*((cs-datasigma[1][i]*ratio*lumi*br)/(datasigma[2][i]*ratio*lumi*br))*((cs-datasigma[1][i]*ratio*lumi*br)/(datasigma[2][i]*ratio*lumi*br)); 
    }
  return factor;
}

// general constribution to chisquare for ratios
// mcb can be defined as the user prefers: M(state) or fixed Mc, Mb
double ratiocont(vector <vector <double> > &datasigma, double *pars, double scale, double ptmmin,
		 double mqq, double Lstate1, double Lstate2, double br,
		 int pos_i, int length)
{
  double m_psip = 3.686, dpt, dy, ratio, factor = 0, cs1, cs2;
  int nmin = 3, nmax = 8, ny = 4, npt, nsteps;
  double dn = (double)(nmax-nmin) / 39.;

  pars[0] = datasigma[9][pos_i]/mqq;
  pars[9] = mqq;

  //cycle over every position in datasigma for that state
  //for each data pt, averages model over pt/M, y
  for(int i = pos_i; i < pos_i+length; i++) {

    dpt = datasigma[6][i]-datasigma[5][i];
    dy = datasigma[8][i]-datasigma[7][i];
    npt = nmin + int(0.5 + dn*(dpt*m_psip-1.));
    nsteps = npt*ny;
    
    cs1 = 0.;
    cs2 = 0.;
    for(int xpt = 0; xpt < npt; xpt++)
      {
	pars[1] = datasigma[5][i]+xpt*dpt/(npt-1.)+1e-10;
	for(int xy = 0; xy < ny; xy++)
	  {
	    pars[2] = datasigma[7][i]+xy*dy/(ny-1.);
	    pars[8] = Lstate1;
	    cs1+=sig(pars);
	    pars[8] = Lstate2;
	    cs2+=sig(pars);
	  }
      }
    cs1/=(scale*nsteps);
    cs2/=(scale*nsteps);
    
    //lambda1=0.;
    //lambda2=0.;
    //ratio=(1+(1-lambda1)/(3+lambda1)*datasigma[4][i])/(1+(1-lambda2)/(3+lambda2)*datasigma[4][i]);
    
    ratio = 1.;
    
    factor-=0.5*((cs1/cs2-datasigma[1][i]*ratio*br)/(datasigma[2][i]*ratio*br))*((cs1/cs2-datasigma[1][i]*ratio*br)/(datasigma[2][i]*ratio*br));
  }
  return factor;
}

//used to obtain estimate for initialization of the L value
double lest(vector <vector <double> > &datasigma, float PTMNORM, int init_i, int len_i)
{
  double min_dpt = 100;
  int min_i = 0;
  
  for(int i = init_i; i < init_i + len_i; i++)
    if(abs(datasigma[0][i] - PTMNORM) < min_dpt) {
      min_i = i;
      min_dpt = abs(datasigma[0][i] - PTMNORM);
    }
  
  return datasigma[1][min_i];
}

///////////////////////////////////////////
// part: auxiliary functions for plotting
///////////////////////////////////////////

//sigma function for plotting
//takes params individually rather than through pointer 
double sigplot(double sqsfm, double pTM, double y, double A, double beta, double tau, double rho, double delta, double L, double M, double b, double c, double d, double e)
{
  double sum = 0;
  
  double ylim = acosh( 1/pTM * ( sqsfm - sqrt(1+pTM*pTM) * cosh(y)));
  double final = 26;
  double dy = 2*ylim/(final-1.);
  double y4 = -ylim;
  double pars[13]={sqsfm, pTM, y, y4, A, beta, tau, rho, delta, b, c, d, e};
  
  for(int i=0; i < final; i++)
    {
      pars[3] = y4;
      sum += ( pTM / pow(sqsfm, 2) ) * integf(pars) * dy * ( L / pow(M,5));
      y4 += dy;
    }
 
  return sum;
}

//calculates the average pT/M of a bin by weighing it with the model cs
//uint and lint are the upper (sigma*pT/M) and lower (sigma) integrals
double avgptm(double lbound, double ubound, double lybound, double uybound, double sqsfm, double A, double rho, double tau, double beta, double delta, double L, double M, double b, double c, double d, double e, int state)
{
  float mass[7]={3.686, 3.556, 3.511, 3.097, 10.355, 10.023, 9.460};
  
  double step=1.;
  double stepy = 0.25;
  double uint=0, lint=0;
  double pars[14] = {sqsfm, 0., 0., A, rho, tau, beta, delta, L, M, b, c, d, e};
  for(double pt=lbound+1e-10; pt<ubound+1.e-5; pt+=step/mass[state])
    for(double y=lybound; y<uybound+1.e-5; y+=stepy)
      {
	pars[1] = pt;
	pars[2] = y;
      	uint+=sig(pars)*pt;
	lint+=sig(pars);
      }
  return uint/lint;
}

// plots data and fit function associated to a specific ratio of states, sqrt(s)
// returns int val at which datasigma readout was left off
// vectors len, lumibr, legtitles must be of size nstates, lumibr = lumi*br
// TODO: Figure out how I'll do the different y bins: all in the same plot, one plot per y bin, only one f for all 5 bins?
int rplot(vector <vector <double> > &datasigma, int init_i, const int nstates, int* len, 
	  double mqq, double *param, int Lpos1, int Lpos2, double *lumibr,
	  double chisquare, int ndf, double chiprob,
	  int *mkrStyle, string* legtitles, string savename) {

  double m_psip = 3.686;

  //defining canvas
  TCanvas *c = new TCanvas("cross section ratio", "cross section ratio", 700, 700);
  
  TH1F *fc = c->DrawFrame(0.01, 0, 9, 1.49);
  fc->SetXTitle("p_{T}/M");
  fc->SetYTitle("ratio");
  fc->GetYaxis()->SetTitleOffset(1);
  c->Modified();
  c->SetTitle("");
  
  //cycle to define TGraphAsymmErrors
  TGraphAsymmErrors** gc = new TGraphAsymmErrors*[nstates];
  int counter = init_i;
 
  for(int ctr = 0; ctr < nstates; ctr++) {
    const int ndata = len[ctr];
    float datapts[5][ndata];

    for(int i = 0; i < ndata; i++) {
      datapts[0][i] = datasigma[0][i+counter];
      datapts[1][i] = datasigma[1][i+counter]*lumibr[ctr];
      datapts[2][i] = datasigma[2][i+counter]*lumibr[ctr];
      datapts[3][i] = datapts[0][i] - datasigma[5][i+counter];
      datapts[4][i] = datasigma[6][i+counter] - datapts[0][i];
    }
    counter += len[ctr];

    gc[ctr] = new TGraphAsymmErrors(ndata, datapts[0], datapts[1], datapts[3], datapts[4], datapts[2], datapts[2]);
    gc[ctr]->SetMarkerStyle(mkrStyle[ctr]);
    gc[ctr]->SetLineColor(kBlack);
    gc[ctr]->SetMarkerColor(kBlack);
    gc[ctr]->SetMarkerSize(.75);
    gc[ctr]->Draw("P");
  }
  
  //plot the fitted function
  TF1 *fit = new TF1("ratio fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [13], [9], [10], [11], [12])/sigplot([0], x, [1], [2], [3], [4], [5], [6], [8], [13], [9], [10], [11], [12])", 0, 10);
  fit->SetParameter(0, datasigma[9][init_i]/mqq);
  fit->SetParameter(1, datasigma[8][init_i]/2);
  for(int i=2; i<7; i++)
    fit->SetParameter(i, param[i+7]);
  fit->SetParameter(7, param[Lpos1]);
  fit->SetParameter(8, param[Lpos2]);
  for(int i=9; i<13; i++)
    fit->SetParameter(i, param[i+5]);
  fit->SetParameter(13, mqq);
  fit->SetLineColor(kBlue);
  fit->Draw("lsame");
  
  //text on the plot
  TLatex lc;
  lc.SetTextSize(0.03);
  lc.DrawLatex(0.5, 1.3, Form("#chi^{2}/ndf = %.0f/%d", chisquare, ndf));
  lc.DrawLatex(0.5, 1.2, Form("P(#chi^{2},ndf) = %.1f%%", 100*chiprob));
  lc.DrawLatex(0.25, 0.05, Form("pp %.0f TeV", datasigma[9][init_i]/1000.));
  
  //draw legend
  double endpt = 0.8 - (double)nstates*0.05;  
  TLegend *leg = new TLegend(0.6, endpt, 0.9, 0.9);
  leg->SetTextSize(0.03);
  for(int i = 0; i < nstates; i++)  {
    const char *st = legtitles[i].c_str();
    leg->AddEntry(gc[i], st, "p");
  }
  leg->AddEntry(fit, "model", "l");
  leg->Draw();
  
  //save plot
  savename = savename+".pdf";
  const char* save = savename.c_str();
  c->SaveAs(save);
  c->Destructor();

  return counter;
}

// plots data and fit function + pulls associated to a specific state, sqrt(s)
// returns int val at which datasigma readout was left off
// vectors len, lumibr, legtitles must be of size nstates, lumibr = lumi*br
// TODO: Figure out how I'll do the different y bins: all in the same plot, one plot per y bin, only one f for all 5 bins?
int csplotpull(vector <vector <double> > &datasigma, int init_i, const int nstates, int* len, 
	   double mqq, double *param, int Lpos, int state, double *lumibr,
	   double ptmmin, double chisquare, int ndf, double chiprob,
	   int *mkrStyle, string* legtitles, string savename) {

  double m_psip = 3.686;
  float xi = -2, xf = 49.9;

  int nmin = 3, nmax = 8, npt, ny = 4, nsteps;
  double cs, pt, y, dpt, dy, dn = (nmax - nmin) / 39.;
  double norm = sigplot(datasigma[9][0]/m_psip, 4, 0, param[9], param[10], param[11], param[12], param[13], 1., m_psip, param[14], param[15], param[16], param[17]);
  
  //defining canvas
  TCanvas *c = new TCanvas("cross section", "cross section", 700, 700);
  c->SetLogy();
  
  TH1F *fc = c->DrawFrame(xi, 1.01e-5, xf, 9.99e2);
  fc->SetXTitle("p_{T}/M");
  fc->SetYTitle("d#sigma / d#xidy (nb/GeV)");
  fc->GetYaxis()->SetTitleOffset(1);
  c->Modified();
  c->SetTitle("");
  
  //cycle to define TGraphAsymmErrors
  TGraphAsymmErrors** gc = new TGraphAsymmErrors*[nstates];
  TGraphAsymmErrors** gp = new TGraphAsymmErrors*[nstates];
  int counter = init_i;
 
  for(int ctr = 0; ctr < nstates; ctr++) {
    const int ndata = len[ctr];
    float datapts[5][ndata];
    float pullpts[2][ndata];
    
    for(int i = 0; i < ndata; i++) {
      //data for cross section plots
      datapts[0][i] = avgptm(datasigma[5][i+counter], datasigma[6][i+counter], datasigma[7][i+counter], datasigma[8][i+counter], datasigma[9][i+counter]/mqq, param[9], param[10], param[11], param[12], param[13], param[Lpos], mqq, param[14], param[15], param[16], param[17], state);
      datapts[1][i] = datasigma[1][i+counter]*lumibr[ctr];
      datapts[2][i] = datasigma[2][i+counter]*lumibr[ctr];
      datapts[3][i] = datapts[0][i] - datasigma[5][i+counter];
      datapts[4][i] = datasigma[6][i+counter] - datapts[0][i];
    
      //calculation of the cross section as an average for higher precision pulls
      cs = 0;
      dpt = datasigma[6][i+counter]-datasigma[5][i+counter];
      dy = datasigma[8][i+counter]-datasigma[7][i+counter];
      npt = nmin + int(0.5 + dn*(dpt*m_psip - 1.));

      for(int xpt = 0; xpt < npt; xpt++)
	{
	  pt = datasigma[5][i+counter]+xpt*dpt/(npt-1.);
	  for(int xy = 0; xy < ny; xy++)
	    {
	      y = datasigma[7][i+counter]+xy*dy/(ny-1.);
	      cs += sigplot(datasigma[9][i+counter]/mqq, pt, y, param[9], param[10], param[11], param[12], param[13], param[Lpos], mqq, param[14], param[15], param[16], param[17]);
	    }
	}
      nsteps = npt*ny;
      cs/=(norm*nsteps);

      pullpts[0][i] = (datapts[1][i] - cs) / cs;
      pullpts[1][i] = datapts[2][i] / cs;
    }
    counter += len[ctr];

    gc[ctr] = new TGraphAsymmErrors(ndata, datapts[0], datapts[1], datapts[3], datapts[4], datapts[2], datapts[2]);
    gc[ctr]->SetMarkerStyle(mkrStyle[ctr]);
    gc[ctr]->SetLineColor(kBlack);
    gc[ctr]->SetMarkerColor(kBlack);
    gc[ctr]->SetMarkerSize(.75);
    gc[ctr]->Draw("P");
    
    gp[ctr] = new TGraphAsymmErrors(ndata, datapts[0], pullpts[0], datapts[3], datapts[4], pullpts[1], pullpts[1]);
    gp[ctr]->SetMarkerStyle(mkrStyle[ctr]);
    gp[ctr]->SetLineColor(kBlack);
    gp[ctr]->SetMarkerColor(kBlack);
    gp[ctr]->SetMarkerSize(.75); 
    for(int i = 0; i < ndata; i++)
      if(datapts[0][i] < ptmmin)
	gp[ctr]->RemovePoint(0);
  }
  
  //plot the fitted function
  TF1 *fit = new TF1("cs fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12])/sigplot([13], 4., 0., [2], [3], [4], [5], [6], 1., [14], [9], [10], [11], [12])", ptmmin, 49.9);
  fit->SetParameter(0, datasigma[9][init_i]/mqq);
  fit->SetParameter(1, datasigma[8][init_i]/2);
  for(int i=2; i<7; i++)
    fit->SetParameter(i, param[i+7]);
  fit->SetParameter(7, param[Lpos]);
  fit->SetParameter(8, mqq);
  for(int i=9; i<13; i++)
    fit->SetParameter(i, param[i+5]);
  fit->SetParameter(13, datasigma[9][0]/(m_psip));
  fit->SetParameter(14, m_psip);
  fit->SetLineColor(kBlue);
  fit->Draw("lsame");
  
  //text on the plot
  TLatex lc;
  lc.SetTextSize(0.03);
  lc.DrawLatex(8, 20, Form("#chi^{2}/ndf = %.0f/%d", chisquare, ndf));
  lc.DrawLatex(8, 7, Form("P(#chi^{2},ndf) = %.1f%%", 100*chiprob));
  lc.DrawLatex(0, 2.e-5, Form("pp %.0f TeV", datasigma[9][init_i]/1000.));

  //draw legend
  double endpt = 0.8 - (double)nstates*0.05;  
  TLegend *leg = new TLegend(0.6, endpt, 0.9, 0.9);
  leg->SetTextSize(0.03);
  for(int i = 0; i < nstates; i++)  {
    const char *st = legtitles[i].c_str();
    leg->AddEntry(gc[i], st, "p");
  }
  leg->AddEntry(fit, "model", "l");
  leg->Draw();
  
  //save cross section plot
  string csname = savename+"_cs.pdf";
  const char* save = csname.c_str();
  c->SaveAs(save);

  //repeat plotting but for pulls
  c->Clear();
  c->SetLogy(0);

  TH1F *fp = c->DrawFrame(xi, -2, xf, 2);
  fp->SetXTitle("p_{T}/M");
  fp->SetYTitle("pulls");
  fp->GetYaxis()->SetTitleOffset(1);
  c->Modified();
  c->SetTitle("");
   
  //plot line at zero
  TF1 *zero = new TF1("zero", "0", xi, xf);
  zero->SetLineColor(kBlue);
  zero->SetLineStyle(7);
  zero->Draw("lsame");

  //draw pulls
  for(int i = 0; i < nstates; i++)
    gp[i]->Draw("P");
  
  //text on the plot
  TLatex lp;
  lp.SetTextSize(0.03);
  lp.DrawLatex(6, 1.7, Form("#chi^{2}/ndf = %.0f/%d", chisquare, ndf));
  lp.DrawLatex(6, 1.3, Form("P(#chi^{2},ndf) = %.1f%%", 100*chiprob));
  lp.DrawLatex(0, -1.7, Form("pp %.0f TeV", datasigma[9][init_i]/1000.));

  //draw legend
  TLegend *legp = new TLegend(0.6, endpt+0.05, 0.9, 0.9);
  legp->SetTextSize(0.03);
  for(int i = 0; i < nstates; i++)  {
    const char *st = legtitles[i].c_str();
    legp->AddEntry(gc[i], st, "p");
  }
  legp->Draw();
  
  //save pulls plot
  string pullname = savename+"_cs_pull.pdf";
  const char* savep = pullname.c_str();
  c->SaveAs(savep);
  c->Destructor();

  return counter;
}
