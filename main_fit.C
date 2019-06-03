#include "TFitter.h"
#include "TMath.h"
#include "TROOT.h"
#include <fstream>
#include <iostream>
using namespace std;

#include "aux_func.C"
#include "aux_input.C"

vector <vector <double> > datasigma(10);
//psi', chi_c2, chi_c1, jpsi, ups(3S), ups(2S), ups(1S)
float mass[7]={3.686, 3.556, 3.511, 3.097, 10.355, 10.023, 9.460};

int ONEC = 1;

//////////////////////////////////////
//part: the log-likelihood function
/////////////////////////////////////

//defining the log-likelihood
double myFunction(double Lpsip, double Lchic2, double Lchic1, double Ljpsi, double Lups3, double Lups2, double Lchib2, double Lchib1, double Lups1, double A, double beta, double tau, double rho, double delta, double b, double c, double d, double e, double brpsipdm, double brpsipdp, double brc2jpsi, double brc1jpsi, double brjpsidm, double brb2ups1, double brb1ups1, double brups3dm, double brups2dm, double brups1dm, double lcms, double latlas, double latlasups, double lrcc, double lrbb, double lcms13, double cc, double cb, double ptmmin)
{
  Lchib2 = Lchic2;
  Lchib1 = Lchic1;
  double factor, cs, cs_aux, lambda, lambda_aux, ratio;
  double loglike = 0.;
  if(ONEC == 1)
    cc = cb;
  double Mb, Mc = cc + mass[0];

  double xinorm = 4, ynorm = 0;
  int ny = 4, npt, nmin = 3, nmax = 8;
  double dn = (double)(nmax-nmin) / 39.;
  double dpt, dy;
  int nsteps;

  //defining the arguments of sig(). Must redefine pars[0] when fitting data at different sqrt(s) (13 TeV, etc)
  double pars[14] = {7000./Mc, xinorm, ynorm, A, beta, tau, rho, delta, 1., Mc, b, c, d, e};
  double signc = sig(pars);
  double signb = sig(pars);

  //cross section data chisquares
  for(int i=0; i<(int)datasigma[0].size(); i++)
    {
      //7 TeV data
      //CMS psiprime cross-section
      if(datasigma[3][i]==0 && datasigma[0][i]>=ptmmin)
	{
	  cs=0.;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
	  Mc = cc + mass[0];
	  pars[0] = datasigma[9][i]/Mc;
	  pars[9] = Mc;
	  pars[8] = Lpsip;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1.);
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1.);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signc*nsteps);
	  
	  lambda=0.;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);
	  factor=-0.5*((cs-datasigma[1][i]*ratio*lcms*brpsipdm)/(datasigma[2][i]*ratio*lcms*brpsipdm))*((cs-datasigma[1][i]*ratio*lcms*brpsipdm)/(datasigma[2][i]*ratio*lcms*brpsipdm));
	  loglike+=factor;
	}

      //ATLAS psiprime cross-section
      else if(datasigma[3][i]==1 && datasigma[0][i]>=ptmmin)
	{
	  cs=0.;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
	  Mc = cc + mass[0];
	  pars[0] = datasigma[9][i]/Mc;
	  pars[9] = Mc;
	  pars[8] = Lpsip;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1.);
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1.);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signc*nsteps);
	  
	  lambda=0.;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);
	  factor=-0.5*((cs-datasigma[1][i]*ratio*latlas*brpsipdp*brjpsidm)/(datasigma[2][i]*ratio*latlas*brpsipdp*brjpsidm))*((cs-datasigma[1][i]*ratio*latlas*brpsipdp*brjpsidm)/(datasigma[2][i]*ratio*latlas*brpsipdp*brjpsidm));
	  loglike+=factor;
	}

      //ATLAS chic2 cross-section
      else if(datasigma[3][i]==2 && datasigma[0][i]>=ptmmin)
	{
	  cs=0.;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
	  Mc = cc + mass[1];
	  pars[0] = datasigma[9][i]/Mc;
	  pars[8] = Lchic2;
	  pars[9] = Mc;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1.);
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1.);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signc*nsteps);
	  
	  lambda=0.;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);
	  factor=-0.5*((cs-datasigma[1][i]*ratio*latlas*brc2jpsi*brjpsidm)/(datasigma[2][i]*ratio*latlas*brc2jpsi*brjpsidm))*((cs-datasigma[1][i]*ratio*latlas*brc2jpsi*brjpsidm)/(datasigma[2][i]*ratio*latlas*brc2jpsi*brjpsidm));
	  loglike+=factor;
	}

      //ATLAS chic1 cross-section
      else if(datasigma[3][i]==3 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mc = cc + mass[2];
	  pars[0] = datasigma[9][i]/Mc;
	  pars[8] = Lchic1;
	  pars[9] = Mc;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1);
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signc*nsteps);
	  	      
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);
	  factor=-0.5*((cs-datasigma[1][i]*ratio*latlas*brc1jpsi*brjpsidm)/(datasigma[2][i]*ratio*latlas*brc1jpsi*brjpsidm))*((cs-datasigma[1][i]*ratio*latlas*brc1jpsi*brjpsidm)/(datasigma[2][i]*ratio*latlas*brc1jpsi*brjpsidm));
	  loglike+=factor;
	}

      //CMS j/psi cross-section
      else if(datasigma[3][i]==4 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mc = cc + mass[3];
	  pars[0] = datasigma[9][i]/Mc;
	  pars[8] = Ljpsi;
	  pars[9] = Mc;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1);
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1.);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signc*nsteps);
	  
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*lcms*brjpsidm)/(datasigma[2][i]*ratio*lcms*brjpsidm))*((cs-datasigma[1][i]*ratio*lcms*brjpsidm)/(datasigma[2][i]*ratio*lcms*brjpsidm));
	  loglike+=factor;
	} 
      
      //CMS chic2/chic1 ratio, using all pT/M
      else if(datasigma[3][i]==5)
	{
	  //cs -> chic2; cs_aux -> chic1
	  cs=0;
	  cs_aux=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mc = cc + mass[3];
	  pars[0] = datasigma[9][i]/Mc;
	  pars[9] = Mc;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1);
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  pars[8] = Lchic2;
		  cs+=sig(pars);
		  pars[8] = Lchic1;
		  cs_aux+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signc*nsteps);
	  cs_aux/=(signc*nsteps);
	  
	  lambda=0;
	  lambda_aux=0;
	  ratio=(1+(1-lambda)/(3+lambda)*datasigma[4][i])/(1+(1-lambda_aux)/(3+lambda_aux)*datasigma[4][i]);
	  factor=-0.5*((cs/cs_aux-datasigma[1][i]*ratio*(brc2jpsi/brc1jpsi))/(datasigma[2][i]*ratio*(brc2jpsi/brc1jpsi)))*((cs/cs_aux-datasigma[1][i]*ratio*(brc2jpsi/brc1jpsi))/(datasigma[2][i]*ratio*(brc2jpsi/brc1jpsi)));
	  loglike+=factor;
	}
      
      //CMS chib2/chib1 ratio, using all pT/M
      else if(datasigma[3][i]==6)
	{
	  cs=0;
	  cs_aux=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mb = cb + mass[6];
	  pars[0] = datasigma[9][i]/Mb;
	  pars[9] = Mb;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1);
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  pars[8] = Lchib2;
		  cs+=sig(pars);
		  pars[8] = Lchib1;
		  cs_aux+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signb*nsteps);
	  cs_aux/=(signb*nsteps);
	  
	  lambda=0;
	  lambda_aux=0;
	  ratio=(1+(1-lambda)/(3+lambda)*datasigma[4][i])/(1+(1-lambda_aux)/(3+lambda_aux)*datasigma[4][i]);
	  factor=-0.5*((cs/cs_aux-datasigma[1][i]*ratio*(brb2ups1/brb1ups1))/(datasigma[2][i]*ratio*(brb2ups1/brb1ups1)))*((cs/cs_aux-datasigma[1][i]*ratio*(brb2ups1/brb1ups1))/(datasigma[2][i]*ratio*(brb2ups1/brb1ups1)));
	  loglike+=factor;
	}

      //CMS ups3 cross-section
      else if(datasigma[3][i]==7 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mb = cb + mass[4];
	  pars[0] = datasigma[9][i]/Mb;
	  pars[9] = Mb;
	  pars[8] = Lups3;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1)+1e-10;
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signb*nsteps);
	  
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*lcms*brups3dm)/(datasigma[2][i]*ratio*lcms*brups3dm))*((cs-datasigma[1][i]*ratio*lcms*brups3dm)/(datasigma[2][i]*ratio*lcms*brups3dm));
	  loglike+=factor;
	}

      //ATLAS ups3 cross-section
      else if(datasigma[3][i]==8 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mb = cb + mass[4];
	  pars[0] = datasigma[9][i]/Mb;
	  pars[9] = Mb;
	  pars[8] = Lups3;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1)+1e-10;
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signb*nsteps);
	  
	  lambda=0;ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*latlasups*brups3dm)/(datasigma[2][i]*ratio*latlasups*brups3dm))*((cs-datasigma[1][i]*ratio*latlasups*brups3dm)/(datasigma[2][i]*ratio*latlasups*brups3dm));
	  loglike+=factor;
	}

      //CMS ups2 cross-section
      else if(datasigma[3][i]==9 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mb = cb + mass[5];
	  pars[0] = datasigma[9][i]/Mb;
	  pars[9] = Mb;
	  pars[8] = Lups2;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1)+1e-10;
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signb*nsteps);
	  
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*lcms*brups2dm)/(datasigma[2][i]*ratio*lcms*brups2dm))*((cs-datasigma[1][i]*ratio*lcms*brups2dm)/(datasigma[2][i]*ratio*lcms*brups2dm));
	  loglike+=factor;
	}

      //ATLAS ups2 cross-section
      else if(datasigma[3][i]==10 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mb = cb + mass[5];
	  pars[0] = datasigma[9][i]/Mb;
	  pars[9] = Mb;
	  pars[8] = Lups2;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1)+1e-10;
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signb*nsteps);
	  
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*latlasups*brups2dm)/(datasigma[2][i]*ratio*latlasups*brups2dm))*((cs-datasigma[1][i]*ratio*latlasups*brups2dm)/(datasigma[2][i]*ratio*latlasups*brups2dm));
	  loglike+=factor;
	}

      //CMS ups1 cross-section
      else if(datasigma[3][i]==11 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mb = cb + mass[6];
	  pars[0] = datasigma[9][i]/Mb;
	  pars[9] = Mb;
	  pars[8] = Lups1;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1)+1e-10;
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signb*nsteps);
	  
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*lcms*brups1dm)/(datasigma[2][i]*ratio*lcms*brups1dm))*((cs-datasigma[1][i]*ratio*lcms*brups1dm)/(datasigma[2][i]*ratio*lcms*brups1dm));
	  loglike+=factor;
	}

      //ATLAS ups1 cross-section
      else if(datasigma[3][i]==12 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mb = cb + mass[6];
	  pars[0] = datasigma[9][i]/Mb;
	  pars[9] = Mb;
	  pars[8] = Lups1;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1)+1e-10;
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  cs+=sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signb*nsteps);
	  
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*latlasups*brups1dm)/(datasigma[2][i]*ratio*latlasups*brups1dm))*((cs-datasigma[1][i]*ratio*latlasups*brups1dm)/(datasigma[2][i]*ratio*latlasups*brups1dm));
	  loglike+=factor;
	}
      
      //13 TeV data
      //CMS psiprime cross-section
      if(datasigma[3][i]==13 && datasigma[0][i]>=ptmmin)
	{
	  cs=0.;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1.));
	  Mc = cc + mass[0];
	  pars[0] = datasigma[9][i]/Mc;
	  pars[9] = Mc;
	  pars[8] = Lpsip;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1.);
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1.);
		  cs+=lrcc*sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signc*nsteps);
	  
	  lambda=0.;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);
	  factor=-0.5*((cs-datasigma[1][i]*ratio*lcms13*brpsipdm)/(datasigma[2][i]*ratio*lcms13*brpsipdm))*((cs-datasigma[1][i]*ratio*lcms13*brpsipdm)/(datasigma[2][i]*ratio*lcms13*brpsipdm));
	  loglike+=factor;
	}

      //CMS j/psi cross-section
      else if(datasigma[3][i]==14 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mc = cc + mass[3];
	  pars[0] = datasigma[9][i]/Mc;
	  pars[8] = Ljpsi;
	  pars[9] = Mc;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1);
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1.);
		  cs+=lrcc*sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signc*nsteps);
	  
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*lcms13*brjpsidm)/(datasigma[2][i]*ratio*lcms13*brjpsidm))*((cs-datasigma[1][i]*ratio*lcms13*brjpsidm)/(datasigma[2][i]*ratio*lcms13*brjpsidm));
	  loglike+=factor;
	} 

      //CMS ups3 cross-section
      else if(datasigma[3][i]==15 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mb = cb + mass[4];
	  pars[0] = datasigma[9][i]/Mb;
	  pars[9] = Mb;
	  pars[8] = Lups3;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1)+1e-10;
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  cs+=lrbb*sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signb*nsteps);
	  
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*lcms13*brups3dm)/(datasigma[2][i]*ratio*lcms13*brups3dm))*((cs-datasigma[1][i]*ratio*lcms13*brups3dm)/(datasigma[2][i]*ratio*lcms13*brups3dm));
	  loglike+=factor;
	}

      //CMS ups2 cross-section
      else if(datasigma[3][i]==16 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mb = cb + mass[5];
	  pars[0] = datasigma[9][i]/Mb;
	  pars[9] = Mb;
	  pars[8] = Lups2;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1)+1e-10;
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  cs+=lrbb*sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signb*nsteps);
	  
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*lcms13*brups2dm)/(datasigma[2][i]*ratio*lcms13*brups2dm))*((cs-datasigma[1][i]*ratio*lcms13*brups2dm)/(datasigma[2][i]*ratio*lcms13*brups2dm));
	  loglike+=factor;
	}

      //CMS ups1 cross-section
      else if(datasigma[3][i]==17 && datasigma[0][i]>=ptmmin)
	{
	  cs=0;
	  dpt = datasigma[6][i]-datasigma[5][i];
	  dy = datasigma[8][i]-datasigma[7][i];
	  npt = nmin + int(0.5 + dn*(dpt*mass[0]-1));
	  Mb = cb + mass[6];
	  pars[0] = datasigma[9][i]/Mb;
	  pars[9] = Mb;
	  pars[8] = Lups1;
	  for(int xpt = 0; xpt < npt; xpt++)
	    {
	      pars[1] = datasigma[5][i]+xpt*dpt/(npt-1)+1e-10;
	      for(int xy = 0; xy < ny; xy++)
		{
		  pars[2] = datasigma[7][i]+xy*dy/(ny-1);
		  cs+=lrbb*sig(pars);
		}
	    }
	  nsteps = npt*ny;
	  cs/=(signb*nsteps);
	  
	  lambda=0;
	  ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);	  
	  factor=-0.5*((cs-datasigma[1][i]*ratio*lcms13*brups1dm)/(datasigma[2][i]*ratio*lcms13*brups1dm))*((cs-datasigma[1][i]*ratio*lcms13*brups1dm)/(datasigma[2][i]*ratio*lcms13*brups1dm));
	  loglike+=factor;
	}
    }
  
  //constraints for the BR nuisance parameters
  loglike-=0.5*((brpsipdm-1)/(0.9/7.9))*((brpsipdm-1)/(0.9/7.9));
  loglike-=0.5*((brpsipdp-1)/(0.3/34.46))*((brpsipdp-1)/(0.3/34.46));
  loglike-=0.5*((brc2jpsi-1)/(0.7/19.2))*((brc2jpsi-1)/(0.7/19.2));
  loglike-=0.5*((brc1jpsi-1)/(1.2/33.9))*((brc1jpsi-1)/(1.2/33.9));
  loglike-=0.5*((brjpsidm-1)/(0.033/5.961))*((brjpsidm-1)/(0.033/5.961));
  loglike-=0.5*((brb2ups1-1)/(1.2/19.1))*((brb2ups1-1)/(1.2/19.1));
  loglike-=0.5*((brb1ups1-1)/(2.2/33.9))*((brb1ups1-1)/(2.2/33.9));
  loglike-=0.5*((brups3dm-1)/(0.21/2.18))*((brups3dm-1)/(0.21/2.18));
  loglike-=0.5*((brups2dm-1)/(0.17/1.93))*((brups2dm-1)/(0.17/1.93));
  loglike-=0.5*((brups1dm-1)/(0.05/2.48))*((brups1dm-1)/(0.05/2.48));
  
  //constraints for the luminosity nuisance parameters
  loglike-=0.5*((lcms-1)/2.2e-2)*((lcms-1)/2.2e-2);
  loglike-=0.5*((latlas-1)/1.8e-2)*((latlas-1)/1.8e-2);
  loglike-=0.5*((latlasups-1)/3.9e-2)*((latlasups-1)/3.9e-2);
  loglike-=0.5*((lcms13-1)/2.2e-2)*((lcms13-1)/2.2e-2);
  
  //I am returning the chisquare because I can't figure out how to change the fit method (must check this)
  return -2*loglike;
}

//Minuit-specific "wrapping"
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg)
{
  result=myFunction(par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9], par[10], par[11], par[12], par[13], par[14], par[15], par[16], par[17], par[18], par[19], par[20], par[21], par[22], par[23], par[24], par[25], par[26], par[27], par[28], par[29], par[30], par[31], par[32], par[33], par[34], par[35], par[36]);
}

//main function: collects data, does the fit, plots the results
void fit()
{
  //define variables
  
  ofstream outf;
  const int nstates = 18;
  int counter = 0, ndf, ns[nstates] = {0};
  double param[36], ptmmin = 2., minimum, chinorm, chiprob;
  
  int pos_psip=0, pos_chic2=0, pos_chic1=0, pos_jpsi=0, pos_ups3=0, pos_ups2=0, pos_ups1=0;
  double pt_psip=100, pt_chic2=100, pt_chic1=100, pt_jpsi=100, pt_ups3=100, pt_ups2=100, pt_ups1=100;

  //////////////////////////////////
  //part: getting data from files
  /////////////////////////////////

  input(datasigma, ns);
  
  ////////////////////////////
  //part: fitting the data
  ///////////////////////////

  //define function to minimize, set initial parameters
  TFitter* fit = new TFitter(37);
  fit->SetFCN(minuitFunction);

  float PTMNORM = 4.;
  //calculate estimates for L_QQ
  for(int i=0; i<(int)datasigma[0].size(); i++)
    {
      if((datasigma[3][i]==0 || datasigma[3][i]==1) && abs(datasigma[0][i]-PTMNORM)<pt_psip)
	{
	  pos_psip=i;
	  pt_psip=abs(datasigma[0][i]-PTMNORM);
	}
      else if(datasigma[3][i]==2 && abs(datasigma[0][i]-PTMNORM)<pt_chic2)
	{
	  pos_chic2=i;
	  pt_chic2=abs(datasigma[0][i]-PTMNORM);
	}
      else if(datasigma[3][i]==3 && abs(datasigma[0][i]-PTMNORM)<pt_chic1)
	{
	  pos_chic1=i;
	  pt_chic1=abs(datasigma[0][i]-PTMNORM);
	}
      else if(datasigma[3][i]==4 && abs(datasigma[0][i]-PTMNORM)<pt_jpsi)
	{
	  pos_jpsi=i;
	  pt_jpsi=abs(datasigma[0][i]-PTMNORM);
	}
      else if((datasigma[3][i]==7 || datasigma[3][i]==8) && abs(datasigma[0][i]-PTMNORM)<pt_ups3)
	{
	  pos_ups3=i;
	  pt_ups3=abs(datasigma[0][i]-PTMNORM);
	}
      else if((datasigma[3][i]==9 || datasigma[3][i]==10) && abs(datasigma[0][i]-PTMNORM)<pt_ups2)
	{
	  pos_ups2=i;
	  pt_ups2=abs(datasigma[0][i]-PTMNORM);
	}
      else if((datasigma[3][i]==11 || datasigma[3][i]==12) && abs(datasigma[0][i]-PTMNORM)<pt_ups1)
	{
	  pos_ups1=i;
	  pt_ups1=abs(datasigma[0][i]-PTMNORM);
	}
    }
  pt_psip=datasigma[1][pos_psip];
  pt_chic2=datasigma[1][pos_chic2];
  pt_chic1=datasigma[1][pos_chic1];
  pt_jpsi=datasigma[1][pos_jpsi];
  pt_ups3=datasigma[1][pos_ups3]*200;
  pt_ups2=datasigma[1][pos_ups2]*200;
  pt_ups1=datasigma[1][pos_ups1]*200;

  //cout << myFunction(pt_psip, pt_chic2, pt_chic1, pt_jpsi, pt_ups3, pt_ups2, pt_chic2, pt_chic1, pt_ups1, 0, 1.3, 3.0, 2.3, 0, 1.9, 6, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2) << endl;
  
  fit->SetParameter(0, "L_psip", pt_psip, pt_psip/10, 0, 0);
  fit->SetParameter(1, "L_chic2", pt_chic2, pt_chic2/10, 0, 0);
  fit->SetParameter(2, "L_chic1", pt_chic1, pt_chic1/10, 0, 0);
  fit->SetParameter(3, "L_jpsi", pt_jpsi, pt_jpsi/10, 0, 0);
  //for(int i=1; i<3; i++)
  //  fit->FixParameter(i);
  fit->SetParameter(4, "L_ups3", pt_ups3, pt_ups3/10, 0, 0);
  fit->SetParameter(5, "L_ups2", pt_ups2, pt_ups2/10, 0, 0);
  fit->SetParameter(6, "L_chib2", pt_ups2, pt_ups2/10, 0, 0);
  fit->FixParameter(6);
  fit->SetParameter(7, "L_chib1", pt_ups2, pt_ups2/10, 0, 0);
  fit->FixParameter(7);
  fit->SetParameter(8, "L_ups1", pt_ups1, pt_ups1/10, 0, 0);

  fit->SetParameter(9, "A", 0., 0.1, 0, 0);
  fit->FixParameter(9);
  fit->SetParameter(10, "beta", 1.5, 0.1, 0, 0);
  //fit->FixParameter(10);
  fit->SetParameter(11, "tau", 1., 0.1, 0, 0);
  //fit->FixParameter(11);
  fit->SetParameter(12, "rho", 2.4, 0.1, 0, 0);
  fit->FixParameter(12);
  fit->SetParameter(13, "delta", 0.1, 0.1, 0, 0);
  fit->FixParameter(13);
  fit->SetParameter(14, "b", 1.8, 0.1, 0, 0);
  //fit->FixParameter(14);
  fit->SetParameter(15, "c", 6, 0.1, 0, 0);
  //fit->FixParameter(15);
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
  //for(int i=0; i<2; i++)
  //  {
  //   fit->FixParameter(i+20);
  //  fit->FixParameter(i+23);
  //  fit->FixParameter(i+29);
  //}
  //fit->FixParameter(19);

  fit->SetParameter(31, "Lr_cc", 1, 0.1, 0, 0);
  fit->SetParameter(32, "Lr_bb", 1, 0.1, 0, 0);
  fit->FixParameter(31);
  fit->FixParameter(32);
  fit->SetParameter(33, "L_CMS_13", 1, 0.1, 0, 0);

  fit->SetParameter(34, "Cc", 0, 0.1, 0, 0);
  fit->SetParameter(35, "Cb", 0, 0.1, 0, 0);
  if(ONEC == 1)
    fit->FixParameter(34);
  fit->FixParameter(35);
  fit->SetParameter(36, "ptm_min", ptmmin, 1, 0, 0);
  fit->FixParameter(36);

  fit->ExecuteCommand("MIGRAD",0,0);
  //fit->ExecuteCommand("MIGRAD",0,0);
  //fit->ExecuteCommand("MIGRAD",0,0);

  //register the fit parameters to an array
  for(int i=0; i<36; i++)
    param[i]=fit->GetParameter(i);
  param[6] = param[1];
  param[7] = param[2];
  if(ONEC == 1) param[34] = param[35];
  
  minimum=myFunction(param[0], param[1], param[2], param[3], param[4], param[5], param[6], param[7], param[8], param[9], param[10], param[11], param[12], param[13], param[14], param[15], param[16], param[17], param[18], param[19], param[20], param[21], param[22], param[23], param[24], param[25], param[26], param[27], param[28], param[29], param[30], param[31], param[32], param[33], param[34], param[35], ptmmin);
  counter=0;
  //calculating the number of data points used in the fit
  for(int i=0; i<(int)datasigma[0].size(); i++)
    {
      if(datasigma[0][i]>=ptmmin)
	counter+=1;
      else if(datasigma[0][i]<ptmmin && (datasigma[3][i]==5 || datasigma[3][i]==6))
      counter+=1;
    }
  cout << counter << " data points" << endl;
  //I am adding one degree of freedom for each of the nuisance parameters, this factor will have to be changed if more are added or removed
  counter+=14;
  //counter+=6;
  //the final ndf is obtained by subtracting the number of free parameters
  ndf=counter-fit->GetNumberFreeParameters();
  cout << "ndf " << ndf << endl;
  chinorm=minimum/ndf;
  chiprob=TMath::Prob(minimum,ndf);
  counter=0;
       
  cout << "chi^2 norm : " << chinorm  << endl;
  cout << "chi^2 prob : " << chiprob  << endl;

  outf.open("fit.txt");
  for(int i=0; i<36; i++)
    outf << param[i] << " ";
  for(int i=0; i<36; i++)
    outf << fit->GetParError(i) << " ";
  outf << ptmmin << " ";
  outf << minimum << " " << ndf << endl;
  outf.close();
}
