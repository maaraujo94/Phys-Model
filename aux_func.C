/////////////////////////////////////////////////////////////////////////////
//auxiliary file that defines the functions used to calculate other functions
/////////////////////////////////////////////////////////////////////////////

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

//sigma function for plotting
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

//calculates the average pT/M of a bin by weighing it with the model cross section
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

// general function for any state
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

// general function for ratios
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
