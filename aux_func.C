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
