//early declaraction of minuitFunction for the doFit call in the main class
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg);

#import "aux_func.C"

// class to use for fitting and plotting
class FitModel {
  // maybe variables private once I'm done with the code I guess
public:
  // states: number of states in the fit
  // chisquare, ndf: results of the fit
  // nparam: number of params of each type (norm, shape, nuis)
  // len: length of sigpar vector
  // ptmmin: minimum pT/M to be included in the fit
  // PTMNORM: pT/M for normalization
  int states, ndf, len;
  int nparam[3] = {0};
  double ptmmin = 2., PTMNORM = 4., chisquare;
  string fitresname = "fit.txt";
  
  // auxiliary params for subdivisions of each y, pT bin
  int nmin = 3, nmax = 8, ny = 4, npt, nsteps;
  double dpt, dy, dn = (double)(nmax-nmin) / 39.;

  // maps between general state names and their masses and latex-formatted names
  map<string, double> QMass;
  map<string, string> QName;

  // ns for nr of pts in each state
  // method for the way to record cross-section uncertainty
  // read to determine whether to read datafile
  // *_att for the norms and nuis of each state
  vector <int> ns, method, read;
  vector <vector <int> > norm_att, nuis_att;
  
  // (e)fitpar for the results of the fit, sigpar for the pars in sig(pars)
  // signorm for the normalization applied on each dataset (y, BR effects)
  // datasigma for storage of all (numerical) info of each state
  vector <double> fitpar, efitpar, signorm, sigpar;
  vector <vector <double> > datasigma;

  // statename saves input files
  // auxnames is for strings needed for legends and savenames
  vector <string> statename;
  vector <vector <string > > auxnames;

  // vectors of each class to define each type of fit param
  vector <NormPar> norm;
  vector <ShapePar> shape;
  vector <NuisPar> nuis;

  // initializes the mass <-> state map
  void init() {
    const int nStates = 17;
  
    string StateNm[nStates] = { "jpsi", "chic0", "chic1", "chic2", "psip", "ups1", "chib0_1P", "chib1_1P", "chib2_1P", "ups2", "chib0_2P", "chib1_2P", "chib2_2P", "ups3", "chib0_3P", "chib1_3P", "chib2_3P" };
    double M[nStates] = { 3.097, 3.415, 3.511, 3.556, 3.686, 9.460, 9.859,  9.893,    9.912,   10.023,  10.233,  10.255,  10.269,  10.355, 10.497,  10.512,    10.523 };
    string LatexNm[nStates] = { "J/#psi", "#chi_{c0}", "#chi_{c1}", "#chi_{c2}", "#psi(2S)", "#Upsilon(1S)", "#chi_{b0}(1P)", "#chi_{b1}(1P)", "#chi_{b2}(1P)", "#Upsilon(2S)", "#chi_{b0}(2P)", "#chi_{b1}(2P)", "#chi_{b2}(2P)", "#Upsilon(3S)", "#chi_{b0}(3P)", "#chi_{b1}(3P)", "#chi_{b2}(3P)" }; 
  
    for(int i = 0; i < nStates; i++) {
      QMass[StateNm[i]] = M[i];
      QName[StateNm[i]] = LatexNm[i];
    }
  }
  
  // method to store all filenames
  void getDataList(const char* input_list) {
    ifstream fin;
    string data;
    int iaux = 0;
    vector <string> aux;
    
    // open file and get number of lines
    fin.open(input_list);
    while(getline(fin, data)) iaux++;

    // set size for statename and ns, fix number of states
    statename.resize(iaux, "");
    ns.resize(iaux, 0);
    read.resize(iaux, 0);
    method.resize(iaux, 0);
    signorm.resize(iaux, 1);
    states = iaux;
    
    // back to the start of the file
    fin.clear();
    fin.seekg(0);
    // for all lines, get stateid and file address
    // read parameter determines whether line is later read
    while(getline(fin, data)) {
      aux = parseString(data, "#");
      aux = parseString(aux[aux.size()-1], " ");
      iaux = stoi(aux[0]);    

      if(data.find("#") == string::npos) read[iaux] = 1;
	
      statename[iaux] = aux[1];
      method[iaux] = stoi(aux[2]);
      aux = parseString(aux[3], "*");
      for(int i = 0; i < aux.size(); i++)
	signorm[iaux] *= stof(aux[i]);
    }
    fin.close();
  }
  
  //method to read data from each filename
  void getData() {
    string data;
    ifstream file;
    double nums[12];
    double aux, mass, sqrts;
    vector <string> saux;
    // can't do datasigma size at initialization bc class reads it as a method
    datasigma.resize(11);

    // vector of strings for legends and savenames
    auxnames.resize(3);
    
    // read and input all datafiles from statename
    for(int id = 0; id < states; id++) {
      saux = parseString(statename[id], "/");
      saux = parseString(saux[saux.size()-1], "_");

      // position code:
      /* [0] = detector, [1] = state, [2] = ybin  */
      auxnames[0].push_back(saux[0]);
      auxnames[1].push_back(saux[1]);
      auxnames[2].push_back(parseString(saux[saux.size()-1], ".")[0]);

      if(read[id] == 1) file.open(statename[id]);
      if(file.is_open()) {
	saux = parseString(statename[id], "/");
	saux = parseString(saux[saux.size()-1], "_");
	mass = QMass[saux[1]];
	sqrts = stof(saux[2])*1000;

	// get nr of pts for each stateid
	while(getline(file,data)) ns[id]+=1;
	
	file.clear();
	file.seekg(0);
	// save info for each stateid
	// position code:
	/* [0] = mass, [1] = yi, [2] = yf, [3] = pT/M, [4] = pTi/M, [5] = pTf/M
	   [6] = sigma*mass/norm, [7] = e(sigma), [8] = K(pol), [9] = sqrt(s)
	   [10] = stateid 	 */
	for(int j = 0; j < ns[id]; j++) {
	  for(int i = 0; i < 12; i++)
	    file >> nums[i];
	  
	  // aux to correct unc: TODO should make the numbers here a variable to be read off somewhere
	  if(method[id] == 1)
	    aux=nums[8]*nums[8]-1.8e-2*nums[5]*1.8e-2*nums[5];
	  
	  datasigma[0].push_back(mass);
	  datasigma[1].push_back(nums[0]);
	  datasigma[2].push_back(nums[1]);
	  datasigma[3].push_back(nums[2]/mass);
	  datasigma[4].push_back(nums[3]/mass);
	  datasigma[5].push_back(nums[4]/mass);
	  datasigma[6].push_back(nums[5]*mass/signorm[id]);

	  //three ways to register uncertainty (see method variable)
	  if(method[id] == 0)
	    datasigma[7].push_back(sqrt(nums[6]*nums[6]+nums[8]*nums[8])*mass/signorm[id]);
	  else if(method[id] == 1)
	    datasigma[7].push_back(sqrt(nums[6]*nums[6]+aux)*mass/signorm[id]);
	  else if(method[id] == 2)
	    datasigma[7].push_back(sqrt(nums[6]*nums[6]*nums[5]*nums[5]+nums[8]*nums[8]*nums[5]*nums[5])*mass/(100*signorm[id]));

	  datasigma[8].push_back(nums[10]/nums[11]-1);
	  datasigma[9].push_back(sqrts);
	  datasigma[10].push_back(id);
	}
      }
      file.close();
    }
  }
  
  //method to read the types and initializations of each fit parameter from file
  void getParams(const char* param_list) {
    ifstream fin;
    string data;
    string names[3] = {"norm", "shape", "nuis"};
    vector <string> aux;
    
    // get nr of params for each of the 3 types
    fin.open(param_list);
    while(getline(fin, data)) {
      aux = parseString(data, " ");
      for(int i = 0; i < 3; i++) 
	if (names[i].compare(aux[0])==0) 
	  nparam[i]++;
    }
    norm.resize(nparam[0]);
    shape.resize(nparam[1]); 
    nuis.resize(nparam[2]);

    fin.clear();
    fin.seekg(0);
    // for each type of par, register associated values
    // MAKE SURE PARAMETERS ARE ORDERED BY TYPE
    for(int i = 0; i < nparam[0]; i++)
      fin >> data >> norm[i].name >> norm[i].mult >> norm[i].fix;
    for(int i = 0; i < nparam[1]; i++)
      fin >> data >> shape[i].name >> shape[i].val >> shape[i].unc >> shape[i].fix;
    for(int i = 0; i < nparam[2]; i++) {
      fin >> data >> nuis[i].name >> data >> nuis[i].fix; 
      aux = parseString(data, "/");
      nuis[i].normunc = stof(aux[0]);
      if(aux.size() > 1) nuis[i].normunc /= stof(aux[1]);
    }
    fin.close();
  }

  // method to get which norm and nuis affect each state
  void getParMethod(const char* method_list) {
    norm_att.resize(states);
    nuis_att.resize(states);    
    ifstream fin;
    int nnorm, nnuis, auxi;
    
    fin.open(method_list);
    for(int id = 0; id < states; id++) {
      fin >> auxi >> nnorm;
      norm_att[id].resize(nnorm);
      for(int n = 0; n < nnorm; n++)
	fin >> norm_att[id][n];
      fin >> nnuis;
      nuis_att[id].resize(nnuis);
      for(int n = 0; n < nnuis; n++)
	fin >> nuis_att[id][n];
    }
    fin.close();
  }

  // auxiliary function that calculates the cross section for a given data point
  double csCalc(int i) {
    
    dpt = datasigma[5][i]-datasigma[4][i];
    dy = datasigma[2][i]-datasigma[1][i];
    npt = nmin + int(0.5 + dn*(dpt*QMass["psip"]-1.));
    nsteps = npt*ny;
    
    double cs = 0.;
    for(int xpt = 0; xpt < npt; xpt++) {
      sigpar[1] = datasigma[4][i]+xpt*dpt/(npt-1.)+1e-10;
      for(int xy = 0; xy < ny; xy++) {
	sigpar[2] = datasigma[1][i]+xy*dy/(ny-1.);
	cs+=sig(sigpar);
      }
    }
    cs/=(nsteps);
    return cs;
  }
  
  // fit function
  double myFunction(double *par) {
    //first calculate const normalization
    len = nparam[1] + 5;
    sigpar.resize(len);

    for(int i = 0; i < len; i++) {
      if(i < 5) sigpar[i] = 1;
      else sigpar[i] = par[i-5+nparam[0]];
    }
    
    double pred, chisq, sum = 0, ratio, lumibr;
    int counter = 0;

    // run over all pts and get chisquare
    for(int id = 0; id < states; id++) {
      //for each state calculate norm factor (product of L terms)
      sigpar[3] = 1.;
      for(int i = 0; i < norm_att[id].size(); i++){
	sigpar[3] *= par[norm_att[id][i]];
      }
      //for each state calculate lumi*br factor
      lumibr = 1.;
      for(int i = 0; i < nuis_att[id].size(); i++)
	lumibr *= par[nuis_att[id][i]+nparam[0]+nparam[1]];
      
      //for each state, assign sqrt(s)/M, M 
      sigpar[0] = datasigma[9][counter]/datasigma[0][counter];
      sigpar[4] = datasigma[0][counter];

      //cycle over each point in the state
      for(int n = 0; n < ns[id]; n++)
	if(datasigma[3][n+counter] >= ptmmin) {
	  //calculate cs prediction
	  pred = csCalc(n+counter);

	  //ratio will be !=1 if we consider polarization
	  //lambda=0.;
	  //ratio=(1+1./3*datasigma[4][i])/(1+(1-lambda)/(3+lambda)*datasigma[4][i]);
	  ratio = 1.;
	  
	  chisq = ((pred - datasigma[6][n+counter]*ratio*lumibr)/(datasigma[7][n+counter]*ratio*lumibr))*((pred - datasigma[6][n+counter]*ratio*lumibr)/(datasigma[7][n+counter]*ratio*lumibr));
	  sum += chisq;
	}
      counter += ns[id];
    }

    // add factors for the nuisance parameters
    for(int i = 0; i < nparam[2]; i++) {
      chisq = ((par[i+nparam[0]+nparam[1]] - 1) / nuis[i].normunc) * ((par[i+nparam[0]+nparam[1]] - 1) / nuis[i].normunc);
      sum += chisq;
    }

    return sum;
  }

  //estimates initialization for L value
  double lest(int init_i, int len_i)
  {
    double min_dpt = 100;
    int min_i = 0;

    for(int i = init_i; i < init_i + len_i; i++)
      if(abs(datasigma[3][i] - PTMNORM) < min_dpt) {
	min_i = i;
	min_dpt = abs(datasigma[3][i] - PTMNORM);
      }
    
    return datasigma[6][min_i];
  }
  
  // method that initializes and runs the fitting
  void doFit() {
    int npartot = nparam[0]+nparam[1]+nparam[2], aux_int = 0;
    char char_array[100];

    //define TFitter
    TFitter* fit = new TFitter(npartot);
    fit->SetFCN(minuitFunction);

    // initialize normalizations
    int len_n[2] = {0};
    double L_est[nparam[0]];
    for(int i = 0; i < 6; i++) len_n[0] += ns[i];
    for(int i = 6; i < 14; i++) len_n[1] += ns[i];
    for(int i = 0; i < 2; i++) {
      L_est[i] = lest(aux_int, len_n[i]);
      aux_int += len_n[i];
    }
    for(int i = 2; i < nparam[0]; i++) L_est[i] = 1;

    // initialize norm parameters
    for(int i = 0; i < nparam[0]; i++) {
      strcpy(char_array, norm[i].name.c_str());
      fit->SetParameter(i, char_array, L_est[i]*norm[i].mult, L_est[i]*norm[i].mult/10, 0, 0);
      if(norm[i].fix == 1) fit->FixParameter(i);
    }
    // initialize shape parameters
    for(int i = 0; i < nparam[1]; i++) {
      strcpy(char_array, shape[i].name.c_str());
      fit->SetParameter(i+nparam[0], char_array, shape[i].val, shape[i].unc, 0, 0);
      if(shape[i].fix == 1) fit->FixParameter(i+nparam[0]);
    }
    // initialize nuisance parameters
    aux_int = 0;
    for(int i = 0; i < nparam[2]; i++) {
      strcpy(char_array, nuis[i].name.c_str());
      fit->SetParameter(i+nparam[0]+nparam[1], char_array, 1, 0.1, 0, 0);
      if(nuis[i].fix == 1) fit->FixParameter(i+nparam[0]+nparam[1]);
      else aux_int++;
    }
    
    fit->ExecuteCommand("MIGRAD",0,0);
    //fit->ExecuteCommand("MIGRAD",0,0);

    fitpar.resize(npartot);
    efitpar.resize(npartot);
    for(int i = 0; i < npartot; i++) {
      fitpar[i] = fit->GetParameter(i);
      efitpar[i] = fit->GetParError(i);
    }

    double chipar[npartot];
    for(int i = 0; i < npartot; i++)
      chipar[i] = fitpar[i];
    chisquare = myFunction(chipar);
  
    ndf = 0;
    //calculating the number of data points used in the fit
    for(int i = 0; i < (int)datasigma[0].size(); i++)
      {
	if(datasigma[3][i] >= ptmmin)
	  ndf += 1;
      }
    cout << ndf << " data points" << endl;
    //I am adding one degree of freedom for each of the nuisance parameters
    ndf += aux_int;
    //the final ndf is obtained by subtracting the number of free parameters
    ndf -= fit->GetNumberFreeParameters();
    cout << "ndf " << ndf << endl;
  
    cout << "chi^2 norm : " << chisquare/ndf  << endl;
    cout << "chi^2 prob : " << TMath::Prob(chisquare,ndf)  << endl;

    //store fit results for posterior plotting if desired
    ofstream fout;
    fout.open(fitresname);
    for(int i = 0; i < npartot; i++)
      fout << fit->GetParName(i) << "   e(" << fit->GetParName(i) << ")   ";
    fout << "chi2   ndf" << endl;
    for(int i = 0; i < npartot; i++)
      fout << fit->GetParameter(i) << "   " << fit->GetParError(i) << "   ";
    fout << chisquare << "   " << ndf << endl;
    fout.close();
  }

  void readRes() {
    ifstream fin;
    string data;
    int npartot = nparam[0]+nparam[1]+nparam[2];

    // read fit results and fill (e)fitpar vectors
    fitpar.resize(npartot);
    efitpar.resize(npartot);
    
    fin.open(fitresname);
    getline(fin, data);
    for(int i = 0; i < npartot; i++)
      fin >> fitpar[i] >> efitpar[i];
    fin >> chisquare >> ndf;
    fin.close();

    // resize sigpar and fill it appropriately
    len = nparam[1] + 5;
    sigpar.resize(len);

    for(int i = 5; i < len; i++)
      sigpar[i] = fitpar[i-5+nparam[0]];
  }
  
  //method that plots everything in each canvas
  void plotcanv(int nsets, int *sets) {
    int id, counter = 0, mkrStyle;
    double lumibr, cs;
    string savename;

    // attribute axis range based on dataset
    double posif[4];
    aRange(auxnames[0][sets[0]], auxnames[1][sets[0]], posif);
    
    //get point of datasigma where current datasets start
    for(int i = 0; i < sets[0]; i++) counter += ns[i];
    
    // plotting starts here
    TCanvas *c = new TCanvas("title", "name", 700, 700);
    c->SetLogy();
    TH1F *fc = c->DrawFrame(posif[0], posif[1], posif[2], posif[3]);
    fc->SetXTitle("p_{T}/M");
    fc->SetYTitle("d#sigma / d#xidy (nb/GeV)");
    fc->GetYaxis()->SetTitleOffset(1);
    c->Modified();
    c->SetTitle("");

    // define data graphs and fit functions
    TGraphAsymmErrors **gd = new TGraphAsymmErrors*[nsets];
    TGraphAsymmErrors **gp = new TGraphAsymmErrors*[nsets];
    TF1 **f = new TF1*[nsets];
    
    //cycle over all states
    for(int i_set = 0; i_set < nsets; i_set++) {
      id = sets[i_set];
      // get normalization and lumi*br factors
      sigpar[3] = 1.;
      for(int i = 0; i < norm_att[id].size(); i++)
	sigpar[3] *= fitpar[norm_att[id][i]];
      lumibr = 1.;
      for(int i = 0; i < nuis_att[id].size(); i++)
	lumibr *= fitpar[nuis_att[id][i]+nparam[0]+nparam[1]];

      // get sqrt(s) / M and MQQ
      sigpar[0] = datasigma[9][counter]/datasigma[0][counter];
      sigpar[4] = datasigma[0][counter];
      
      float datapts[5][ns[id]];
      float pullpts[2][ns[id]];
      
      // cycle over all points of [id]-state, to plot data and later pulls
      for(int j = 0; j < ns[id]; j++) {
	// filling arrays for data plotting
	datapts[0][j] = avgptm(datasigma[4][j+counter], datasigma[5][j+counter], datasigma[1][j+counter], datasigma[2][j+counter], sigpar);
	datapts[1][j] = datasigma[6][j+counter]*lumibr;
	datapts[2][j] = datasigma[7][j+counter]*lumibr;
	datapts[3][j] = datapts[0][j] - datasigma[4][j+counter];
	datapts[4][j] = datasigma[5][j+counter] - datapts[0][j];

	// filling arrays for pulls plotting
	cs = csCalc(j+counter);
	pullpts[0][j] = (datapts[1][j] - cs) / datapts[2][j];
	pullpts[1][j] = 0.;
      }

      mkrStyle = getStyle(auxnames[0][id]);
      
      // data plots defined and drawn
      gd[i_set] = new TGraphAsymmErrors(ns[id], datapts[0], datapts[1], datapts[3], datapts[4], datapts[2], datapts[2]);
      gd[i_set]->SetLineColor(i_set+1);
      gd[i_set]->SetMarkerColor(i_set+1);
      gd[i_set]->SetMarkerStyle(mkrStyle);
      gd[i_set]->SetMarkerSize(.75);
      gd[i_set]->Draw("P");

      // pulls plots defined (but not drawn)
      // all pts below ptmmin aren't drawn
      gp[i_set] = new TGraphAsymmErrors(ns[id], datapts[0], pullpts[0], datapts[3], datapts[4], pullpts[1], pullpts[1]);
      gp[i_set]->SetLineColor(i_set+1);
      gp[i_set]->SetMarkerColor(i_set+1);
      gp[i_set]->SetMarkerStyle(mkrStyle);
      gp[i_set]->SetMarkerSize(.75);
      for(int i_data = 0; i_data < ns[id]; i_data++)
	if(datapts[0][i_data] < ptmmin)
	  gp[i_set]->RemovePoint(0);
      
      // fit model can't be made as flexible regarding param number
      // TODO think a bit abt how (if?) this could be improved
      f[i_set] = new TF1("cs fit", "sigplot([0], x, [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12], [13])", ptmmin, 49.9);
      f[i_set]->SetParameter(0, sigpar[0]);
      f[i_set]->SetParameter(1, (datasigma[1][counter]+datasigma[2][counter])/2);
      for(int j=2; j<len-1; j++)
	f[i_set]->SetParameter(j, sigpar[j+1]);
      f[i_set]->SetLineColor(i_set+1);
      f[i_set]->Draw("lsame");

      counter += ns[id]; 
    }

    // text on the plot
    TLatex lc;
    double xpos = getPos(posif[0], posif[2], 0.625, 0);
    lc.SetTextSize(0.03);
    lc.DrawLatex(xpos, getPos(posif[3], posif[1], 0.5*(nsets+2)/8, 1), Form("#chi^{2}/ndf = %.0f/%d", chisquare, ndf));
    lc.DrawLatex(xpos, getPos(posif[3], posif[1], 0.5*(nsets+3)/8, 1), Form("P(#chi^{2},ndf) = %.1f%%", 100*TMath::Prob(chisquare, ndf)));
    xpos = getPos(posif[0], posif[2], 1./20, 0);
    lc.DrawLatex(xpos, getPos(posif[1], posif[3], 1./20, 1), Form("pp %.0f TeV", datasigma[9][counter-1]/1000.));

    // draw legend
    double endpt = 0.85 - (double)nsets*0.05;  
    TLegend *leg = new TLegend(0.6, endpt, 0.9, 0.9);
    leg->SetTextSize(0.03);
    for(int i = 0; i < nsets; i++)  {
      savename = auxnames[0][sets[i]]+" "+QName[auxnames[1][sets[i]]]+" "+auxnames[2][sets[i]];
      const char *st = savename.c_str();
      leg->AddEntry(gd[i], st, "p");
    }
    leg->Draw();
      
    // save data plot
    savename = "plots/"+auxnames[1][sets[0]]+"_"+Form("%.0f", datasigma[9][counter-1]/1000.)+"_"+auxnames[0][sets[0]]+"_cs.pdf";
    const char *st = savename.c_str();
    c->SaveAs(st);

    // redo plotting for pulls
    c->Clear();
    c->SetLogy(0);

    TH1F *fp = c->DrawFrame(posif[0], -10, posif[2], 10);
    fp->SetXTitle("p_{T}/M");
    fp->SetYTitle("pulls");
    fp->GetYaxis()->SetTitleOffset(1);
    c->Modified();
    c->SetTitle("");
   
    //plot line at zero
    TF1 *zero = new TF1("zero", "0", posif[0], posif[2]);
    zero->SetLineColor(kBlue);
    zero->SetLineStyle(7);
    zero->Draw("lsame");

    // plot pt/M cutoff
    TLine *ptm = new TLine(ptmmin, -10, ptmmin, 10);
    ptm->SetLineStyle(7);
    ptm->Draw("lsame");
  
    //draw pulls
    for(int i = 0; i < nsets; i++) {
      gp[i]->Draw("P"); 
    }

    //text on the plot
    TLatex lp;
    lp.SetTextSize(0.03);
    lp.DrawLatex(xpos, getPos(10, -10, 1.5/20, 0), Form("#chi^{2}/ndf = %.0f/%d", chisquare, ndf));
    lp.DrawLatex(xpos, getPos(10, -10, 3./20, 0), Form("P(#chi^{2},ndf) = %.1f%%", 100*TMath::Prob(chisquare, ndf)));
    lp.DrawLatex(xpos, getPos(-10, 10, 1./20, 0), Form("pp %.0f TeV", datasigma[9][counter-1]/1000.));
    
    //draw legend
    TLegend *legp = new TLegend(0.6, endpt, 0.9, 0.9);
    legp->SetTextSize(0.03);
    for(int i = 0; i < nsets; i++)  {
      savename = auxnames[0][sets[i]]+" "+QName[auxnames[1][sets[i]]]+" "+auxnames[2][sets[i]];
      const char *st = savename.c_str();
      legp->AddEntry(gp[i], st, "p");
    }
    legp->Draw();
    
    //save pulls plot
    savename = "plots/"+auxnames[1][sets[0]]+"_"+Form("%.0f", datasigma[9][counter-1]/1000.)+"_"+auxnames[0][sets[0]]+"_cs_pulls.pdf";
    const char* savep = savename.c_str();
    c->SaveAs(savep);
    
    c->Destructor();
  }
  
  //method that does the plotting
  void doPlot(const char* state_div) {
    gROOT->SetBatch();
    
    ifstream fin;
    ofstream tex;
    int nsets, nval;

    int npartot = nparam[0]+nparam[1]+nparam[2];
    vector <string> names = {"L_{J/\\psi}", "L_{\\Upsilon(1S)}", "L_{LHCb,\\psi}(y1)", "L_{LHCb,\\psi}(y2)", "L_{LHCb,\\psi}(y3)", "L_{LHCb,\\psi}(y4)", "L_{LHCb,\\psi}(y5)","L_{LHCb,\\Upsilon}(y1)", "L_{LHCb,\\Upsilon}(y2)", "L_{LHCb,\\Upsilon}(y3)", "L_{LHCb,\\Upsilon}(y4)", "\\beta", "\\tau", "\\rho", "\\delta", "\\gamma_1", "\\gamma_2", "b", "c", "d", "e", "BR_{jpsidm}", "BR_{ups1dm}", "\\mathcal L_{CMS}", "\\mathcal L_{ATLAS}", "\\mathcal L_{LHCb}(\\psi)", "\\mathcal L_{LHCb}(\\Upsilon)"};
    
    // make latex file with fit parameters
    tex.open("plots/fitp.tex");

    tex << "\\begin{table}" << endl;
    tex << "\\centering" << endl;
    tex << "\\begin{tabular}{c|c|c}" << endl;
    tex << "Parameter & Value & Uncertainty \\\\" << endl;
    tex << "\\hline" << endl;
    for(int i = 0; i < nparam[0]; i++) {
      tex << "$" << names[i] << "$ & ";
      tex << fitpar[i] << " & ";
      if(efitpar[i] == 0) tex << "fixed \\\\" << endl;
      else tex << efitpar[i] << " \\\\" << endl;
    }
    for(int i = 0; i < nparam[1]; i++) {
      tex << "$" << names[i+nparam[0]] << "$ & ";
      tex << fitpar[i+nparam[0]] << " & ";
      if(efitpar[i+nparam[0]] == 0) tex << "fixed \\\\" << endl;
      else tex << efitpar[i+nparam[0]] << " \\\\" << endl;
    }
    for(int i = 0; i < nparam[2]; i++) {
      tex << "$" << names[i+nparam[0]+nparam[1]] << "$ & ";
      tex << fitpar[i+nparam[0]+nparam[1]] << " & ";
      if(efitpar[i+nparam[0]+nparam[1]] == 0) tex << "fixed \\\\" << endl;
      else tex << efitpar[i+nparam[0]+nparam[1]] << " \\\\" << endl;
      }
    tex << "min $p_T/M$ & " << ptmmin << " & fixed \\\\" << endl;
    tex << "\\end{tabular}" << endl;
    tex << "\\caption{Fit parameters ($\\chi^2$ / ndf = " << chisquare << " / " << ndf << " = " << chisquare/ndf << ")}" << endl;
    tex << "\\end{table}" << endl;
    tex.close();
    
    // cycle over each line of the plotting division file
    fin.open(state_div);
    while(1) {
      nval = 0;
      fin >> nsets;
      int sets[nsets];

      for(int i = 0; i < nsets; i++) {
	fin >> sets[i];
	nval += ns[sets[i]];
      }

      // plot only if there are no points in dataset being considered
      if(nval != 0) plotcanv(nsets, sets);

      //break cycle when we get to the last dataset
      if(sets[nsets-1] == states-1)
	break;
    }
    fin.close();
  }
};

// initialization of fit model outside main so it can be accessed by the tfitter wrapper
FitModel fitclass;

//Minuit-specific "wrapping"
void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg)
{
  result = fitclass.myFunction(par);
}
