// function to parse a string into components separated by "deli"
vector< string > parseString( string line, string deli) {
  vector< string > out;
  string aux = line;
  size_t pos = 0;
  
  while( (pos = aux.find(deli)) != string::npos)
    {
      out.push_back( aux.substr( 0, pos ) );
      aux.erase( 0, pos + deli.length() );	    
    }
  out.push_back( aux );
  
  return out;
}


// comparing the fit results 
void plotComp()
{
  const int n_inp = 3;
  const int n_state = 7;

  double chinorm[n_inp], chiProb[n_inp], chisquare[n_inp], ndf[n_inp];
  double L_QQ[n_state][n_inp], eL_QQ[n_state][n_inp];

  double beta[n_state][n_inp], e_beta[n_state][n_inp];
  
  string parL[] = {"check_0_universal", "check_1_beta", "check_2_fbeta"};
  ifstream fin;
  string data;
  string names[3] = {"norm", "shape", "nuis"};
  vector <string> aux;
  for(int i = 0; i < n_inp; i++) {
  
    const char* param_list = Form("%s/param_list.txt", parL[i].c_str());
    
    int nparam[4] = {0};
    vector <double> fitpar, efitpar;

    // get nr of params for each of the 3 types
    fin.open(param_list);
    while(getline(fin, data)) {
      aux = parseString(data, " ");
      for(int i = 0; i < 3; i++) 
	if (names[i].compare(aux[0])==0) 
	  nparam[i]++;
    }
    if(i!=1)
      nparam[3] = n_state;
    fin.close();
    
    ifstream fin2;
    int npartot = nparam[0]+nparam[1]+nparam[2]+nparam[3];

    // read fit results and fill (e)fitpar vectors
    fitpar.resize(npartot);
    efitpar.resize(npartot);
    
    // fill the arrays of results we want
    string fitresname = Form("%s/fit.txt", parL[i].c_str());
  
    fin2.open(fitresname);
    getline(fin2, data);
    for(int j = 0; j < npartot; j++)
      fin2 >> fitpar[j] >> efitpar[j];
    fin2 >> chisquare[i] >> ndf[i];
    fin2.close();

    cout << chisquare[i] << "  "<< ndf[i] << endl;
    
    chinorm[i] = chisquare[i]/ndf[i];
    chiProb[i] = TMath::Prob(chisquare[i], ndf[i]);
    for(int j = 0; j < n_state; j++) {
      L_QQ[j][i] = fitpar[j];
      eL_QQ[j][i] = efitpar[j];
      // beta value depends on case
      // 0 -> universal beta
      if(i==0) {
	if(j==1) {
	  beta[j][i] = fitpar[nparam[0]];
	  e_beta[j][i] = efitpar[nparam[0]];
	}
	else { 
	  beta[j][i] = 0;
	  e_beta[j][i] = 0;
	}
      }
      // 1 -> beta per state
      else if(i==1) {
	if(j!=1 && j != 2) {
	  beta[j][i] = fitpar[nparam[0]+2+j];
	  e_beta[j][i] = efitpar[nparam[0]+2+j];
	}
	else { 
	  beta[j][i] = 0;
	  e_beta[j][i] = 0;
	}
      }
      // 2 -> just the two betas
      else if(i==2) {
	if(j == 1 || j == 2) {
	  beta[j][i] = fitpar[nparam[0]+j-1];
	  e_beta[j][i] = efitpar[nparam[0]+j-1];
	}
      }
      else { 
	beta[j][i] = 0;
	e_beta[j][i] = 0;
      }
    }

  }
  
  // make the TGraphs from the results
  string pdf_n[] = {"single #beta", "state #beta", "two #beta"};

  TH1F *h_chi = new TH1F("h_chi", "#chi^{2}/ndf per case", n_inp, 0, n_inp);
  TH1F *h_chiP = new TH1F("h_chiP", "P(#chi^{2},ndf) per case", n_inp, 0, n_inp);
 
  TH1F **h_LQQ = new TH1F*[n_state];
  TH1F **h_beta = new TH1F*[n_state];
  for(int i = 0; i < n_state; i++) {
    h_LQQ[i] = new TH1F(Form("h_LQQ_%d", i), "L_{QQ} per case", n_inp, 0, n_inp);
    h_beta[i] = new TH1F(Form("h_beta_%d",i), "#beta per case", n_inp, 0, n_inp);
  }

  for(int i = 0; i < n_inp; i++) {
    h_chi->Fill(pdf_n[i].c_str(), chinorm[i]);
    h_chi->SetBinError(i+1, 0);

    h_chiP->Fill(pdf_n[i].c_str(), chiProb[i]);
    h_chiP->SetBinError(i+1, 0);

    
    for(int j = 0; j < n_state; j++) {
      h_LQQ[j]->Fill(pdf_n[i].c_str(), L_QQ[j][i]);
      h_LQQ[j]->SetBinError(i+1, eL_QQ[j][i]);

      h_beta[j]->Fill(pdf_n[i].c_str(), beta[j][i]);
      h_beta[j]->SetBinError(i+1, e_beta[j][i]);
    }
  }
    
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);

  h_chi->SetStats(0);
  h_chi->GetXaxis()->SetTitle("case");
  h_chi->GetYaxis()->SetTitle("#chi^{2}/ndf");
  h_chi->GetYaxis()->SetTitleOffset(1.3);
  h_chi->GetYaxis()->SetRangeUser(0, 2);
  h_chi->LabelsDeflate("X");
  h_chi->SetMarkerStyle(20);
  h_chi->SetLineColor(kBlack);
  h_chi->SetMarkerColor(kBlack);
  h_chi->Draw("P");

  TLine *l_z = new TLine(0, 1, n_inp, 1);
  l_z->SetLineStyle(kDashed);
  l_z->SetLineColor(kBlack);
  l_z->Draw();

  TLatex lc;
  lc.SetTextSize(0.03);
  for(int i = 0; i < n_inp; i++) {
    lc.DrawLatex(i+0.1, h_chi->GetBinContent(i+1)-0.1, Form("%.0f/%.0f", chisquare[i], ndf[i]));
  }  

  c->SaveAs("chi2.pdf");
  c->Clear();


  TLegend *leg = new TLegend(0.46, 0.65, 0.61, 0.9);
  leg->SetTextSize(0.03);
  
  int ct = 1;
  int c_col[] = {kBlack, kRed, kGreen+2, kBlue, kYellow+1};
  string c_lbl[] = {"J/#psi", "#psi(2S)", "#Upsilon(1S)", "#Upsilon(2S)", "#Upsilon(3S)"};
  h_beta[0]->SetStats(0);
  h_beta[0]->GetXaxis()->SetTitle("case");
  h_beta[0]->GetYaxis()->SetTitle("#beta");
  h_beta[0]->GetYaxis()->SetTitleOffset(1.3);
  h_beta[0]->GetYaxis()->SetRangeUser(2.0, 3.5);
  h_beta[0]->LabelsDeflate("X");
  h_beta[0]->SetMarkerStyle(20);
  h_beta[0]->SetMarkerSize(.75);
  h_beta[0]->SetLineColor(kBlack);
  h_beta[0]->SetMarkerColor(kBlack);
  h_beta[0]->Draw("error");
  leg->AddEntry(h_beta[0], c_lbl[0].c_str(), "pl");

  h_beta[1]->SetMarkerStyle(20);
  h_beta[1]->SetMarkerSize(.75);
  h_beta[1]->SetLineColor(c_col[0]);
  h_beta[1]->SetLineStyle(kDashed);
  h_beta[1]->SetMarkerColor(c_col[0]);
  h_beta[1]->Draw("same");
  h_beta[2]->SetMarkerStyle(20);
  h_beta[2]->SetMarkerSize(.75);
  h_beta[2]->SetLineStyle(kDashed);
  h_beta[2]->SetLineColor(c_col[0]);
  h_beta[2]->SetMarkerColor(c_col[0]);
  h_beta[2]->Draw("same");
  
  for(int i = 3; i < n_state; i++) {
    h_beta[i]->SetMarkerStyle(20);
    h_beta[i]->SetMarkerSize(.75);
    h_beta[i]->SetLineColor(c_col[ct]);
    h_beta[i]->SetMarkerColor(c_col[ct]);
    h_beta[i]->Draw("same");

    leg->AddEntry(h_beta[i], c_lbl[ct].c_str(), "pl");
    ct++;
  }

  leg->Draw();

  c->SaveAs("beta.pdf");
  c->Clear();

  h_LQQ[0]->SetStats(0);
  h_LQQ[0]->GetXaxis()->SetTitle("case");
  h_LQQ[0]->GetYaxis()->SetTitle("L_{QQ}");
  h_LQQ[0]->GetYaxis()->SetTitleOffset(1.3);
  h_LQQ[0]->GetYaxis()->SetRangeUser(0, 2.6);
  h_LQQ[0]->LabelsDeflate("X");
  h_LQQ[0]->SetMarkerStyle(20);
  h_LQQ[0]->SetMarkerSize(.75);
  h_LQQ[0]->SetLineColor(kBlack);
  h_LQQ[0]->SetMarkerColor(kBlack);
  h_LQQ[0]->Draw("error");

  ct = 1;
  for(int i = 3; i < n_state; i++) {
    h_LQQ[i]->SetMarkerStyle(20);
    h_LQQ[i]->SetMarkerSize(.75);
    h_LQQ[i]->SetLineColor(c_col[ct]);
    h_LQQ[i]->SetMarkerColor(c_col[ct]);
    h_LQQ[i]->Draw("same");

    ct++;
  }

  leg->Draw();
  
  c->SaveAs("L_QQ.pdf");
  c->Clear();

  c->SetLogy();

  h_chiP->SetStats(0);
  h_chiP->GetXaxis()->SetTitle("case");
  h_chiP->GetYaxis()->SetTitle("P(#chi^{2},ndf)");
  h_chiP->GetYaxis()->SetTitleOffset(1.3);
  h_chiP->GetYaxis()->SetRangeUser(1e-25, 1);
  h_chiP->LabelsDeflate("X");
  h_chiP->SetMarkerStyle(20);
  h_chiP->SetLineColor(kBlack);
  h_chiP->SetMarkerColor(kBlack);
  h_chiP->Draw("P");

  for(int i = 0; i < n_inp; i++) {
    lc.DrawLatex(i+0.5, h_chiP->GetBinContent(i+1)*2, Form("%.2e", chiProb[i]));
  }  

  c->SaveAs("chiP.pdf");
  c->Clear();

 
  c->Destructor();
}
