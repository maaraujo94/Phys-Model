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
  const int n_inp = 12;
  const int n_state = 7;

  const char* param_list = "check_0/param_list.txt"; //fit parameters. set last value to 1 to fix in the fit. norm parameters have relative uncertainty, shape have central value and uncertainty, nuis have uncertainty

  ifstream fin;
  string data;
  string names[3] = {"norm", "shape", "nuis"};
  vector <string> aux;
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
  nparam[3] = n_state;
  fin.close();
  
  ifstream fin2;
  int npartot = nparam[0]+nparam[1]+nparam[2]+nparam[3];

  // read fit results and fill (e)fitpar vectors
  fitpar.resize(npartot);
  efitpar.resize(npartot);
  
  double chinorm[n_inp], chiProb[n_inp], chisquare[n_inp], ndf[n_inp], beta2[n_inp], e_beta2[n_inp], xi_err[n_inp];
  double L_QQ[n_state][n_inp], eL_QQ[n_state][n_inp], f_QQ[n_state][n_inp], ef_QQ[n_state][n_inp];

  // fill the arrays of results we want
  for(int i = 0; i < n_inp; i++) {
    string fitresname = Form("check_%d/fit.txt", i);
  
    fin2.open(fitresname);
    getline(fin2, data);
    for(int j = 0; j < npartot; j++)
      fin2 >> fitpar[j] >> efitpar[j];
    fin2 >> chisquare[i] >> ndf[i];
    fin2.close();

    chinorm[i] = chisquare[i]/ndf[i];
    chiProb[i] = TMath::Prob(chisquare[i], ndf[i]);
    beta2[i] = fitpar[nparam[0]+1];
    e_beta2[i] = efitpar[nparam[0]+1];
    for(int j = 0; j < n_state; j++) {
      L_QQ[j][i] = fitpar[j];
      f_QQ[j][i] = fitpar[nparam[0]+nparam[1]+nparam[2]+j];
      eL_QQ[j][i] = efitpar[j];
      ef_QQ[j][i] = efitpar[nparam[0]+nparam[1]+nparam[2]+j];
    }
    xi_err[i] = 0.25;
  }

  // make the TGraphs from the results  
  TH1F *h_chi = new TH1F("h_chi", "#chi^{2}/ndf per g(cos#alpha)", n_inp, 0, n_inp);
  TH1F *h_chiP = new TH1F("h_chiP", "P(#chi^{2},ndf) per g(cos#alpha)", n_inp, 0, n_inp);
  TH1F *h_beta = new TH1F("h_beta", "#beta_{2} per g(cos#alpha)", n_inp, 0, n_inp);

  TH1F **h_LQQ = new TH1F*[n_state];
  TH1F **h_fQQ = new TH1F*[n_state];
  for(int i = 0; i < n_state; i++) {
    h_LQQ[i] = new TH1F(Form("h_LQQ_%d", i), "L_{QQ} per g(cos#alpha)", n_inp, 0, n_inp);
    h_fQQ[i] = new TH1F(Form("h_fQQ_%d", i), "f_{#beta_{1}} per g(cos#alpha)", n_inp, 0, n_inp);
  }

  string cosa_n[n_inp];
  int rhov[] = {2,3}, delv[] = {0,1};
  for(int i = 0; i < n_inp; i++) {
    cosa_n[i] = Form("g_{%d} (%d, %d)", i/4+1, rhov[(i/2)%2], delv[i%2]);
    
    h_chi->Fill(cosa_n[i].c_str(), chinorm[i]);
    h_chi->SetBinError(i+1, 0);

    h_chiP->Fill(cosa_n[i].c_str(), chiProb[i]);
    h_chiP->SetBinError(i+1, 0);

    h_beta->Fill(cosa_n[i].c_str(), beta2[i]);
    h_beta->SetBinError(i+1, e_beta2[i]);
    
    for(int j = 0; j < n_state; j++) {
      h_LQQ[j]->Fill(cosa_n[i].c_str(), L_QQ[j][i]);
      h_LQQ[j]->SetBinError(i+1, eL_QQ[j][i]);

      h_fQQ[j]->Fill(cosa_n[i].c_str(), f_QQ[j][i]);
      h_fQQ[j]->SetBinError(i+1, ef_QQ[j][i]);
    }
  }
    
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);

  h_chi->SetStats(0);
  h_chi->GetXaxis()->SetTitle("g(cos#alpha) form");
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

  TLine *l_11 = new TLine(4, h_chi->GetMinimum(), 4, h_chi->GetMaximum());
  l_11->SetLineStyle(kDashed);
  l_11->SetLineColor(kBlack);
  l_11->Draw();
  TLine *l_12 = new TLine(8, h_chi->GetMinimum(), 8, h_chi->GetMaximum());
  l_12->SetLineStyle(kDashed);
  l_12->SetLineColor(kBlack);
  l_12->Draw();

  TLatex lc;
  lc.SetTextSize(0.03);
  for(int i = 0; i < n_inp; i++) {
    //lc.DrawLatex(i+0.5, h_chi->GetBinContent(i+1)+0.1, Form("%.0f/%.0f", chisquare[i], ndf[i]));
  }  

  c->SaveAs("chi2.pdf");
  c->Clear();

  h_beta->SetStats(0);
  h_beta->GetXaxis()->SetTitle("g(cos#alpha) form");
  h_beta->GetYaxis()->SetTitle("#beta_{2}");
  h_beta->GetYaxis()->SetTitleOffset(1.3);
  h_beta->GetYaxis()->SetRangeUser(2.5, 3.5);
  h_beta->LabelsDeflate("X");
  h_beta->SetMarkerStyle(20);
  h_beta->SetMarkerSize(.75);
  h_beta->SetLineColor(kBlack);
  h_beta->SetMarkerColor(kBlack);
  h_beta->Draw("error");

  TLine *l_21 = new TLine(4, h_beta->GetMinimum(), 4, h_beta->GetMaximum());
  l_21->SetLineStyle(kDashed);
  l_21->SetLineColor(kBlack);
  l_21->Draw();
  TLine *l_22 = new TLine(8, h_beta->GetMinimum(), 8, h_beta->GetMaximum());
  l_22->SetLineStyle(kDashed);
  l_22->SetLineColor(kBlack);
  l_22->Draw();

  c->SaveAs("beta.pdf");
  c->Clear();

  TLegend *leg = new TLegend(0.82, 0.65, 0.97, 0.9);
  leg->SetTextSize(0.03);
  
  int ct = 1;
  int c_col[] = {kBlack, kRed, kGreen+2, kBlue, kYellow+1};
  string c_lbl[] = {"J/#psi", "#psi(2S)", "#Upsilon(1S)", "#Upsilon(2S)", "#Upsilon(3S)"};
  h_fQQ[0]->SetStats(0);
  h_fQQ[0]->GetXaxis()->SetTitle("g(cos#alpha) form");
  h_fQQ[0]->GetYaxis()->SetTitle("f_{#beta_{1}}");
  h_fQQ[0]->GetYaxis()->SetTitleOffset(1);
  h_fQQ[0]->GetYaxis()->SetRangeUser(0, 0.55);
  h_fQQ[0]->LabelsDeflate("X");
  h_fQQ[0]->SetMarkerStyle(20);
  h_fQQ[0]->SetMarkerSize(.75);
  h_fQQ[0]->SetLineColor(kBlack);
  h_fQQ[0]->SetMarkerColor(kBlack);
  h_fQQ[0]->Draw("error");
  leg->AddEntry(h_fQQ[0], c_lbl[0].c_str(), "pl");

  for(int i = 3; i < n_state; i++) {
    h_fQQ[i]->SetMarkerStyle(20);
    h_fQQ[i]->SetMarkerSize(.75);
    h_fQQ[i]->SetLineColor(c_col[ct]);
    h_fQQ[i]->SetMarkerColor(c_col[ct]);
    h_fQQ[i]->Draw("same");

    leg->AddEntry(h_fQQ[i], c_lbl[ct].c_str(), "pl");
    ct++;
  }

  leg->Draw();
  
  TLine *l_31 = new TLine(4, h_fQQ[0]->GetMinimum(), 4, h_fQQ[0]->GetMaximum());
  l_31->SetLineStyle(kDashed);
  l_31->SetLineColor(kBlack);
  l_31->Draw();
  TLine *l_32 = new TLine(8, h_fQQ[0]->GetMinimum(), 8, h_fQQ[0]->GetMaximum());
  l_32->SetLineStyle(kDashed);
  l_32->SetLineColor(kBlack);
  l_32->Draw();

  c->SaveAs("f_beta.pdf");
  c->Clear();

  h_LQQ[0]->SetStats(0);
  h_LQQ[0]->GetXaxis()->SetTitle("g(cos#alpha) form");
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

  TLine *l_41 = new TLine(4, h_LQQ[0]->GetMinimum(), 4, h_LQQ[0]->GetMaximum());
  l_41->SetLineStyle(kDashed);
  l_41->SetLineColor(kBlack);
  l_41->Draw();
  TLine *l_42 = new TLine(8, h_LQQ[0]->GetMinimum(), 8, h_LQQ[0]->GetMaximum());
  l_42->SetLineStyle(kDashed);
  l_42->SetLineColor(kBlack);
  l_42->Draw();

  c->SaveAs("L_QQ.pdf");
  c->Clear();

  c->SetLogy();

  h_chiP->SetStats(0);
  h_chiP->GetXaxis()->SetTitle("g(cos#alpha) form");
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
