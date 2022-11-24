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
  const int n_inp = 6;
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
  double xi_cut[] = {1.5, 2, 2.5, 3, 3.5, 4};
  TGraph *g_chi = new TGraph(n_inp, xi_cut, chinorm);
  TGraph *g_chiP = new TGraph(n_inp, xi_cut, chiProb);
  TGraphErrors *g_beta = new TGraphErrors(n_inp, xi_cut, beta2, xi_err, e_beta2);
  TGraphErrors **g_LQQ = new TGraphErrors*[n_state];
  TGraphErrors **g_fQQ = new TGraphErrors*[n_state];
  for(int i = 0; i < n_state; i++) {
    g_LQQ[i] = new TGraphErrors(n_inp, xi_cut, L_QQ[i], xi_err, eL_QQ[i]);
    g_fQQ[i] = new TGraphErrors(n_inp, xi_cut, f_QQ[i], xi_err, ef_QQ[i]);
  }
  
  TCanvas *c = new TCanvas("", "", 700, 700);
  c->SetRightMargin(0.03);
  
  TH1F *fc = c->DrawFrame(1, 0, 4.5, 5);
  fc->SetTitle("#chi^{2}/ndf vs #xi cutoff");
  fc->SetXTitle("#xi cutoff");
  fc->SetYTitle("#chi^{2}/ndf");
  fc->GetYaxis()->SetTitleOffset(1.3);

  g_chi->SetMarkerStyle(20);
  g_chi->Draw("psame");

  TLine *l_z = new TLine(1, 1, 4.5, 1);
  l_z->SetLineStyle(kDashed);
  l_z->SetLineColor(kBlack);
  l_z->Draw();

  TLatex lc;
  lc.SetTextSize(0.03);
  for(int i = 0; i < g_chi->GetN(); i++) {
    lc.DrawLatex(g_chi->GetX()[i], g_chi->GetY()[i]+0.1, Form("%.0f/%.0f", chisquare[i], ndf[i]));
  }  

  c->SaveAs("chi2.pdf");
  c->Clear();

  TH1F *fc2 = c->DrawFrame(1, 2.5, 4.5, 3.5);
  fc2->SetTitle("#beta_{2} vs #xi cutoff");
  fc2->SetXTitle("#xi cutoff");
  fc2->SetYTitle("#beta_{2}");
  fc2->GetYaxis()->SetTitleOffset(1.3);

  g_beta->SetMarkerStyle(20);
  g_beta->SetMarkerSize(.75);
  g_beta->Draw("psame");

  c->SaveAs("beta.pdf");
  c->Clear();

  TH1F *fc3 = c->DrawFrame(1, 0, 4.5, 0.55);
  fc3->SetTitle("f_{#beta_{1}} vs #xi cutoff");
  fc3->SetXTitle("#xi cutoff");
  fc3->SetYTitle("f_{#beta_{1}}");
  fc3->GetYaxis()->SetTitleOffset(1);

  TLegend *leg = new TLegend(0.82, 0.65, 0.97, 0.9);
  leg->SetTextSize(0.03);
  
  int ct = 0;
  int c_col[] = {kBlack, kRed, kGreen+2, kBlue, kYellow+1};
  string c_lbl[] = {"J/#psi", "#psi(2S)", "#Upsilon(1S)", "#Upsilon(2S)", "#Upsilon(3S)"};
  for(int i = 0; i < n_state; i++)
    if(i != 1 && i != 2) {
      g_fQQ[i]->SetMarkerStyle(20);
      g_fQQ[i]->SetMarkerSize(.75);
      g_fQQ[i]->SetMarkerColor(c_col[ct]);
      g_fQQ[i]->SetLineColor(c_col[ct]);
      g_fQQ[i]->Draw("psame");
      leg->AddEntry(g_fQQ[i], c_lbl[ct].c_str(), "pl");
      ct++;
    }

  leg->Draw();

  c->SaveAs("f_beta.pdf");
  c->Clear();

  TH1F *fc4 = c->DrawFrame(1, 0, 4.5, 2.6);
  fc4->SetTitle("L_{QQ} vs #xi cutoff");
  fc4->SetXTitle("#xi cutoff");
  fc4->SetYTitle("L_{QQ}");
  fc4->GetYaxis()->SetTitleOffset(1.3);

  ct = 0;
  for(int i = 0; i < n_state; i++)
    if(i != 1 && i != 2) {
      g_LQQ[i]->SetMarkerStyle(20);
      g_LQQ[i]->SetMarkerSize(.75);
      g_LQQ[i]->SetMarkerColor(c_col[ct]);
      g_LQQ[i]->SetLineColor(c_col[ct]);
      g_LQQ[i]->Draw("psame");
      ct++;
    }

  leg->Draw();
  
  c->SaveAs("L_QQ.pdf");
  c->Clear();

  c->SetLogy();
  TH1F *fc1 = c->DrawFrame(1, 1e-51, 4.5, 1);
  fc1->SetTitle("P(#chi^{2},ndf) vs #xi cutoff");
  fc1->SetXTitle("#xi cutoff");
  fc1->SetYTitle("P(#chi^{2},ndf)");
  fc1->GetYaxis()->SetTitleOffset(1.3);

  g_chiP->SetMarkerStyle(20);
  g_chiP->Draw("psame");

  for(int i = 0; i < g_chiP->GetN(); i++) {
    lc.DrawLatex(g_chiP->GetX()[i]-0.5, g_chiP->GetY()[i]*0.5, Form("%.2e", chiProb[i]));
  }  

  
  c->SaveAs("chiP.pdf");
  c->Clear();

 
  c->Destructor();
}
