using namespace PlotUtils;
using namespace std;

double Chi2(MnvH1D* m1, MnvH1D* m2){
  if (m1->GetNbinsX() != m2->GetNbinsX()){
    cout << "Histos to compare don't have the same number of bins! Returning -999!!!" << endl;
    return -999.0;
  }
  TH1D* h1 = (TH1D*)m1->GetCVHistoWithError().Clone();
  TH1D* h2 = (TH1D*)m2->GetCVHistoWithError().Clone();
  double chi2 = 0.0;
  for (int whichBin = 1; whichBin <= h1->GetNbinsX(); ++whichBin){
    double h1Content = h1->GetBinContent(whichBin);
    double h2Content = h2->GetBinContent(whichBin);
    double h2Err = h2->GetBinError(whichBin);
    double diff = h1Content-h2Content;
    if(h2Err > 1e-10) chi2 += (diff*diff)/(h2Err*h2Err);
  }

  delete h1;
  delete h2;
  return chi2;
}

void PlotCHProjX(TFile* dataFile, TFile* mcFile, TString USDSIn, TString region, TString mat, bool proj=true, int loBin=0, int hiBin=-1){
  if (region != "signal" && region != "sideband"){
    cout << "Not a valid region. Exiting" << endl;
    return;
  }

  if (mat != "Fe" && mat != "Pb" && mat != "C" && mat != "Water"){
    cout << "Not a valid material. Exiting" << endl;
    return;
  }

  if (USDSIn != "US" && USDSIn != "DS" && USDSIn != "Sum"){
    cout << "Not a valid plastic region. Exiting" << endl;
    return;
  }

  TString tag = (region == "sideband") ? "_PreRecoilCut_" : "_";

  TString USDS = USDSIn;

  if (proj && loBin == 0 && hiBin == -1){
    if (USDS=="US"){
      loBin = 5;
      hiBin = 8;
    }
    else{
      loBin = 2;
      hiBin = 5;
    }
  }
  else if (!proj){
    loBin = 0;
    hiBin = -1;
  }

  double USDS_ratio = 1.0;

  if (USDSIn == "Sum") {
    MnvH1D* USBkg = (MnvH1D*)mcFile->Get("pTmu"+tag+mat+"_background_USPlastic");
    MnvH1D* DSBkg = (MnvH1D*)mcFile->Get("pTmu"+tag+mat+"_background_DSPlastic");

    USDS_ratio = USBkg->Integral()/DSBkg->Integral();

    USDS = "US";
    loBin = 5;
    hiBin = 8;
  }

  MnvH2D* m2D_data = (MnvH2D*)dataFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_data");
  MnvH2D* m2D_sig = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_selected_signal_reco");
  MnvH2D* m2D_1PiC = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_1chargePi");
  MnvH2D* m2D_1Pi0 = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_1neutPi");
  MnvH2D* m2D_NPi = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_NPi");
  MnvH2D* m2D_Other = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_Other");
  MnvH2D* m2D_Tgt = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_Wrong_Nucleus");
  MnvH1D* pT_data = m2D_data->ProjectionX("pT_data",loBin, hiBin);
  MnvH1D* pT_sig = m2D_sig->ProjectionX("pT_sig",loBin, hiBin);
  MnvH1D* pT_1PiC = m2D_1PiC->ProjectionX("pT_1PiC",loBin, hiBin);
  MnvH1D* pT_1Pi0 = m2D_1Pi0->ProjectionX("pT_1Pi0",loBin, hiBin);
  MnvH1D* pT_NPi = m2D_NPi->ProjectionX("pT_NPi",loBin, hiBin);
  MnvH1D* pT_Other = m2D_Other->ProjectionX("pT_Other",loBin, hiBin);
  MnvH1D* pT_Tgt = m2D_Tgt->ProjectionX("pT_Tgt",loBin, hiBin);

  if (USDSIn == "Sum"){
    USDS = "DS";
    loBin = 2;
    hiBin = 5;
    MnvH2D* m2D_data_2 = (MnvH2D*)dataFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_data");
    MnvH2D* m2D_sig_2 = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_selected_signal_reco");
    MnvH2D* m2D_1PiC_2 = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_1chargePi");
    MnvH2D* m2D_1Pi0_2 = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_1neutPi");
    MnvH2D* m2D_NPi_2 = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_NPi");
    MnvH2D* m2D_Other_2 = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_Other");
    MnvH2D* m2D_Tgt_2 = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_Wrong_Nucleus");
    pT_data->Add(m2D_data_2->ProjectionX("pT_data_2",loBin, hiBin),USDS_ratio);
    pT_sig->Add(m2D_sig_2->ProjectionX("pT_sig_2",loBin, hiBin),USDS_ratio);
    pT_1PiC->Add(m2D_1PiC_2->ProjectionX("pT_1PiC_2",loBin, hiBin),USDS_ratio);
    pT_1Pi0->Add(m2D_1Pi0_2->ProjectionX("pT_1Pi0_2",loBin, hiBin),USDS_ratio);
    pT_NPi->Add(m2D_NPi_2->ProjectionX("pT_NPi_2",loBin, hiBin),USDS_ratio);
    pT_Other->Add(m2D_Other_2->ProjectionX("pT_Other_2",loBin, hiBin),USDS_ratio);
    pT_Tgt->Add(m2D_Tgt_2->ProjectionX("pT_Tgt_2",loBin, hiBin),USDS_ratio);
    USDS = "Sum";
  }

  MnvH1D* pT_sum = (MnvH1D*)pT_sig->Clone();
  pT_sum->Add(pT_1PiC);
  pT_sum->Add(pT_1Pi0);
  pT_sum->Add(pT_NPi);
  pT_sum->Add(pT_Other);
  pT_sum->Add(pT_Tgt);

  pT_data->AddMissingErrorBandsAndFillWithCV(*pT_sig);

  cout << "Chi 2 for this distribution: " << Chi2(pT_sum,pT_data) << endl;

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();
  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();
  double areaScale = topArea/bottomArea;
  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);
  top->cd();
  pT_sig->SetLineColor(TColor::GetColor("#999933"));
  pT_sig->SetFillColor(TColor::GetColor("#999933"));
  pT_1PiC->SetFillColor(TColor::GetColor("#88CCEE"));
  pT_1PiC->SetLineColor(TColor::GetColor("#88CCEE"));
  pT_1Pi0->SetLineColor(TColor::GetColor("#117733"));
  pT_1Pi0->SetFillColor(TColor::GetColor("#117733"));
  pT_NPi->SetFillColor(TColor::GetColor("#CC6677"));
  pT_NPi->SetLineColor(TColor::GetColor("#CC6677"));
  pT_Other->SetLineColor(TColor::GetColor("#882255"));
  pT_Other->SetFillColor(TColor::GetColor("#882255"));
  pT_Tgt->SetLineColor(TColor::GetColor("#909497"));
  pT_Tgt->SetFillColor(TColor::GetColor("#909497"));

  TH1D* h_pT_data = (TH1D*)((pT_data->GetBinNormalizedCopy()).GetCVHistoWithError().Clone());
  THStack* h = new THStack();
  h_pT_data->SetLineColor(kBlack);
  h_pT_data->SetLineWidth(3);
  h->Add((TH1D*)((pT_Tgt->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->Add((TH1D*)((pT_Other->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->Add((TH1D*)((pT_NPi->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->Add((TH1D*)((pT_1Pi0->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->Add((TH1D*)((pT_1PiC->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->Add((TH1D*)((pT_sig->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->SetMaximum(h_pT_data->GetMaximum()*1.05);
  h->Draw("hist");
  h_pT_data->Draw("same");
  bottom->cd();
  gStyle->SetOptStat(0);
  MnvH1D* ratio = (MnvH1D*)pT_data->Clone();
  ratio->Divide(ratio,pT_sum);
  TH1D* mcRatio = new TH1D(pT_sum->GetTotalError(false, true, false));
    for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
      mcRatio->SetBinError(iBin, std::max(mcRatio->GetBinContent(iBin),1.0e-9));
      mcRatio->SetBinContent(iBin, 1);
    }
  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(3);
  ratio->SetTitle("");
  mcRatio->SetLineColor(kRed);
  mcRatio->SetLineWidth(3);
  mcRatio->SetFillColorAlpha(kPink+1, 0.4);
  TH1D* straightLine = (TH1D*)mcRatio->Clone();
  straightLine->SetFillStyle(0);
  ratio->SetMinimum(0.5);
  ratio->SetMaximum(1.5);
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetTitleSize(0.05*areaScale);
  ratio->GetYaxis()->SetTitleOffset(0.75/areaScale);
  ratio->GetYaxis()->SetLabelSize(ratio->GetYaxis()->GetLabelSize()*areaScale);
  ratio->GetXaxis()->SetLabelSize(ratio->GetXaxis()->GetLabelSize()*areaScale);
  ratio->GetXaxis()->SetTitleSize(0.04*areaScale);
  ratio->Draw();
  mcRatio->Draw("E2, same");
  straightLine->Draw("hist, same");
  ratio->Draw("same");
  top->cd();
  TLegend* leg = new TLegend(0.6,0.6,0.9,0.9);
  leg->AddEntry(h_pT_data,"DATA");
  leg->AddEntry(pT_sig,"Signal");
  leg->AddEntry(pT_1PiC,"single #pi^{#pm}");
  leg->AddEntry(pT_1Pi0,"single #pi^{0}");
  leg->AddEntry(pT_NPi,"N#pi");
  leg->AddEntry(pT_Other,"Other");
  leg->AddEntry(pT_Tgt,"In Tgt");
  leg->Draw();
  (proj) ? c1->Print("pTmu_"+USDSIn+"_"+region+"_region_"+mat+"_loBin_"+to_string(loBin)+"_hiBin_"+to_string(hiBin)+".pdf") : c1->Print("pTmu_"+USDSIn+"_"+region+"_region_"+mat+"_noProjection.pdf");
  (proj) ? c1->Print("pTmu_"+USDSIn+"_"+region+"_region_"+mat+"_loBin_"+to_string(loBin)+"_hiBin_"+to_string(hiBin)+".png") : c1->Print("pTmu_"+USDSIn+"_"+region+"_region_"+mat+"_noProjection.png");
  (proj) ? c1->Print("pTmu_"+USDSIn+"_"+region+"_region_"+mat+"_loBin_"+to_string(loBin)+"_hiBin_"+to_string(hiBin)+".C") : c1->Print("pTmu_"+USDSIn+"_"+region+"_region_"+mat+"_noProjection.C");
}

void PlotCHProjY(TFile* dataFile, TFile* mcFile, TString USDS, TString region, TString mat){
  if (region != "signal" && region != "sideband"){
    cout << "Not a valid region. Exiting" << endl;
    return;
  }

  if (mat != "Fe" && mat != "Pb" && mat != "C" && mat != "Water"){
    cout << "Not a valid material. Exiting" << endl;
    return;
  }

  if (USDS != "US" && USDS != "DS"){
    cout << "Not a valid plastic region. Exiting" << endl;
    return;
  }

  TString tag = (region == "sideband") ? "_PreRecoilCut_" : "_";

  MnvH2D* m2D_data = (MnvH2D*)dataFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_data");
  MnvH2D* m2D_sig = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_selected_signal_reco");
  MnvH2D* m2D_1PiC = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_1chargePi");
  MnvH2D* m2D_1Pi0 = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_1neutPi");
  MnvH2D* m2D_NPi = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_NPi");
  MnvH2D* m2D_Other = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_Other");
  MnvH2D* m2D_Tgt = (MnvH2D*)mcFile->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_Wrong_Nucleus");
  MnvH1D* pT_data = m2D_data->ProjectionY("pT_data");
  MnvH1D* pT_sig = m2D_sig->ProjectionY("pT_sig");
  MnvH1D* pT_1PiC = m2D_1PiC->ProjectionY("pT_1PiC");
  MnvH1D* pT_1Pi0 = m2D_1Pi0->ProjectionY("pT_1Pi0");
  MnvH1D* pT_NPi = m2D_NPi->ProjectionY("pT_NPi");
  MnvH1D* pT_Other = m2D_Other->ProjectionY("pT_Other");
  MnvH1D* pT_Tgt = m2D_Tgt->ProjectionY("pT_Tgt");

  MnvH1D* pT_sum = (MnvH1D*)pT_sig->Clone();
  pT_sum->Add(pT_1PiC);
  pT_sum->Add(pT_1Pi0);
  pT_sum->Add(pT_NPi);
  pT_sum->Add(pT_Other);
  pT_sum->Add(pT_Tgt);

  pT_data->AddMissingErrorBandsAndFillWithCV(*pT_sig);

  cout << "Chi 2 for this distribution: " << Chi2(pT_sum,pT_data) << endl;

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();
  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();
  double areaScale = topArea/bottomArea;
  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);
  top->cd();
  pT_sig->SetLineColor(TColor::GetColor("#999933"));
  pT_sig->SetFillColor(TColor::GetColor("#999933"));
  pT_1PiC->SetFillColor(TColor::GetColor("#88CCEE"));
  pT_1PiC->SetLineColor(TColor::GetColor("#88CCEE"));
  pT_1Pi0->SetLineColor(TColor::GetColor("#117733"));
  pT_1Pi0->SetFillColor(TColor::GetColor("#117733"));
  pT_NPi->SetFillColor(TColor::GetColor("#CC6677"));
  pT_NPi->SetLineColor(TColor::GetColor("#CC6677"));
  pT_Other->SetLineColor(TColor::GetColor("#882255"));
  pT_Other->SetFillColor(TColor::GetColor("#882255"));
  pT_Tgt->SetLineColor(TColor::GetColor("#909497"));
  pT_Tgt->SetFillColor(TColor::GetColor("#909497"));
  TH1D* h_pT_data = (TH1D*)((pT_data->GetBinNormalizedCopy()).GetCVHistoWithError().Clone());
  THStack* h = new THStack();
  h_pT_data->SetLineColor(kBlack);
  h_pT_data->SetLineWidth(3);
  h->Add((TH1D*)((pT_Tgt->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->Add((TH1D*)((pT_Other->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->Add((TH1D*)((pT_NPi->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->Add((TH1D*)((pT_1Pi0->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->Add((TH1D*)((pT_1PiC->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->Add((TH1D*)((pT_sig->GetBinNormalizedCopy()).GetCVHistoWithError().Clone()));
  h->SetMaximum(h_pT_data->GetMaximum()*1.05);
  h->Draw("hist");
  h_pT_data->Draw("same");
  bottom->cd();
  gStyle->SetOptStat(0);
  MnvH1D* ratio = (MnvH1D*)pT_data->Clone();
  ratio->Divide(ratio,pT_sum);
  TH1D* mcRatio = new TH1D(pT_sum->GetTotalError(false, true, false));
    for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
      mcRatio->SetBinError(iBin, std::max(mcRatio->GetBinContent(iBin),1.0e-9));
      mcRatio->SetBinContent(iBin, 1);
    }
  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(3);
  ratio->SetTitle("");
  mcRatio->SetLineColor(kRed);
  mcRatio->SetLineWidth(3);
  mcRatio->SetFillColorAlpha(kPink+1, 0.4);
  TH1D* straightLine = (TH1D*)mcRatio->Clone();
  straightLine->SetFillStyle(0);
  ratio->SetMinimum(0.5);
  ratio->SetMaximum(1.5);
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetTitleSize(0.05*areaScale);
  ratio->GetYaxis()->SetTitleOffset(0.75/areaScale);
  ratio->GetYaxis()->SetLabelSize(ratio->GetYaxis()->GetLabelSize()*areaScale);
  ratio->GetXaxis()->SetLabelSize(ratio->GetXaxis()->GetLabelSize()*areaScale);
  ratio->GetXaxis()->SetTitleSize(0.04*areaScale);
  ratio->Draw();
  mcRatio->Draw("E2, same");
  straightLine->Draw("hist, same");
  ratio->Draw("same");
  top->cd();

  TLegend* leg = nullptr;
  if (USDS == "DS"){
    leg = new TLegend(0.75,0.6,0.9,0.9);
  }
  else if (USDS == "US"){
    leg = new TLegend(0.1,0.6,0.35,0.9);
  }
  leg->AddEntry(h_pT_data,"DATA");
  leg->AddEntry(pT_sig,"Signal");
  leg->AddEntry(pT_1PiC,"single #pi^{#pm}");
  leg->AddEntry(pT_1Pi0,"single #pi^{0}");
  leg->AddEntry(pT_NPi,"N#pi");
  leg->AddEntry(pT_Other,"Other");
  leg->AddEntry(pT_Tgt,"In Tgt");
  leg->Draw();

  cout << "Chi 2 for this distribution: " << Chi2(pT_sum,pT_data) << endl;

  c1->Print("NPlanes_"+USDS+"_"+region+"_region_"+mat+".pdf");
  c1->Print("NPlanes_"+USDS+"_"+region+"_region_"+mat+".png");
  c1->Print("NPlanes_"+USDS+"_"+region+"_region_"+mat+".C");
}


void runAllPlotCHProj(TFile* dataFile, TFile* mcFile, TString mat){
  vector<TString> USDS_tags = {"US", "DS"};
  vector<TString> regions = {"signal", "sideband"};
  for (auto region : regions){
    cout << "Region: " << region << endl;
    cout << "USDS: Sum" << endl;
    PlotCHProjX(dataFile, mcFile, "Sum", region, mat);
    /*
    for (auto USDS : USDS_tags){
      cout << "USDS: " << USDS << endl;
      cout << "pTmu with projections." << endl;
      PlotCHProjX(dataFile, mcFile, USDS, region, mat);
      cout << "pTmu without chopped projections." << endl;
      PlotCHProjX(dataFile, mcFile, USDS, region, mat, false);
      cout << "NPlanes without chopped projections." << endl;
      PlotCHProjY(dataFile, mcFile, USDS, region, mat);      
      /*
      for (int i=1; i<=10; ++i){
	cout << "Bin by bin, bin: " << i << endl;
	PlotCHProjX(dataFile, mcFile, USDS, region, mat, true, i, i);
      }

    }      */
    cout << "" << endl;
  }
}
