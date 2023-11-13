using namespace PlotUtils;
using namespace std;

void PlotCHProjX(TFile* dataFile, TFile* mcFile, TString USDS, TString region, TString mat, bool proj=true, int loBin=0, int hiBin=-1){
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
  else{
    loBin = 0;
    hiBin = -1;
  }

  MnvH2D* m2D_data = (MnvH2D*)_file1->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_data");
  MnvH2D* m2D_sig = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_selected_signal_reco");
  MnvH2D* m2D_1PiC = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_1chargePi");
  MnvH2D* m2D_1Pi0 = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_1neutPi");
  MnvH2D* m2D_NPi = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_NPi");
  MnvH2D* m2D_Other = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_Other");
  MnvH2D* m2D_Tgt = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_Wrong_Nucleus");
  MnvH1D* pT_data = m2D_data->ProjectionX("pT_data",loBin, hiBin);
  MnvH1D* pT_sig = m2D_sig->ProjectionX("pT_sig",loBin, hiBin);
  MnvH1D* pT_1PiC = m2D_1PiC->ProjectionX("pT_1PiC",loBin, hiBin);
  MnvH1D* pT_1Pi0 = m2D_1Pi0->ProjectionX("pT_1Pi0",loBin, hiBin);
  MnvH1D* pT_NPi = m2D_NPi->ProjectionX("pT_NPi",loBin, hiBin);
  MnvH1D* pT_Other = m2D_Other->ProjectionX("pT_Other",loBin, hiBin);
  MnvH1D* pT_Tgt = m2D_Tgt->ProjectionX("pT_Tgt",loBin, hiBin);
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
  MnvH1D* pT_sum = (MnvH1D*)pT_sig->Clone();
  pT_sum->Add(pT_1PiC);
  pT_sum->Add(pT_1Pi0);
  pT_sum->Add(pT_NPi);
  pT_sum->Add(pT_Other);
  pT_sum->Add(pT_Tgt);
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
  pT_data->AddMissingErrorBandsAndFillWithCV(*pT_sig);
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
  (proj) ? c1->Print("pTmu_"+USDS+"_"+region+"_region_"+mat+"_loBin_"+to_string(loBin)+"_hiBin_"+to_string(hiBin)+".pdf") : c1->Print("pTmu_"+USDS+"_"+region+"_region_"+mat+"_noProjection.pdf");
  (proj) ? c1->Print("pTmu_"+USDS+"_"+region+"_region_"+mat+"_loBin_"+to_string(loBin)+"_hiBin_"+to_string(hiBin)+".png") : c1->Print("pTmu_"+USDS+"_"+region+"_region_"+mat+"_noProjection.png");
  (proj) ? c1->Print("pTmu_"+USDS+"_"+region+"_region_"+mat+"_loBin_"+to_string(loBin)+"_hiBin_"+to_string(hiBin)+".C") : c1->Print("pTmu_"+USDS+"_"+region+"_region_"+mat+"_noProjection.C");
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

  MnvH2D* m2D_data = (MnvH2D*)_file1->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_data");
  MnvH2D* m2D_sig = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_selected_signal_reco");
  MnvH2D* m2D_1PiC = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_1chargePi");
  MnvH2D* m2D_1Pi0 = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_1neutPi");
  MnvH2D* m2D_NPi = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_NPi");
  MnvH2D* m2D_Other = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_Other");
  MnvH2D* m2D_Tgt = (MnvH2D*)_file0->Get("TwoD/TwoD_vtxZ_v_pT_Outer"+USDS+"Plastic"+tag+mat+"_background_Wrong_Nucleus");
  MnvH1D* pT_data = m2D_data->ProjectionY("pT_data");
  MnvH1D* pT_sig = m2D_sig->ProjectionY("pT_sig");
  MnvH1D* pT_1PiC = m2D_1PiC->ProjectionY("pT_1PiC");
  MnvH1D* pT_1Pi0 = m2D_1Pi0->ProjectionY("pT_1Pi0");
  MnvH1D* pT_NPi = m2D_NPi->ProjectionY("pT_NPi");
  MnvH1D* pT_Other = m2D_Other->ProjectionY("pT_Other");
  MnvH1D* pT_Tgt = m2D_Tgt->ProjectionY("pT_Tgt");
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
  MnvH1D* pT_sum = (MnvH1D*)pT_sig->Clone();
  pT_sum->Add(pT_1PiC);
  pT_sum->Add(pT_1Pi0);
  pT_sum->Add(pT_NPi);
  pT_sum->Add(pT_Other);
  pT_sum->Add(pT_Tgt);
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
  pT_data->AddMissingErrorBandsAndFillWithCV(*pT_sig);
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
  c1->Print("NPlanes_"+USDS+"_"+region+"_region_"+mat+".pdf");
  c1->Print("NPlanes_"+USDS+"_"+region+"_region_"+mat+".png");
  c1->Print("NPlanes_"+USDS+"_"+region+"_region_"+mat+".C");
}
