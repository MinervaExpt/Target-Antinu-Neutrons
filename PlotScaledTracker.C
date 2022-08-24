//File: QuickPlotMacro.C
//Info: Quick macro to make initial pmu sideband/signal region plots.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

//C++ includes
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <bitset>
#include <time.h>
#include <sys/stat.h>

//ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"
#include "TParameter.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"

using namespace std;
using namespace PlotUtils;

map<int, TString> colors = {{0,"#117733"},{1,"#CC6677"},{2,"#88CCEE"}};

void DrawFromMnvH1Ds(MnvH1D* h_data, map<TString, MnvH1D*> hFit, map<TString, MnvH1D*> hUnfit, bool sigFit, TString nameToSave){

  if (hFit.size()==0){
    cout << "This script should not be used to plot histograms not involved with fitting at this time." << endl;
    return;
  }

  MnvH1D* mcSum;
  MnvH1D* unfitSum;
  MnvH1D* h_sig;

  bool allFit = ((sigFit && hUnfit.size()==0) || (!sigFit && hUnfit.size()==1));

  //cout << "Getting signal hist." << endl;                                                                                                                                                                        

  if (sigFit) h_sig = hFit["Signal"];
  else h_sig = hUnfit["Signal"];

  mcSum = h_sig->Clone();

  //cout << "Summing hists." << endl;                                                                                                                                                                              

  for (auto hist:hFit){
    if (hist.first != "Signal") mcSum->Add(hist.second);
  }

  int iHist = 0;
  for (auto hist:hUnfit){
    if (hist.first != "Signal"){
      mcSum->Add(hist.second);
      if (iHist!=0)unfitSum->Add(hist.second);
      else unfitSum = hist.second->Clone();
      iHist++;
    }
  }

  //cout << "Coloring hists." << endl;                                                                                                                                                                             
  h_sig->SetLineColor(TColor::GetColor("#999933"));
  h_sig->SetFillColor(TColor::GetColor("#999933"));
  if (!allFit){
    unfitSum->SetLineColor(TColor::GetColor("#882255"));
    unfitSum->SetFillColor(TColor::GetColor("#882255"));
  }

  iHist = 0;
  for (auto hist:hFit){
    if (hist.first != "Signal"){
      hist.second->SetLineColor(TColor::GetColor(colors[iHist%3]));
      hist.second->SetFillColor(TColor::GetColor(colors[iHist%3]));
      iHist++;
    }
  }
  //cout << "Stacking hists." << endl;                                                                                                                                                                             

  THStack* h = new THStack();
  if(!allFit) h->Add((TH1D*)unfitSum->GetCVHistoWithError().Clone());
  for (auto hist:hFit){
    if(hist.first != "Signal") h->Add((TH1D*)hist.second->GetCVHistoWithError().Clone());
  }
  h->Add((TH1D*)h_sig->GetCVHistoWithError().Clone());

  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone();
  dataHist->SetLineColor(kBlack);
  dataHist->SetLineWidth(3);

  //cout << "Drawing hists." << endl;                                                                                                                                                                              

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  //cout << "areaScale: " << areaScale << endl;                                                                                                                                                                    

  h->Draw("hist");
  h->SetMaximum((dataHist->GetMaximum())*1.05);
  c1->Update();
  dataHist->Draw("same");
  c1->Update();

  //cout << "Legend time." << endl;                                                                                                                                                                                
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry(h_sig,"Signal");
  for (auto hist = hFit.rbegin(); hist != hFit.rend(); ++hist){
    if (hist->first != "Signal") leg->AddEntry(hist->second,hist->first);
  }
  if (!allFit) leg->AddEntry(unfitSum,"Not fit BKGs");

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  //cout << "Ratio business." << endl;                                                                                                                                                                             
  MnvH1D* ratio = (MnvH1D*)h_data->Clone();
  ratio->Divide(ratio,mcSum);

  TH1D* mcRatio = new TH1D(mcSum->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin,1);
  }

  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(3);
  ratio->SetTitle("");
  //ratio->SetTitleSize(0);                                                                                                                                                                                        
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetTitleSize(0.05*areaScale);
  ratio->GetYaxis()->SetTitleOffset(0.75/areaScale);
  ratio->GetYaxis()->SetLabelSize(ratio->GetYaxis()->GetLabelSize()*areaScale);

  ratio->GetXaxis()->SetLabelSize(ratio->GetXaxis()->GetLabelSize()*areaScale);
  ratio->GetXaxis()->SetTitleSize(0.04*areaScale);
  ratio->SetMinimum(0.5);
  ratio->SetMaximum(1.5);
  ratio->Draw();

  mcRatio->SetLineColor(kRed);
  mcRatio->SetLineWidth(3);
  mcRatio->SetFillColorAlpha(kPink + 1, 0.4);
  mcRatio->Draw("E2 SAME");

  TH1D* straightLine = (TH1D*)mcRatio->Clone();
  straightLine->SetFillStyle(0);
  straightLine->Draw("HIST SAME");

  ratio->Draw("SAME");

  c1->Update();
  c1->Print(nameToSave+".pdf");
  c1->Print(nameToSave+".png");
  top->SetLogy();
  c1->Update();
  c1->Print(nameToSave+"_log.pdf");
  c1->Print(nameToSave+"_log.png");

  delete c1;

  //Does deleting these fix my seg fault issues?                                                                                                                                                                   
  delete mcSum;
  if (!allFit) delete unfitSum;
  delete h;
  delete ratio;
  delete mcRatio;
  delete straightLine;

  return;
}

void PlotScaledTracker(TString fileName, TString dataName, TString scaleName, TString outDir, TString tag){

  gStyle->SetOptStat(0);

  TFile* mcFile = new TFile(fileName,"READ");
  TFile* dataFile = new TFile(dataName,"READ");
  TFile* scaleFile = new TFile(scaleName,"READ");

  TParameter<double>* dataPOT = (TParameter<double>*)dataFile->Get("POTUsed");
  TParameter<double>* mcPOT = (TParameter<double>*)mcFile->Get("POTUsed");

  double scale = dataPOT->GetVal()/mcPOT->GetVal();

  MnvH1D* dataHist_Top = (MnvH1D*)(dataFile->Get("pTmu"+tag+"_data"));
  MnvH1D* dataHist = new MnvH1D(dataHist_Top->GetBinNormalizedCopy());

  MnvH1D* sigHist_Top = (MnvH1D*)(mcFile->Get("pTmu"+tag+"_selected_signal_reco"))->Clone();
  sigHist_Top->Scale(scale);
  MnvH1D* sigHist = new MnvH1D(sigHist_Top->GetBinNormalizedCopy());

  MnvH1D* totHist_Top = (MnvH1D*)(mcFile->Get("pTmu"+tag+"_data"))->Clone();
  totHist_Top->Scale(scale);
  MnvH1D* totHist = new MnvH1D(totHist_Top->GetBinNormalizedCopy());

  MnvH1D* PiCHist_Top = (MnvH1D*)(mcFile->Get("pTmu"+tag+"_background_1chargePi"))->Clone();
  PiCHist_Top->Scale(scale);
  MnvH1D* PiCHist = new MnvH1D(PiCHist_Top->GetBinNormalizedCopy());

  MnvH1D* Pi0Hist_Top = (MnvH1D*)(mcFile->Get("pTmu"+tag+"_background_1neutPi"))->Clone();
  Pi0Hist_Top->Scale(scale);
  MnvH1D* Pi0Hist = new MnvH1D(Pi0Hist_Top->GetBinNormalizedCopy());

  MnvH1D* OnePiHist = (MnvH1D*)PiCHist->Clone();
  OnePiHist->Add(Pi0Hist);

  MnvH1D* NPiHist_Top = (MnvH1D*)(mcFile->Get("pTmu"+tag+"_background_NPi"))->Clone();
  NPiHist_Top->Scale(scale);
  MnvH1D* NPiHist = new MnvH1D(NPiHist_Top->GetBinNormalizedCopy());

  MnvH1D* RESHist_Top = (MnvH1D*)(mcFile->Get("pTmu"+tag+"_bkg_IntType_RES"))->Clone();
  RESHist_Top->Scale(scale);
  MnvH1D* RESHist = new MnvH1D(RESHist_Top->GetBinNormalizedCopy());

  MnvH1D* DISHist_Top = (MnvH1D*)(mcFile->Get("pTmu"+tag+"_bkg_IntType_DIS"))->Clone();
  DISHist_Top->Scale(scale);
  MnvH1D* DISHist = new MnvH1D(DISHist_Top->GetBinNormalizedCopy());

  MnvH1D* bkgHist = totHist->Clone();
  bkgHist->Add(sigHist,-1.0);

  MnvH1D* bkgHist_Other3B = bkgHist->Clone();
  bkgHist_Other3B->Add(OnePiHist,-1.0);
  bkgHist_Other3B->Add(NPiHist,-1.0);

  MnvH1D* bkgHist_Other6B = bkgHist->Clone();
  bkgHist_Other6B->Add(RESHist,-1.0);
  bkgHist_Other6B->Add(DISHist,-1.0);

  MnvH1D* bkgHist_scaled = bkgHist->Clone();
  MnvH1D* scaleHist1B_BKG = (MnvH1D*)(scaleFile->Get("pTmu_PreRecoilCut_fit_recoilE_PreRecoilCut_fit1B_low_26_hi_50_BKG"))->Clone();
  bkgHist_scaled->Multiply(bkgHist_scaled,scaleHist1B_BKG);

  MnvH1D* OnePiHist_scaled = OnePiHist->Clone();
  MnvH1D* scaleHist3B_1Pi = (MnvH1D*)(scaleFile->Get("pTmu_PreRecoilCut_fit_recoilE_PreRecoilCut_fit3B_low_26_hi_50_single\ #pi"))->Clone();
  OnePiHist_scaled->Multiply(OnePiHist_scaled,scaleHist3B_1Pi);
  MnvH1D* NPiHist_scaled = NPiHist->Clone();
  MnvH1D* scaleHist3B_NPi = (MnvH1D*)(scaleFile->Get("pTmu_PreRecoilCut_fit_recoilE_PreRecoilCut_fit3B_low_26_hi_50_N#pi"))->Clone();
  NPiHist_scaled->Multiply(NPiHist_scaled,scaleHist3B_NPi);

  MnvH1D* RESHist_scaled = RESHist->Clone();
  MnvH1D* scaleHist6B_RES = (MnvH1D*)(scaleFile->Get("pTmu_PreRecoilCut_fit_recoilE_PreRecoilCut_fit6B_low_26_hi_50_RES"))->Clone();
  RESHist_scaled->Multiply(RESHist_scaled,scaleHist6B_RES);
  MnvH1D* DISHist_scaled = DISHist->Clone();
  MnvH1D* scaleHist6B_DIS = (MnvH1D*)(scaleFile->Get("pTmu_PreRecoilCut_fit_recoilE_PreRecoilCut_fit6B_low_26_hi_50_DIS"))->Clone();
  DISHist_scaled->Multiply(DISHist_scaled,scaleHist6B_DIS);

  map<TString,MnvH1D*> fitHists1B, unFitHists1B, fitHists1B_scaled;
  map<TString,MnvH1D*> fitHists3B, unFitHists3B, fitHists3B_scaled;
  map<TString,MnvH1D*> fitHists6B, unFitHists6B, fitHists6B_scaled;

  unFitHists1B["Signal"] = sigHist;
  fitHists1B["BKG"] = bkgHist;
  fitHists1B_scaled["BKG"] = bkgHist_scaled;

  DrawFromMnvH1Ds(dataHist,fitHists1B,unFitHists1B,false,outDir+"pTmu"+tag+"_preFit_recoil_1B_low_26_hi_50");
  DrawFromMnvH1Ds(dataHist,fitHists1B_scaled,unFitHists1B,false,outDir+"pTmu"+tag+"_postFit_recoil_1B_low_26_hi_50");

  unFitHists3B["Signal"] = sigHist;
  unFitHists3B["Unfit BKGs"] = bkgHist_Other3B;
  fitHists3B["single #pi"] = OnePiHist;
  fitHists3B_scaled["single #pi"] = OnePiHist_scaled;
  fitHists3B["N#pi"] = NPiHist;
  fitHists3B_scaled["N#pi"] = NPiHist_scaled;

  DrawFromMnvH1Ds(dataHist,fitHists3B,unFitHists3B,false,outDir+"pTmu"+tag+"_preFit_recoil_3B_low_26_hi_50");
  DrawFromMnvH1Ds(dataHist,fitHists3B_scaled,unFitHists3B,false,outDir+"pTmu"+tag+"_postFit_recoil_3B_low_26_hi_50");

  unFitHists6B["Signal"] = sigHist;
  unFitHists6B["Unfit BKGs"] = bkgHist_Other6B;
  fitHists6B["RES"] = RESHist;
  fitHists6B_scaled["RES"] = RESHist_scaled;
  fitHists6B["DIS"] = DISHist;
  fitHists6B_scaled["DIS"] = DISHist_scaled;

  DrawFromMnvH1Ds(dataHist,fitHists6B,unFitHists6B,false,outDir+"pTmu"+tag+"_preFit_recoil_6B_low_26_hi_50");
  DrawFromMnvH1Ds(dataHist,fitHists6B_scaled,unFitHists6B,false,outDir+"pTmu"+tag+"_postFit_recoil_6B_low_26_hi_50");

  cout << "Deleting" << endl;

  fitHists1B.clear();
  unFitHists1B.clear();
  fitHists1B_scaled.clear();

  fitHists3B.clear();
  unFitHists3B.clear();
  fitHists3B_scaled.clear();

  fitHists6B.clear();
  unFitHists6B.clear();
  fitHists6B_scaled.clear();

  delete dataHist;
  delete dataHist_Top;
  delete sigHist;
  delete sigHist_Top;
  delete totHist;
  delete totHist_Top;
  delete PiCHist;
  delete PiCHist_Top;
  delete Pi0Hist;
  delete Pi0Hist_Top;
  delete bkgHist_Other3B;
  delete bkgHist_Other6B;
  delete bkgHist;
  delete bkgHist_scaled;
  delete OnePiHist;
  delete OnePiHist_scaled;
  delete NPiHist;
  delete NPiHist_Top;
  delete NPiHist_scaled;
  delete RESHist;
  delete RESHist_Top;
  delete RESHist_scaled;
  delete DISHist;
  delete DISHist_Top;
  delete DISHist_scaled;
  delete scaleHist1B_BKG;
  delete scaleHist3B_1Pi;
  delete scaleHist3B_NPi;
  delete scaleHist6B_RES;
  delete scaleHist6B_DIS;

  cout << "Closing Files" << endl;

  mcFile->Close();
  dataFile->Close();
  scaleFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return;
}
