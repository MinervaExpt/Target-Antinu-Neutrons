//File: CategoryPlots.h
//Definitions of useful plotting functions to be used by multiple plotting scripts
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

//Don't know what needs including... can go through at some point and test I guess...

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
#include "PlotUtils/MnvPlotter.h"

//My includes
#include "plotTools/plot.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

using namespace std;
using namespace PlotUtils;

void TESTDraw2DBKGCategLines(string name, TFile* mcFile, TFile* dataFile, TString sample, double scale, TString nameToSave){

  bool primPar = false;

  TString sampleName = sample;

  bool isTracker = sampleName.Contains("Tracker") ? true : false;

  MnvH2D* h_Sig_Top = (MnvH2D*)mcFile->Get((TString)name);
  //MnvH2D* h_Sig = new MnvH1D(h_Sig_Top->GetBinNormalizedCopy());
  TH2* h_Sig = (TH2*)h_Sig_Top->GetCVHistoWithError().Clone();
  h_Sig->Scale(scale);
  h_Sig->SetLineColor(TColor::GetColor("#999933"));
  h_Sig->SetFillColor(TColor::GetColor("#999933"));

  TH2* mcSum = (TH2*)h_Sig->Clone();
  mcSum->SetLineColor(kRed);
  mcSum->SetLineWidth(3);

  cout << "Handling: " << name << endl;
  string title = (string)h_Sig->GetTitle();
  TString Xtitle = h_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_Sig->GetYaxis()->GetTitle();
  string units = Xtitle.Data();
  units.erase(0,units.find("["));

  string name_bkg = name;
  name_bkg.erase(name_bkg.length()-21,name_bkg.length());

  MnvH2D* h_1PiC_Bkg_Top = (MnvH2D*)mcFile->Get((TString)name_bkg+"_by_BKG_Label_1chargePi");
  //MnvH2D* h_1PiC_Bkg = new MnvH1D(h_1PiC_Bkg_Top->GetBinNormalizedCopy());
  TH2* h_1PiC_Bkg = (TH2*)h_1PiC_Bkg_Top->GetCVHistoWithError().Clone();
  h_1PiC_Bkg->Scale(scale);
  mcSum->Add(h_1PiC_Bkg);
  h_1PiC_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_1PiC_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));

  MnvH2D* h_1Pi0_Bkg_Top = (MnvH2D*)mcFile->Get((TString)name_bkg+"_by_BKG_Label_1neutPi");
  //MnvH2D* h_1Pi0_Bkg = new MnvH1D(h_1Pi0_Bkg_Top->GetBinNormalizedCopy());
  TH2* h_1Pi0_Bkg = (TH2*)h_1Pi0_Bkg_Top->GetCVHistoWithError().Clone();
  h_1Pi0_Bkg->Scale(scale);
  mcSum->Add(h_1Pi0_Bkg);
  h_1Pi0_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_1Pi0_Bkg->SetFillColor(TColor::GetColor("#117733"));

  MnvH2D* h_NPi_Bkg_Top = (MnvH2D*)mcFile->Get((TString)name_bkg+"_by_BKG_Label_NPi");
  //MnvH2D* h_NPi_Bkg = new MnvH1D(h_NPi_Bkg_Top->GetBinNormalizedCopy());
  TH2* h_NPi_Bkg = (TH2*)h_NPi_Bkg_Top->GetCVHistoWithError().Clone();
  h_NPi_Bkg->Scale(scale);
  mcSum->Add(h_NPi_Bkg);
  h_NPi_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_NPi_Bkg->SetFillColor(TColor::GetColor("#CC6677"));

  /*
  MnvH2D* h_USPlastic_Bkg_Top = (MnvH2D*)mcFile->Get((TString)name_bkg+"_background_USPlastic");
  //MnvH2D* h_USPlastic_Bkg = new MnvH1D(h_USPlastic_Bkg_Top->GetBinNormalizedCopy());
  TH2* h_USPlastic_Bkg = (TH2*)h_USPlastic_Bkg_Top->GetCVHistoWithError().Clone();
  h_USPlastic_Bkg->Scale(scale);
  mcSum->Add(h_USPlastic_Bkg);
  h_USPlastic_Bkg->SetLineColor(TColor::GetColor("#A26E1C"));
  h_USPlastic_Bkg->SetFillColor(TColor::GetColor("#A26E1C"));

  MnvH2D* h_DSPlastic_Bkg_Top = (MnvH2D*)mcFile->Get((TString)name_bkg+"_background_DSPlastic");
  //MnvH2D* h_DSPlastic_Bkg = new MnvH1D(h_DSPlastic_Bkg_Top->GetBinNormalizedCopy());
  TH2* h_DSPlastic_Bkg = (TH2*)h_DSPlastic_Bkg_Top->GetCVHistoWithError().Clone();
  h_DSPlastic_Bkg->Scale(scale);
  mcSum->Add(h_DSPlastic_Bkg);
  h_DSPlastic_Bkg->SetLineColor(TColor::GetColor("#C1B185"));
  h_DSPlastic_Bkg->SetFillColor(TColor::GetColor("#C1B185"));

  MnvH2D* h_Wrong_Nucleus_Bkg_Top = (MnvH2D*)mcFile->Get((TString)name_bkg+"_background_Wrong_Nucleus");
  //MnvH2D* h_Wrong_Nucleus_Bkg = new MnvH1D(h_Wrong_Nucleus_Bkg_Top->GetBinNormalizedCopy());
  TH2* h_Wrong_Nucleus_Bkg = (TH2*)h_Wrong_Nucleus_Bkg_Top->GetCVHistoWithError().Clone();
  h_Wrong_Nucleus_Bkg->Scale(scale);
  mcSum->Add(h_Wrong_Nucleus_Bkg);
  h_Wrong_Nucleus_Bkg->SetLineColor(TColor::GetColor("#909497"));
  h_Wrong_Nucleus_Bkg->SetFillColor(TColor::GetColor("#909497"));
  */

  MnvH2D* h_Other_Bkg_Top = (MnvH2D*)mcFile->Get((TString)name_bkg+"_by_BKG_Label_Other");
  //MnvH2D* h_Other_Bkg = new MnvH1D(h_Other_Bkg_Top->GetBinNormalizedCopy());
  TH2* h_Other_Bkg = (TH2*)h_Other_Bkg_Top->GetCVHistoWithError().Clone();
  h_Other_Bkg->Scale(scale);
  mcSum->Add(h_Other_Bkg);
  h_Other_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#882255"));

  MnvH2D* h_data_Top = (MnvH2D*)dataFile->Get((TString)name_bkg+"_data");
  //MnvH2D* h_data = new MnvH1D(h_data_Top->GetBinNormalizedCopy());
  TH2* dataHist = (TH2*)h_data_Top->GetCVHistoWithError().Clone();
  dataHist->SetLineColor(kBlack);
  dataHist->SetLineWidth(3);

  vector<pair<TH2*, const char*>> hVec;
  //hVec.push_back(make_pair(h_Wrong_Nucleus_Bkg,"hists"));
  //hVec.push_back(make_pair(h_USPlastic_Bkg,"hists"));
  //hVec.push_back(make_pair(h_DSPlastic_Bkg,"hists"));
  hVec.push_back(make_pair(h_Other_Bkg,"hists"));
  hVec.push_back(make_pair(h_NPi_Bkg,"hists"));
  hVec.push_back(make_pair(h_1Pi0_Bkg,"hists"));
  hVec.push_back(make_pair(h_1PiC_Bkg,"hists"));
  hVec.push_back(make_pair(h_Sig,"hists"));
  hVec.push_back(make_pair(mcSum,"E2"));
  hVec.push_back(make_pair(dataHist,""));

  double multipliers[]={1,1,1,1,
			1,1,1,1,
			1,1,1,1};

  GridCanvas* gc=plotYAxis1D(hVec, "Eh","Eh-err", multipliers);
  gc->Print(nameToSave+"_BKG_stacked.pdf");
  gc->Print(nameToSave+"_BKG_stacked.png");

  /*
  THStack* h = new THStack();
  h->Add((TH2D*)h_Wrong_Nucleus_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH2D*)h_USPlastic_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH2D*)h_DSPlastic_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH2D*)h_Other_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH2D*)h_NPi_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH2D*)h_1Pi0_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH2D*)h_1PiC_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH2D*)h_Sig->GetCVHistoWithError().Clone());

  //TLegend* leg = new TLegend(1.0-0.6,1.0-0.5,1.0-0.9,1.0-0.9);
  TLegend* leg = new TLegend(0.1,0.5,0.4,0.9);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry(h_Sig,"Signal");
  leg->AddEntry(h_1PiC_Bkg,"single #pi^{#pm}");
  leg->AddEntry(h_1Pi0_Bkg,"single #pi^{0}");
  leg->AddEntry(h_NPi_Bkg,"N#pi");
  leg->AddEntry(h_Other_Bkg,"Other");
  if (!isTracker) leg->AddEntry(h_DSPlastic_Bkg,"DS Plastic");
  if (!isTracker) leg->AddEntry(h_USPlastic_Bkg,"US Plastic");
  if (!isTracker) leg->AddEntry(h_Wrong_Nucleus_Bkg,"Wrong Nucleus");

  leg->Draw();

  TH2D* mcRatio = new TH2D(mcSum->GetTotalError(false, true, false));
  */
  
  delete mcSum;
  delete dataHist;
  //delete straightLine;
  //delete h_data;
  delete h_Sig;
  delete h_1PiC_Bkg;
  delete h_1Pi0_Bkg;
  delete h_NPi_Bkg;
  //delete h_USPlastic_Bkg;
  //delete h_DSPlastic_Bkg;
  //delete h_Wrong_Nucleus_Bkg;
  delete h_Other_Bkg;
  delete gc;

  return;
}

void DrawBKGCategTEST(string name, TFile* mcFile, TFile* dataFile, TString sample, double scale, TString nameToSave){

  bool primPar = false;

  TString sampleName = sample;

  bool isTracker = sampleName.Contains("Tracker") ? true : false;

  MnvH1D* h_Sig_Top = (MnvH1D*)mcFile->Get((TString)name);
  MnvH1D* h_Sig = new MnvH1D(h_Sig_Top->GetBinNormalizedCopy());
  h_Sig->Scale(scale);
  MnvH1D* mcSum = (MnvH1D*)h_Sig->Clone();
  h_Sig->SetLineColor(TColor::GetColor("#999933"));
  h_Sig->SetFillColor(TColor::GetColor("#999933"));

  cout << "Handling: " << name << endl;
  string title = (string)h_Sig->GetTitle();
  TString Xtitle = h_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_Sig->GetYaxis()->GetTitle();
  string units = Xtitle.Data();
  units.erase(0,units.find("["));

  string name_bkg = name;
  name_bkg.erase(name_bkg.length()-21,name_bkg.length());

  MnvH1D* h_1PiC_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_1chargePi");
  MnvH1D* h_1PiC_Bkg = new MnvH1D(h_1PiC_Bkg_Top->GetBinNormalizedCopy());
  h_1PiC_Bkg->Scale(scale);
  mcSum->Add(h_1PiC_Bkg);
  h_1PiC_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_1PiC_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));

  MnvH1D* h_1Pi0_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_1neutPi");
  MnvH1D* h_1Pi0_Bkg = new MnvH1D(h_1Pi0_Bkg_Top->GetBinNormalizedCopy());
  h_1Pi0_Bkg->Scale(scale);
  mcSum->Add(h_1Pi0_Bkg);
  h_1Pi0_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_1Pi0_Bkg->SetFillColor(TColor::GetColor("#117733"));

  MnvH1D* h_NPi_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_NPi");
  MnvH1D* h_NPi_Bkg = new MnvH1D(h_NPi_Bkg_Top->GetBinNormalizedCopy());
  h_NPi_Bkg->Scale(scale);
  mcSum->Add(h_NPi_Bkg);
  h_NPi_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_NPi_Bkg->SetFillColor(TColor::GetColor("#CC6677"));

  MnvH1D* h_USPlastic_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_USPlastic");
  MnvH1D* h_USPlastic_Bkg = new MnvH1D(h_USPlastic_Bkg_Top->GetBinNormalizedCopy());
  h_USPlastic_Bkg->Scale(scale);
  mcSum->Add(h_USPlastic_Bkg);
  h_USPlastic_Bkg->SetLineColor(TColor::GetColor("#A26E1C"));
  h_USPlastic_Bkg->SetFillColor(TColor::GetColor("#A26E1C"));

  MnvH1D* h_DSPlastic_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_DSPlastic");
  MnvH1D* h_DSPlastic_Bkg = new MnvH1D(h_DSPlastic_Bkg_Top->GetBinNormalizedCopy());
  h_DSPlastic_Bkg->Scale(scale);
  mcSum->Add(h_DSPlastic_Bkg);
  h_DSPlastic_Bkg->SetLineColor(TColor::GetColor("#C1B185"));
  h_DSPlastic_Bkg->SetFillColor(TColor::GetColor("#C1B185"));

  MnvH1D* h_Wrong_Nucleus_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_Wrong_Nucleus");
  MnvH1D* h_Wrong_Nucleus_Bkg = new MnvH1D(h_Wrong_Nucleus_Bkg_Top->GetBinNormalizedCopy());
  h_Wrong_Nucleus_Bkg->Scale(scale);
  mcSum->Add(h_Wrong_Nucleus_Bkg);
  h_Wrong_Nucleus_Bkg->SetLineColor(TColor::GetColor("#909497"));
  h_Wrong_Nucleus_Bkg->SetFillColor(TColor::GetColor("#909497"));

  MnvH1D* h_Other_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_Other");
  MnvH1D* h_Other_Bkg = new MnvH1D(h_Other_Bkg_Top->GetBinNormalizedCopy());
  h_Other_Bkg->Scale(scale);
  mcSum->Add(h_Other_Bkg);
  h_Other_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_data_Top = (MnvH1D*)dataFile->Get((TString)name_bkg+"_data");
  MnvH1D* h_data = new MnvH1D(h_data_Top->GetBinNormalizedCopy());
  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone();
  dataHist->SetLineColor(kBlack);
  dataHist->SetLineWidth(3);
  h_data->AddMissingErrorBandsAndFillWithCV(*h_Sig);

  THStack* h = new THStack();
  h->Add((TH1D*)h_Wrong_Nucleus_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_USPlastic_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_DSPlastic_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Other_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_NPi_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_1Pi0_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_1PiC_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Sig->GetCVHistoWithError().Clone());

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  cout << "areaScale: " << areaScale << endl;

  h->Draw("hist");
  c1->Update();

  size_t pos = 0;
  if ((pos=name.find("pTmu_")) != string::npos){
    cout << "Fixing Axis?" << endl;
    h->GetXaxis()->SetRangeUser(0,2.5);
  }

  pos = 0;
  if ((pos=name.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"Other");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.04);
    primPar = true;
  }

  h->SetTitle(sampleName);//+" "+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.04);
  if (units.length() != 0) h->GetYaxis()->SetTitle("Events / "+(TString)units);
  else h->GetYaxis()->SetTitle("Events / bin");
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.75);
  if (!primPar)h->SetMaximum((dataHist->GetMaximum())*1.25);

  pos=0;
  if ((pos=name.find("_ENHitSB")) != string::npos){
    sampleName = "EM Blob E/NHit SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name.find("_NBlobsSB")) != string::npos){
    sampleName = "N EM Blob SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name.find("_MichelSB")) != string::npos){
    sampleName = "Michel SideBand " + sample;
    h->SetTitle(sampleName);
  }

  /*
  if (Xtitle.Contains("pmu")){
    h->GetXaxis()->SetTitle("p_{#mu} [GeV/c]");
    h->SetTitle("Muon Momentum "+sampleName);
  }

  if (Xtitle.Contains("vtxZ")){
    h->GetXaxis()->SetTitle("vtx. Z [mm]");
    h->SetTitle("Vertex Z "+sampleName);
  }
  */

  //h->Draw("hist");
  c1->Update();
  /* This is useful for debugging whether systematics actuall changed.
  TH1D* h_Tot = (TH1D*)h->GetStack()->Last()->Clone();
  h_Tot->SetLineColor(kRed);
  h_Tot->SetFillColorAlpha(kPink + 1, 0.4);
  h_Tot->Draw("E2 SAME");
  */
  dataHist->Draw("same");
  c1->Update();

  //TLegend* leg = new TLegend(1.0-0.6,1.0-0.5,1.0-0.9,1.0-0.9);
  TLegend* leg = new TLegend(0.1,0.5,0.4,0.9);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry(h_Sig,"Signal");
  leg->AddEntry(h_1PiC_Bkg,"single #pi^{#pm}");
  leg->AddEntry(h_1Pi0_Bkg,"single #pi^{0}");
  leg->AddEntry(h_NPi_Bkg,"N#pi");
  leg->AddEntry(h_Other_Bkg,"Other");
  if (!isTracker) leg->AddEntry(h_DSPlastic_Bkg,"DS Plastic");
  if (!isTracker) leg->AddEntry(h_USPlastic_Bkg,"US Plastic");
  if (!isTracker) leg->AddEntry(h_Wrong_Nucleus_Bkg,"Wrong Nucleus");

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  MnvH1D* ratio = (MnvH1D*)h_data->Clone();
  ratio->Divide(ratio,mcSum);

  TH1D* mcRatio = new TH1D(mcSum->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin, 1);
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

  pos=0;
  if ((pos=name.find("pTmu_")) != string::npos){
    cout << "Fixing Axis?" << endl;
    ratio->GetXaxis()->SetRangeUser(0,2.5);
  }
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

  c1->Print(nameToSave+"_BKG_stacked.pdf");
  c1->Print(nameToSave+"_BKG_stacked.png");
  top->SetLogy();
  c1->Update();
  c1->Print(nameToSave+"_BKG_stacked_log.pdf");
  c1->Print(nameToSave+"_BKG_stacked_log.png");         
  
  delete mcSum;
  delete dataHist;
  delete h;
  delete ratio;
  delete straightLine;
  delete h_data;
  delete h_Sig;
  delete h_1PiC_Bkg;
  delete h_1Pi0_Bkg;
  delete h_NPi_Bkg;
  delete h_USPlastic_Bkg;
  delete h_DSPlastic_Bkg;
  delete h_Wrong_Nucleus_Bkg;
  delete h_Other_Bkg;
  delete c1;

  return;
}

void DrawIntTypeTEST(string name_QE, TFile* mcFile, TFile* dataFile, TString sample, double scale, TString nameToSave){

  bool primPar = false;

  TString sampleName = sample;

  bool isTracker = sampleName.Contains("Tracker") ? true : false;

  MnvH1D* h_QE_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_QE);
  MnvH1D* h_QE_Sig = new MnvH1D(h_QE_Sig_Top->GetBinNormalizedCopy());
  h_QE_Sig->Scale(scale);
  MnvH1D* mcSum = (MnvH1D*)h_QE_Sig->Clone();
  h_QE_Sig->SetLineColor(TColor::GetColor("#88CCEE"));
  h_QE_Sig->SetFillColor(TColor::GetColor("#88CCEE"));

  //  string name_sig = (string)h_QE_Sig->GetName();
  string name_sig = name_QE;
  name_sig.erase(name_sig.length()-3,name_sig.length());
  string name_bkg = name_sig;
  name_bkg.erase(name_bkg.length()-12,name_bkg.length());

  cout << "Handling: " << name_sig << endl;
  string title = (string)h_QE_Sig->GetTitle();
  TString Xtitle = h_QE_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_QE_Sig->GetYaxis()->GetTitle();
  string units = Xtitle.Data();
  units.erase(0,units.find("["));
  /*
  cout << title << endl;
  title.erase(0, 8);
  cout << title << endl;
  cout << "" << endl;
  */

  MnvH1D* h_RES_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_RES");
  MnvH1D* h_RES_Sig = new MnvH1D(h_RES_Sig_Top->GetBinNormalizedCopy());
  h_RES_Sig->Scale(scale);
  mcSum->Add(h_RES_Sig);
  h_RES_Sig->SetLineColor(TColor::GetColor("#117733"));
  h_RES_Sig->SetFillColor(TColor::GetColor("#117733"));

  MnvH1D* h_DIS_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_DIS");
  MnvH1D* h_DIS_Sig = new MnvH1D(h_DIS_Sig_Top->GetBinNormalizedCopy());
  h_DIS_Sig->Scale(scale);
  mcSum->Add(h_DIS_Sig);
  h_DIS_Sig->SetLineColor(TColor::GetColor("#CC6677"));
  h_DIS_Sig->SetFillColor(TColor::GetColor("#CC6677"));

  MnvH1D* h_2p2h_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_2p2h");
  MnvH1D* h_2p2h_Sig = new MnvH1D(h_2p2h_Sig_Top->GetBinNormalizedCopy());
  h_2p2h_Sig->Scale(scale);
  mcSum->Add(h_2p2h_Sig);
  h_2p2h_Sig->SetLineColor(TColor::GetColor("#44AA99"));
  h_2p2h_Sig->SetFillColor(TColor::GetColor("#44AA99"));

  MnvH1D* h_Other_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_Other");
  MnvH1D* h_Other_Sig = new MnvH1D(h_Other_Sig_Top->GetBinNormalizedCopy());
  h_Other_Sig->Scale(scale);
  mcSum->Add(h_Other_Sig);
  h_Other_Sig->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_QE_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_QE"); 
  MnvH1D* h_QE_Bkg = new MnvH1D(h_QE_Bkg_Top->GetBinNormalizedCopy());
  h_QE_Bkg->Scale(scale);
  mcSum->Add(h_QE_Bkg);
  h_QE_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_QE_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_QE_Bkg->SetFillStyle(3003);

  MnvH1D* h_RES_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_RES");
  MnvH1D* h_RES_Bkg = new MnvH1D(h_RES_Bkg_Top->GetBinNormalizedCopy());
  h_RES_Bkg->Scale(scale);
  mcSum->Add(h_RES_Bkg);
  h_RES_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_RES_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_RES_Bkg->SetFillStyle(3003);

  MnvH1D* h_DIS_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_DIS");
  MnvH1D* h_DIS_Bkg = new MnvH1D(h_DIS_Bkg_Top->GetBinNormalizedCopy());
  h_DIS_Bkg->Scale(scale);
  mcSum->Add(h_DIS_Bkg);
  h_DIS_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_DIS_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_DIS_Bkg->SetFillStyle(3003);

  MnvH1D* h_2p2h_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_2p2h");
  MnvH1D* h_2p2h_Bkg = new MnvH1D(h_2p2h_Bkg_Top->GetBinNormalizedCopy());
  h_2p2h_Bkg->Scale(scale);
  mcSum->Add(h_2p2h_Bkg);
  h_2p2h_Bkg->SetLineColor(TColor::GetColor("#44AA99"));
  h_2p2h_Bkg->SetFillColor(TColor::GetColor("#44AA99"));
  h_2p2h_Bkg->SetFillStyle(3003);

  MnvH1D* h_USPlastic_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_USPlastic");
  MnvH1D* h_USPlastic_Bkg = new MnvH1D(h_USPlastic_Bkg_Top->GetBinNormalizedCopy());
  h_USPlastic_Bkg->Scale(scale);
  mcSum->Add(h_USPlastic_Bkg);
  h_USPlastic_Bkg->SetLineColor(TColor::GetColor("#A26E1C"));
  h_USPlastic_Bkg->SetFillColor(TColor::GetColor("#A26E1C"));
  h_USPlastic_Bkg->SetFillStyle(3003);

  MnvH1D* h_DSPlastic_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_DSPlastic");
  MnvH1D* h_DSPlastic_Bkg = new MnvH1D(h_DSPlastic_Bkg_Top->GetBinNormalizedCopy());
  h_DSPlastic_Bkg->Scale(scale);
  mcSum->Add(h_DSPlastic_Bkg);
  h_DSPlastic_Bkg->SetLineColor(TColor::GetColor("#C1B185"));
  h_DSPlastic_Bkg->SetFillColor(TColor::GetColor("#C1B185"));
  h_DSPlastic_Bkg->SetFillStyle(3003);

  MnvH1D* h_Wrong_Nucleus_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_Wrong_Nucleus");
  MnvH1D* h_Wrong_Nucleus_Bkg = new MnvH1D(h_Wrong_Nucleus_Bkg_Top->GetBinNormalizedCopy());
  h_Wrong_Nucleus_Bkg->Scale(scale);
  mcSum->Add(h_Wrong_Nucleus_Bkg);
  h_Wrong_Nucleus_Bkg->SetLineColor(TColor::GetColor("#909497"));
  h_Wrong_Nucleus_Bkg->SetFillColor(TColor::GetColor("#909497"));
  h_Wrong_Nucleus_Bkg->SetFillStyle(3003);

  MnvH1D* h_Other_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_Other");
  MnvH1D* h_Other_Bkg = new MnvH1D(h_Other_Bkg_Top->GetBinNormalizedCopy());
  h_Other_Bkg->Scale(scale);
  mcSum->Add(h_Other_Bkg);
  h_Other_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillStyle(3003);

  MnvH1D* h_data_Top = (MnvH1D*)dataFile->Get((TString)name_bkg+"_data");
  MnvH1D* h_data = new MnvH1D(h_data_Top->GetBinNormalizedCopy());
  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone();
  dataHist->SetLineColor(kBlack);
  dataHist->SetLineWidth(3);
  h_data->AddMissingErrorBandsAndFillWithCV(*h_QE_Sig);

  THStack* h = new THStack();
  h->Add((TH1D*)h_Wrong_Nucleus_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_USPlastic_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_DSPlastic_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Other_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_2p2h_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_DIS_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_RES_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_QE_Bkg->GetCVHistoWithError().Clone());

  h->Add((TH1D*)h_Other_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_2p2h_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_DIS_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_RES_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_QE_Sig->GetCVHistoWithError().Clone());

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  h->Draw("hist");
  c1->Update();

  size_t pos = 0;
  if ((pos=name_sig.find("pTmu_")) != string::npos){
    cout << "Fixing Axis?" << endl;
    h->GetXaxis()->SetRangeUser(0,2.5);
  }

  pos = 0;
  if ((pos=name_sig.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"Other");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.04);
    primPar = true;
  }

  h->SetTitle(sampleName);//+" "+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.04);
  if (units.length() != 0) h->GetYaxis()->SetTitle("Events / "+(TString)units);
  else h->GetYaxis()->SetTitle("Events / bin");
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.75);
  if (!primPar) h->SetMaximum((dataHist->GetMaximum())*1.25);
  
  pos=0;
  if ((pos=name_sig.find("_ENHitSB")) != string::npos){
    sampleName = "EM Blob E/NHit SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_NBlobsSB")) != string::npos){
    sampleName = "N EM Blob SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_MichelSB")) != string::npos){
    sampleName = "Michel SideBand " + sample;
    h->SetTitle(sampleName);
  }

  /*
  if (Xtitle.Contains("pmu")){
    h->GetXaxis()->SetTitle("p_{#mu} [GeV/c]");
    h->SetTitle("Muon Momentum "+sampleName);
  }

  if (Xtitle.Contains("vtxZ")){
    h->GetXaxis()->SetTitle("vtx. Z [mm]");
    h->SetTitle("Vertex Z "+sampleName);
  }
  */

  //h->Draw("hist");
  c1->Update();
  /* This is useful for debugging whether systematics actuall changed.
  TH1D* h_Tot = (TH1D*)h->GetStack()->Last()->Clone();
  h_Tot->SetLineColor(kRed);
  h_Tot->SetFillColorAlpha(kPink + 1, 0.4);
  h_Tot->Draw("E2 SAME");
  */
  dataHist->Draw("same");
  c1->Update();

  TH1D* hTmp = new TH1D();
  hTmp->SetFillColor(kWhite);
  hTmp->SetLineColor(kWhite);

  //egend* leg = new TLegend(0.6,0.5,0.9,0.9);
  TLegend* leg = new TLegend(0.1,0.5,0.4,0.9);

  leg->SetNColumns(2);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry((TObject*)0,"","");

  leg->AddEntry(h_QE_Sig,"Sig. + QE");
  leg->AddEntry(h_QE_Bkg,"Bkg. + QE");

  leg->AddEntry(h_RES_Sig,"Sig. + RES");
  leg->AddEntry(h_RES_Bkg,"Bkg. + RES");

  leg->AddEntry(h_DIS_Sig,"Sig. + DIS");
  leg->AddEntry(h_DIS_Bkg,"Bkg. + DIS");

  leg->AddEntry(h_2p2h_Sig,"Sig. + 2p2h");
  leg->AddEntry(h_2p2h_Bkg,"Bkg. + 2p2h");

  leg->AddEntry(h_Other_Sig,"Sig. + Other");
  leg->AddEntry(h_Other_Bkg,"Bkg. + Other");

  if (!isTracker){
    leg->AddEntry(hTmp,"","f");
    leg->AddEntry(h_DSPlastic_Bkg,"DS Plastic");

    leg->AddEntry(hTmp,"","f");
    leg->AddEntry(h_USPlastic_Bkg,"US Plastic");

    leg->AddEntry(hTmp,"","f");
    leg->AddEntry(h_Wrong_Nucleus_Bkg,"Wrong Nucleus");
  }

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  MnvH1D* ratio = (MnvH1D*)h_data->Clone();
  ratio->Divide(ratio,mcSum);

  TH1D* mcRatio = new TH1D(mcSum->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin, 1);
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

  pos=0;
  if ((pos=name_sig.find("pTmu_")) != string::npos){
    cout << "Fixing Axis?" << endl;
    ratio->GetXaxis()->SetRangeUser(0,2.5);
  }
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

  c1->Print(nameToSave+"_IntType_stacked.pdf");
  c1->Print(nameToSave+"_IntType_stacked.png");
  top->SetLogy();
  c1->Update();
  c1->Print(nameToSave+"_IntType_stacked_log.pdf");
  c1->Print(nameToSave+"_IntType_stacked_log.png");     

  delete mcSum;
  delete dataHist;
  delete h;
  delete ratio;
  delete straightLine;
  delete h_data;
  delete h_Wrong_Nucleus_Bkg;
  delete h_DSPlastic_Bkg;
  delete h_USPlastic_Bkg;
  delete h_Other_Bkg;
  delete h_2p2h_Bkg;
  delete h_DIS_Bkg;
  delete h_RES_Bkg;
  delete h_QE_Bkg;
  delete h_Other_Sig;
  delete h_2p2h_Sig;
  delete h_DIS_Sig;
  delete h_RES_Sig;
  delete h_QE_Sig;
  delete c1;

  delete hTmp;

  return;
}

void DrawTargetTypeTEST(string name_Plastic, TFile* mcFile, TFile* dataFile, TString sample, double scale, TString nameToSave){

  bool primPar = false;

  TString sampleName = sample;

  MnvH1D* h_Plastic_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_Plastic);
  MnvH1D* h_Plastic_Sig = new MnvH1D(h_Plastic_Sig_Top->GetBinNormalizedCopy());
  h_Plastic_Sig->Scale(scale);
  MnvH1D* mcSum = (MnvH1D*)h_Plastic_Sig->Clone();
  h_Plastic_Sig->SetLineColor(TColor::GetColor("#DDCC77"));
  h_Plastic_Sig->SetFillColor(TColor::GetColor("#DDCC77"));

  //string name_sig = (string)h_Plastic_Sig->GetName();
  string name_sig = name_Plastic;
  name_sig.erase(name_sig.length()-8,name_sig.length());
  string name_bkg = name_sig;
  name_bkg.erase(name_bkg.length()-15,name_bkg.length());

  cout << "Handling: " << name_sig << endl;
  string title = (string)h_Plastic_Sig->GetTitle();
  TString Xtitle = h_Plastic_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_Plastic_Sig->GetYaxis()->GetTitle();
  string units = Xtitle.Data();
  units.erase(0,units.find("["));
  /*
  cout << title << endl;
  title.erase(0, 7);
  cout << title << endl;
  cout << "" << endl;
  */

  MnvH1D* h_US_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_USPlastic");
  MnvH1D* h_US_Sig = new MnvH1D(h_US_Sig_Top->GetBinNormalizedCopy());
  h_US_Sig->Scale(scale);
  mcSum->Add(h_US_Sig);
  h_US_Sig->SetLineColor(TColor::GetColor("#A26E1C"));
  h_US_Sig->SetFillColor(TColor::GetColor("#A26E1C"));

  MnvH1D* h_DS_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_DSPlastic");
  MnvH1D* h_DS_Sig = new MnvH1D(h_DS_Sig_Top->GetBinNormalizedCopy());
  h_DS_Sig->Scale(scale);
  mcSum->Add(h_DS_Sig);
  h_DS_Sig->SetLineColor(TColor::GetColor("#C1B185"));
  h_DS_Sig->SetFillColor(TColor::GetColor("#C1B185"));

  MnvH1D* h_Fe_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_Fe");
  MnvH1D* h_Fe_Sig = new MnvH1D(h_Fe_Sig_Top->GetBinNormalizedCopy());
  h_Fe_Sig->Scale(scale);
  mcSum->Add(h_Fe_Sig);
  h_Fe_Sig->SetLineColor(TColor::GetColor("#882255"));
  h_Fe_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_Pb_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_Pb");
  MnvH1D* h_Pb_Sig = new MnvH1D(h_Pb_Sig_Top->GetBinNormalizedCopy());
  h_Pb_Sig->Scale(scale);
  mcSum->Add(h_Pb_Sig);
  h_Pb_Sig->SetLineColor(TColor::GetColor("#117733"));
  h_Pb_Sig->SetFillColor(TColor::GetColor("#117733"));

  MnvH1D* h_Water_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_Water");
  MnvH1D* h_Water_Sig = new MnvH1D(h_Water_Sig_Top->GetBinNormalizedCopy());
  h_Water_Sig->Scale(scale);
  mcSum->Add(h_Water_Sig);
  h_Water_Sig->SetLineColor(TColor::GetColor("#332288"));
  h_Water_Sig->SetFillColor(TColor::GetColor("#332288"));

  MnvH1D* h_C_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_C");
  MnvH1D* h_C_Sig = new MnvH1D(h_C_Sig_Top->GetBinNormalizedCopy());
  h_C_Sig->Scale(scale);
  mcSum->Add(h_C_Sig);
  h_C_Sig->SetLineColor(TColor::GetColor("#88CCEE"));
  h_C_Sig->SetFillColor(TColor::GetColor("#88CCEE"));

  //h_Prot_Sig->SetFillColor(TColor::GetColor("#999933"));

  MnvH1D* h_Other_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_Other");
  MnvH1D* h_Other_Sig = new MnvH1D(h_Other_Sig_Top->GetBinNormalizedCopy());
  h_Other_Sig->Scale(scale);
  mcSum->Add(h_Other_Sig);
  h_Other_Sig->SetLineColor(TColor::GetColor("#CC6677"));
  h_Other_Sig->SetFillColor(TColor::GetColor("#CC6677"));

  //h_None_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_Plastic_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_Plastic");
  MnvH1D* h_Plastic_Bkg = new MnvH1D(h_Plastic_Bkg_Top->GetBinNormalizedCopy());
  h_Plastic_Bkg->Scale(scale);
  mcSum->Add(h_Plastic_Bkg);
  h_Plastic_Bkg->SetLineColor(TColor::GetColor("#DDCC77"));
  h_Plastic_Bkg->SetFillColor(TColor::GetColor("#DDCC77"));
  h_Plastic_Bkg->SetFillStyle(3003);

  MnvH1D* h_US_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_USPlastic");
  MnvH1D* h_US_Bkg = new MnvH1D(h_US_Bkg_Top->GetBinNormalizedCopy());
  h_US_Bkg->Scale(scale);
  mcSum->Add(h_US_Bkg);
  h_US_Bkg->SetLineColor(TColor::GetColor("#A26E1C"));
  h_US_Bkg->SetFillColor(TColor::GetColor("#A26E1C"));
  h_US_Bkg->SetFillStyle(3003);

  MnvH1D* h_DS_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_DSPlastic");
  MnvH1D* h_DS_Bkg = new MnvH1D(h_DS_Bkg_Top->GetBinNormalizedCopy());
  h_DS_Bkg->Scale(scale);
  mcSum->Add(h_DS_Bkg);
  h_DS_Bkg->SetLineColor(TColor::GetColor("#C1B185"));
  h_DS_Bkg->SetFillColor(TColor::GetColor("#C1B185"));
  h_DS_Bkg->SetFillStyle(3003);

  MnvH1D* h_Fe_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_Fe");
  MnvH1D* h_Fe_Bkg = new MnvH1D(h_Fe_Bkg_Top->GetBinNormalizedCopy());
  h_Fe_Bkg->Scale(scale);
  mcSum->Add(h_Fe_Bkg);
  h_Fe_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Fe_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_Fe_Bkg->SetFillStyle(3003);

  MnvH1D* h_Pb_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_Pb");
  MnvH1D* h_Pb_Bkg = new MnvH1D(h_Pb_Bkg_Top->GetBinNormalizedCopy());
  h_Pb_Bkg->Scale(scale);
  mcSum->Add(h_Pb_Bkg);
  h_Pb_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_Pb_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_Pb_Bkg->SetFillStyle(3003);

  MnvH1D* h_Water_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_Water");
  MnvH1D* h_Water_Bkg = new MnvH1D(h_Water_Bkg_Top->GetBinNormalizedCopy());
  h_Water_Bkg->Scale(scale);
  mcSum->Add(h_Water_Bkg);
  h_Water_Bkg->SetLineColor(TColor::GetColor("#332288"));
  h_Water_Bkg->SetFillColor(TColor::GetColor("#332288"));
  h_Water_Bkg->SetFillStyle(3003);

  MnvH1D* h_C_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_C");
  MnvH1D* h_C_Bkg = new MnvH1D(h_C_Bkg_Top->GetBinNormalizedCopy());
  h_C_Bkg->Scale(scale);
  mcSum->Add(h_C_Bkg);
  h_C_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_C_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_C_Bkg->SetFillStyle(3003);

  MnvH1D* h_Other_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_Other");
  MnvH1D* h_Other_Bkg = new MnvH1D(h_Other_Bkg_Top->GetBinNormalizedCopy());
  h_Other_Bkg->Scale(scale);
  mcSum->Add(h_Other_Bkg);
  h_Other_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_Other_Bkg->SetFillStyle(3003);

  MnvH1D* h_data_Top = (MnvH1D*)dataFile->Get((TString)name_bkg+"_data");
  MnvH1D* h_data = new MnvH1D(h_data_Top->GetBinNormalizedCopy());
  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone(); 
  dataHist->SetLineColor(kBlack);
  dataHist->SetLineWidth(3);
  h_data->AddMissingErrorBandsAndFillWithCV(*h_Plastic_Sig);

  THStack* h = new THStack();
  h->Add((TH1D*)h_Other_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_C_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Water_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Pb_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Fe_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_DS_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_US_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Plastic_Bkg->GetCVHistoWithError().Clone());

  h->Add((TH1D*)h_Other_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_C_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Water_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Pb_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Fe_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_DS_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_US_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Plastic_Sig->GetCVHistoWithError().Clone());

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  h->Draw("hist");
  c1->Update();

  size_t pos = 0;
  if ((pos=name_sig.find("pTmu_")) != string::npos){
    cout << "Fixing Axis?" << endl;
    h->GetXaxis()->SetRangeUser(0,2.5);
  }

  pos=0;
  if ((pos=name_sig.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"None");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.04);
    primPar = true;
  }

  h->SetTitle(sampleName);//+" "+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.04);
  if (units.length() != 0) h->GetYaxis()->SetTitle("Events / "+(TString)units);
  else h->GetYaxis()->SetTitle("Events / bin");
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.75);
  if (!primPar)h->SetMaximum((dataHist->GetMaximum())*1.25);
  
  pos=0;
  if ((pos=name_sig.find("_ENHitSB")) != string::npos){
    sampleName = "EM Blob E/NHit SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_NBlobsSB")) != string::npos){
    sampleName = "N EM Blob SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_MichelSB")) != string::npos){
    sampleName = "Michel SideBand " + sample;
    h->SetTitle(sampleName);
  }

  /*
  if (Xtitle.Contains("pmu")){
    h->GetXaxis()->SetTitle("p_{#mu} [GeV/c]");
    h->SetTitle("Muon Momentum "+sampleName);
  }

  if (Xtitle.Contains("vtxZ")){
    h->GetXaxis()->SetTitle("vtx. Z [mm]");
    h->SetTitle("Vertex Z "+sampleName);
  }
  */

  //h->Draw("hist");
  c1->Update();
  /* This is useful for debugging whether systematics actuall changed.
  TH1D* h_Tot = (TH1D*)h->GetStack()->Last()->Clone();
  h_Tot->SetLineColor(kRed);
  h_Tot->SetFillColorAlpha(kPink + 1, 0.4);
  h_Tot->Draw("E2 SAME");
  */
  dataHist->Draw("same");
  c1->Update();

  //TLegend* leg = new TLegend(0.6,0.5,0.9,0.9);
  TLegend* leg = new TLegend(0.1,0.5,0.4,0.9);

  leg->SetNColumns(2);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry((TObject*)0,"","");

  leg->AddEntry(h_Plastic_Sig,"Sig. + Plastic");
  leg->AddEntry(h_Plastic_Bkg,"Bkg. + Plastic");

  leg->AddEntry(h_US_Sig,"Sig. + US Plastic");
  leg->AddEntry(h_US_Bkg,"Bkg. + US Plastic");

  leg->AddEntry(h_DS_Sig,"Sig. + DS Plastic");
  leg->AddEntry(h_DS_Bkg,"Bkg. + DS Plastic");

  leg->AddEntry(h_Fe_Sig,"Sig. + Fe");
  leg->AddEntry(h_Fe_Bkg,"Bkg. + Fe");

  leg->AddEntry(h_Pb_Sig,"Sig. + Pb");  
  leg->AddEntry(h_Pb_Bkg,"Bkg. + Pb");
  
  leg->AddEntry(h_Water_Sig,"Sig. + Water");
  leg->AddEntry(h_Water_Bkg,"Bkg. + Water");

  leg->AddEntry(h_C_Sig,"Sig. + C");
  leg->AddEntry(h_C_Bkg,"Bkg. + C");

  leg->AddEntry(h_Other_Sig,"Sig. + Other");
  leg->AddEntry(h_Other_Bkg,"Bkg. + Other");

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  MnvH1D* ratio = (MnvH1D*)h_data->Clone();
  ratio->Divide(ratio,mcSum);

  TH1D* mcRatio = new TH1D(mcSum->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin, 1);
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

  pos=0;
  if ((pos=name_sig.find("pTmu_")) != string::npos){
    cout << "Fixing Axis?" << endl;
    ratio->GetXaxis()->SetRangeUser(0,2.5);
  }
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

  c1->Print(nameToSave+"_TargetType_stacked.pdf");
  c1->Print(nameToSave+"_TargetType_stacked.png");
  top->SetLogy();
  c1->Update();
  c1->Print(nameToSave+"_TargetType_stacked_log.pdf");
  c1->Print(nameToSave+"_TargetType_stacked_log.png");  

  delete mcSum;
  delete dataHist;
  delete h;
  delete ratio;
  delete straightLine;
  delete h_data;
  delete h_Other_Bkg;
  delete h_C_Bkg;
  delete h_Water_Bkg;
  delete h_Pb_Bkg;
  delete h_Fe_Bkg;
  delete h_DS_Bkg;
  delete h_US_Bkg;
  delete h_Plastic_Bkg;
  delete h_Other_Sig;
  delete h_C_Sig;
  delete h_Water_Sig;
  delete h_Pb_Sig;
  delete h_Fe_Sig;
  delete h_DS_Sig;
  delete h_US_Sig;
  delete h_Plastic_Sig;
  delete c1;

  return;
}

void DrawLeadBlobTypeTEST(string name_Neut, TFile* mcFile, TFile* dataFile, TString sample, double scale, TString nameToSave){

  bool primPar = false;

  TString sampleName = sample;

  MnvH1D* h_Neut_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_Neut);
  MnvH1D* h_Neut_Sig = new MnvH1D(h_Neut_Sig_Top->GetBinNormalizedCopy());
  h_Neut_Sig->Scale(scale);
  MnvH1D* mcSum = (MnvH1D*)h_Neut_Sig->Clone();
  h_Neut_Sig->SetLineColor(TColor::GetColor("#88CCEE"));
  h_Neut_Sig->SetFillColor(TColor::GetColor("#88CCEE"));

  //string name_sig = (string)h_Neut_Sig->GetName();
  string name_sig = name_Neut;
  name_sig.erase(name_sig.length()-5,name_sig.length());
  string name_bkg = name_sig;
  name_bkg.erase(name_bkg.length()-17,name_bkg.length());
  cout << "Handling: " << name_sig << endl;

  string title = (string)h_Neut_Sig->GetTitle();
  TString Xtitle = h_Neut_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_Neut_Sig->GetYaxis()->GetTitle();
  string units = Xtitle.Data();
  units.erase(0,units.find("["));

  /*
  cout << title << endl;
  title.erase(0, 10);
  cout << title << endl;
  cout << "" << endl;
  */

  MnvH1D* h_Mu_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_mu");
  MnvH1D* h_Mu_Sig = new MnvH1D(h_Mu_Sig_Top->GetBinNormalizedCopy());
  h_Mu_Sig->Scale(scale);
  mcSum->Add(h_Mu_Sig);
  h_Mu_Sig->SetLineColor(TColor::GetColor("#44AA99"));
  h_Mu_Sig->SetFillColor(TColor::GetColor("#44AA99"));

  MnvH1D* h_Pi0_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_pi0");
  MnvH1D* h_Pi0_Sig = new MnvH1D(h_Pi0_Sig_Top->GetBinNormalizedCopy());
  h_Pi0_Sig->Scale(scale);
  mcSum->Add(h_Pi0_Sig);
  h_Pi0_Sig->SetLineColor(TColor::GetColor("#117733"));
  h_Pi0_Sig->SetFillColor(TColor::GetColor("#117733"));

  MnvH1D* h_PiM_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_pim");
  MnvH1D* h_PiM_Sig = new MnvH1D(h_PiM_Sig_Top->GetBinNormalizedCopy());
  h_PiM_Sig->Scale(scale);
  mcSum->Add(h_PiM_Sig);
  h_PiM_Sig->SetLineColor(TColor::GetColor("#332288"));
  h_PiM_Sig->SetFillColor(TColor::GetColor("#332288"));

  MnvH1D* h_PiP_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_pip");
  MnvH1D* h_PiP_Sig = new MnvH1D(h_PiP_Sig_Top->GetBinNormalizedCopy());
  h_PiP_Sig->Scale(scale);
  mcSum->Add(h_PiP_Sig);
  h_PiP_Sig->SetLineColor(TColor::GetColor("#DDCC77"));
  h_PiP_Sig->SetFillColor(TColor::GetColor("#DDCC77"));

  MnvH1D* h_Prot_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_prot");
  MnvH1D* h_Prot_Sig = new MnvH1D(h_Prot_Sig_Top->GetBinNormalizedCopy());
  h_Prot_Sig->Scale(scale);
  mcSum->Add(h_Prot_Sig);
  h_Prot_Sig->SetLineColor(TColor::GetColor("#999933"));
  h_Prot_Sig->SetFillColor(TColor::GetColor("#999933"));

  MnvH1D* h_Other_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_Other");
  MnvH1D* h_Other_Sig = new MnvH1D(h_Other_Sig_Top->GetBinNormalizedCopy());
  h_Other_Sig->Scale(scale);
  mcSum->Add(h_Other_Sig);
  h_Other_Sig->SetLineColor(TColor::GetColor("#CC6677"));
  h_Other_Sig->SetFillColor(TColor::GetColor("#CC6677"));

  MnvH1D* h_None_Sig_Top = (MnvH1D*)mcFile->Get((TString)name_sig+"_None");
  MnvH1D* h_None_Sig = new MnvH1D(h_None_Sig_Top->GetBinNormalizedCopy());
  h_None_Sig->Scale(scale);
  mcSum->Add(h_None_Sig);
  h_None_Sig->SetLineColor(TColor::GetColor("#882255"));
  h_None_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_Neut_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_neut");
  MnvH1D* h_Neut_Bkg = new MnvH1D(h_Neut_Bkg_Top->GetBinNormalizedCopy());
  h_Neut_Bkg->Scale(scale);
  mcSum->Add(h_Neut_Bkg);
  h_Neut_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_Neut_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_Neut_Bkg->SetFillStyle(3003);

  MnvH1D* h_Mu_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_mu");
  MnvH1D* h_Mu_Bkg = new MnvH1D(h_Mu_Bkg_Top->GetBinNormalizedCopy());
  h_Mu_Bkg->Scale(scale);
  mcSum->Add(h_Mu_Bkg);
  h_Mu_Bkg->SetLineColor(TColor::GetColor("#44AA99"));
  h_Mu_Bkg->SetFillColor(TColor::GetColor("#44AA99"));
  h_Mu_Bkg->SetFillStyle(3003);

  MnvH1D* h_Pi0_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_pi0");
  MnvH1D* h_Pi0_Bkg = new MnvH1D(h_Pi0_Bkg_Top->GetBinNormalizedCopy());
  h_Pi0_Bkg->Scale(scale);
  mcSum->Add(h_Pi0_Bkg);
  h_Pi0_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_Pi0_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_Pi0_Bkg->SetFillStyle(3003);

  MnvH1D* h_PiM_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_pim");
  MnvH1D* h_PiM_Bkg = new MnvH1D(h_PiM_Bkg_Top->GetBinNormalizedCopy());
  h_PiM_Bkg->Scale(scale);
  mcSum->Add(h_PiM_Bkg);
  h_PiM_Bkg->SetLineColor(TColor::GetColor("#332288"));
  h_PiM_Bkg->SetFillColor(TColor::GetColor("#332288"));
  h_PiM_Bkg->SetFillStyle(3003);

  MnvH1D* h_PiP_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_pip");
  MnvH1D* h_PiP_Bkg = new MnvH1D(h_PiP_Bkg_Top->GetBinNormalizedCopy());
  h_PiP_Bkg->Scale(scale);
  mcSum->Add(h_PiP_Bkg);
  h_PiP_Bkg->SetLineColor(TColor::GetColor("#DDCC77"));
  h_PiP_Bkg->SetFillColor(TColor::GetColor("#DDCC77"));
  h_PiP_Bkg->SetFillStyle(3003);

  MnvH1D* h_Prot_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_prot");
  MnvH1D* h_Prot_Bkg = new MnvH1D(h_Prot_Bkg_Top->GetBinNormalizedCopy());
  h_Prot_Bkg->Scale(scale);
  mcSum->Add(h_Prot_Bkg);
  h_Prot_Bkg->SetLineColor(TColor::GetColor("#999933"));
  h_Prot_Bkg->SetFillColor(TColor::GetColor("#999933"));
  h_Prot_Bkg->SetFillStyle(3003);

  MnvH1D* h_Other_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_Other");
  MnvH1D* h_Other_Bkg = new MnvH1D(h_Other_Bkg_Top->GetBinNormalizedCopy());
  h_Other_Bkg->Scale(scale);
  mcSum->Add(h_Other_Bkg);
  h_Other_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_Other_Bkg->SetFillStyle(3003);

  MnvH1D* h_None_Bkg_Top = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_None");
  MnvH1D* h_None_Bkg = new MnvH1D(h_None_Bkg_Top->GetBinNormalizedCopy());
  h_None_Bkg->Scale(scale);
  mcSum->Add(h_None_Bkg);
  h_None_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_None_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_None_Bkg->SetFillStyle(3003);

  MnvH1D* h_data_Top = (MnvH1D*)dataFile->Get((TString)name_bkg+"_data");
  MnvH1D* h_data = new MnvH1D(h_data_Top->GetBinNormalizedCopy());
  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone();
  dataHist->SetLineColor(kBlack);
  dataHist->SetLineWidth(3);
  h_data->AddMissingErrorBandsAndFillWithCV(*h_Neut_Sig);

  THStack* h = new THStack();
  h->Add((TH1D*)h_None_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Other_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Prot_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_PiP_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_PiM_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Pi0_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Mu_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Neut_Bkg->GetCVHistoWithError().Clone());

  h->Add((TH1D*)h_None_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Other_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Prot_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_PiP_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_PiM_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Pi0_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Mu_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Neut_Sig->GetCVHistoWithError().Clone());

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  h->Draw("hist");
  c1->Update();

  //ToDo: Get the naming of the axes fixed to be what it needs to be/make it easier to automate. This will need to be accompanied by a change to the Variable class to get the labels correct there. Current changes temporary in the interest of making plots for a talk on Oct. 14, 2021 in the exclusives meeting.

  size_t pos = 0;
  if ((pos=name_sig.find("pTmu_")) != string::npos){
    cout << "Fixing Axis?" << endl;
    h->GetXaxis()->SetRangeUser(0,2.5);
  }

  pos=0;
  if ((pos=name_sig.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"None");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.04);
    primPar = true;
  }

  h->SetTitle(sampleName);//+" "+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.04);
  if (units.length() != 0) h->GetYaxis()->SetTitle("Events / "+(TString)units);
  else h->GetYaxis()->SetTitle("Events / bin");
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.75);
  if (!primPar)h->SetMaximum((dataHist->GetMaximum())*1.25);
  
  pos=0;
  if ((pos=name_sig.find("_ENHitSB")) != string::npos){
    sampleName = "EM Blob E/NHit SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_NBlobsSB")) != string::npos){
    sampleName = "N EM Blob SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_MichelSB")) != string::npos){
    sampleName = "Michel SideBand " + sample;
    h->SetTitle(sampleName);
  }

  /*
  if (Xtitle.Contains("pmu")){
    h->GetXaxis()->SetTitle("p_{#mu} [GeV/c]");
    h->SetTitle("Muon Momentum "+sampleName);
  }

  if (Xtitle.Contains("vtxZ")){
    h->GetXaxis()->SetTitle("vtx. Z [mm]");
    h->SetTitle("Vertex Z "+sampleName);
  }
  */

  //  h->Draw("hist");
  c1->Update();
  /* This is useful for debugging whether systematics actuall changed.
  TH1D* h_Tot = (TH1D*)h->GetStack()->Last()->Clone();
  h_Tot->SetLineColor(kRed);
  h_Tot->SetFillColorAlpha(kPink + 1, 0.4);
  h_Tot->Draw("E2 SAME");
  */
  dataHist->Draw("same");
  c1->Update();

  //TLegend* leg = new TLegend(0.6,0.5,0.9,0.9);
  TLegend* leg = new TLegend(0.1,0.5,0.4,0.9);

  leg->SetNColumns(2);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry((TObject*)0,"","");

  leg->AddEntry(h_Neut_Sig,"Sig. + n");
  leg->AddEntry(h_Neut_Bkg,"Bkg. + n");

  leg->AddEntry(h_Mu_Sig,"Sig. + #mu");
  leg->AddEntry(h_Mu_Bkg,"Bkg. + #mu");

  leg->AddEntry(h_Pi0_Sig,"Sig. + #pi^{0}");  
  leg->AddEntry(h_Pi0_Bkg,"Bkg. + #pi^{0}");
  
  leg->AddEntry(h_PiM_Sig,"Sig. + #pi^{-}");
  leg->AddEntry(h_PiM_Bkg,"Bkg. + #pi^{-}");

  leg->AddEntry(h_PiP_Sig,"Sig. + #pi^{+}");
  leg->AddEntry(h_PiP_Bkg,"Bkg. + #pi^{+}");

  leg->AddEntry(h_Prot_Sig,"Sig. + p");
  leg->AddEntry(h_Prot_Bkg,"Bkg. + p");

  leg->AddEntry(h_Other_Sig,"Sig. + Other");
  leg->AddEntry(h_Other_Bkg,"Bkg. + Other");

  leg->AddEntry(h_None_Sig,"Sig. + None");
  leg->AddEntry(h_None_Bkg,"Bkg. + None");

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  MnvH1D* ratio = (MnvH1D*)h_data->Clone();
  ratio->Divide(ratio,mcSum);

  TH1D* mcRatio = new TH1D(mcSum->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin, 1);
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

  pos=0;
  if ((pos=name_sig.find("pTmu_")) != string::npos){
    cout << "Fixing Axis?" << endl;
    ratio->GetXaxis()->SetRangeUser(0,2.5);
  }
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

  c1->Print(nameToSave+"_LeadBlobType_stacked.pdf");
  c1->Print(nameToSave+"_LeadBlobType_stacked.png");
  top->SetLogy();
  c1->Update();
  c1->Print(nameToSave+"_LeadBlobType_stacked_log.pdf");
  c1->Print(nameToSave+"_LeadBlobType_stacked_log.png");

  delete mcSum;
  delete dataHist;
  delete h;
  delete ratio;
  delete straightLine;
  delete h_data;
  delete h_None_Bkg;
  delete h_Other_Bkg;
  delete h_Prot_Bkg;
  delete h_PiP_Bkg;
  delete h_PiM_Bkg;
  delete h_Pi0_Bkg;
  delete h_Mu_Bkg;
  delete h_Neut_Bkg;
  delete h_None_Sig;
  delete h_Other_Sig;
  delete h_Prot_Sig;
  delete h_PiP_Sig;
  delete h_PiM_Sig;
  delete h_Pi0_Sig;
  delete h_Mu_Sig;
  delete h_Neut_Sig;
  delete c1;

  return;
}
