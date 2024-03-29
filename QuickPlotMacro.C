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

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"

using namespace std;
using namespace PlotUtils;

TCanvas* DrawStack(TString name, TFile* inFile, MnvH1D* h_data, double scale){

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();
  gStyle->SetOptStat(0);

  MnvH1D* h_sig = (MnvH1D*)inFile->Get(name+"_selected_signal_reco");
  h_sig->SetLineColor(TColor::GetColor("#999933"));
  h_sig->SetFillColor(TColor::GetColor("#999933"));
  h_sig->Scale(scale);

  MnvH1D* h_1PiC_Bkg = (MnvH1D*)inFile->Get(name+"_background_1chargePi");
  h_1PiC_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_1PiC_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_1PiC_Bkg->Scale(scale);

  MnvH1D* h_1Pi0_Bkg = (MnvH1D*)inFile->Get((TString)name+"_background_1neutPi");
  h_1Pi0_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_1Pi0_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_1Pi0_Bkg->Scale(scale);

  MnvH1D* h_NPi_Bkg = (MnvH1D*)inFile->Get((TString)name+"_background_NPi");
  h_NPi_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_NPi_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_NPi_Bkg->Scale(scale);

  MnvH1D* h_Other_Bkg = (MnvH1D*)inFile->Get((TString)name+"_background_Other");
  h_Other_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_Other_Bkg->Scale(scale);

  THStack* h = new THStack();
  h->Add(h_Other_Bkg);
  h->Add(h_NPi_Bkg);
  h->Add(h_1Pi0_Bkg);
  h->Add(h_1PiC_Bkg);

  h->Add(h_sig);
  h_data->Draw();
  h->Draw("hist, same");
  h_data->Draw("same");
  c1->Update();

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  leg->AddEntry(h_data, "DATA");
  leg->AddEntry(h_sig,"Signal");
  leg->AddEntry(h_1PiC_Bkg,"single #pi^{+/-}");
  leg->AddEntry(h_1Pi0_Bkg,"single #pi^{0}");
  leg->AddEntry(h_NPi_Bkg,"N#pi");
  leg->AddEntry(h_Other_Bkg,"Other");

  leg->Draw();
  c1->Update();
  return c1;
}

void QuickPlotMacro() {

  TString fileName="/minerva/data/users/dlast/TargetNeutronsAna/NewStructure/EMSideBandsTrackerMnvTune1_STDplaylist_6A/runEventLoopMC.root";
  TString dataName="/minerva/data/users/dlast/TargetNeutronsAna/NewStructure/EMSideBandsTrackerMnvTune1_STDplaylist_6A/runEventLoopData.root";
  TString outDir="/minerva/data/users/dlast/TargetNeutronsAna/NewStructure/EMSideBandsTrackerMnvTune1_STDplaylist_6A/testPlots/";

  TFile* inFile = new TFile(fileName,"READ");
  TFile* dataFile = new TFile(dataName,"READ");
  MnvH1D* h_data = (MnvH1D*)dataFile->Get("pmu_data");
  MnvH1D* h_data_ENHitSB = (MnvH1D*)dataFile->Get("pmu_ENHitSB_data");
  MnvH1D* h_data_NBlobsSB = (MnvH1D*)dataFile->Get("pmu_NBlobsSB_data");

  TParameter<double>* dataPOT = (TParameter<double>*)dataFile->Get("POTUsed");
  TParameter<double>* mcPOT = (TParameter<double>*)inFile->Get("POTUsed");

  double scale = dataPOT->GetVal()/mcPOT->GetVal();

  cout << scale << endl;

  TCanvas* c1 = DrawStack("pmu",inFile,h_data,scale);
  c1->Print(outDir+"pmu_signal_region.pdf");
  c1->Print(outDir+"pmu_signal_region.png");
  delete c1;
  c1 = DrawStack("pmu_ENHitSB",inFile,h_data_ENHitSB,scale);
  c1->Print(outDir+"pmu_ENHitSB.pdf");
  c1->Print(outDir+"pmu_ENHitSB.png");
  delete c1;
  c1 = DrawStack("pmu_NBlobsSB",inFile,h_data_NBlobsSB,scale);
  c1->Print(outDir+"pmu_NBlobsSB.pdf");
  c1->Print(outDir+"pmu_NBlobsSB.png");
  delete c1;

  cout << "HEY YOU DID IT!!!" << endl;
  return;
}
