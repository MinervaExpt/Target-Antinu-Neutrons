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

TCanvas* DrawStack(TString tag, TString Side, TFile* inFile, MnvH1D* h_data, double scale){

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();
  gStyle->SetOptStat(0);

  MnvH1D* h_sig = (MnvH1D*)inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Plastic/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Plastic_selected_signal_reco");
  h_sig->SetLineColor(TColor::GetColor("#999933"));
  h_sig->SetFillColor(TColor::GetColor("#999933"));
  h_sig->Scale(scale);

  MnvH1D* h_1PiC_Bkg = (MnvH1D*)inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Plastic/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Plastic_background_1chargePi");
  h_1PiC_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_1PiC_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_1PiC_Bkg->Scale(scale);

  MnvH1D* h_1Pi0_Bkg = (MnvH1D*)inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Plastic/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Plastic_background_1neutPi");
  h_1Pi0_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_1Pi0_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_1Pi0_Bkg->Scale(scale);

  MnvH1D* h_NPi_Bkg = (MnvH1D*)inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Plastic/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Plastic_background_NPi");
  h_NPi_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_NPi_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_NPi_Bkg->Scale(scale);

  MnvH1D* h_Other_Bkg = (MnvH1D*)(inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Plastic/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Plastic_background_Other"))->Clone();
  h_Other_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Plastic/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Plastic_background_Wrong_Nucleus"))->Clone());
  h_Other_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Plastic/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Plastic_background_DSPlastic"))->Clone());
  h_Other_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Plastic/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Plastic_background_USPlastic"))->Clone());
  h_Other_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_Other_Bkg->Scale(scale);

  MnvH1D* h_Mat_Bkg = (MnvH1D*)inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_"+tag+"/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_"+tag+"_data");
  h_Mat_Bkg->SetLineColor(TColor::GetColor("#909497"));
  h_Mat_Bkg->SetFillColor(TColor::GetColor("#909497"));
  h_Mat_Bkg->Scale(scale);

  MnvH1D* h_US_Bkg = (MnvH1D*)inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_USPlastic/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_USPlastic_data");
  h_US_Bkg->SetLineColor(TColor::GetColor("#A26E1C"));
  h_US_Bkg->SetFillColor(TColor::GetColor("#A26E1C"));
  h_US_Bkg->Scale(scale);

  MnvH1D* h_DS_Bkg = (MnvH1D*)inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_DSPlastic/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_DSPlastic_data");
  h_DS_Bkg->SetLineColor(TColor::GetColor("#C1B185"));
  h_DS_Bkg->SetFillColor(TColor::GetColor("#C1B185"));
  h_DS_Bkg->Scale(scale);

  TString otherMat = "Pb";
  if (tag == "Fe") otherMat = "Pb";
  else if (tag == "Pb" ) otherMat = "Fe";

  MnvH1D* h_All_Bkg = (MnvH1D*)(inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_"+otherMat+"/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_"+otherMat+"_data"))->Clone();
  h_All_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_C/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_C_data"))->Clone());
  h_All_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Water/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Water_data"))->Clone());
  h_All_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Other/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Other_data"))->Clone());

  h_All_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_All_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_All_Bkg->Scale(scale);
  h_All_Bkg->SetFillStyle(3003);

  THStack* h = new THStack();
  h->Add(h_All_Bkg);
  h->Add(h_Mat_Bkg);
  h->Add(h_DS_Bkg);
  h->Add(h_US_Bkg);
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
  leg->AddEntry(h_US_Bkg,"USPlastic?");
  leg->AddEntry(h_DS_Bkg,"DSPlastic?");
  leg->AddEntry(h_Mat_Bkg,tag);
  leg->AddEntry(h_All_Bkg,"Other Materials");

  leg->Draw();
  c1->Update();
  return c1;
}

void PlotVtxDemonstration(TString tag = "Fe", TString Side = "US") {

  TString fileName="runEventLoopMC_SkippedSyst_MnvTuneV1_FVregion_SingleTarget_Tgt1_wNeutCuts_neutKE_10.000000.root";
  TString dataName="runEventLoopData_SkippedSyst_MnvTuneV1_FVregion_SingleTarget_Tgt1_wNeutCuts_neutKE_10.000000.root";
  TString outDir="Tgt1_VtxDemonstration/";

  TFile* inFile = new TFile(fileName,"READ");
  TFile* dataFile = new TFile(dataName,"READ");
  MnvH1D* h_data = (MnvH1D*)dataFile->Get("ByTgt_Tgt1_"+tag+"/"+Side+"_ByType_Other/vtxZ_ByTgt_Tgt1_"+tag+"_PreRecoilCut_Other_data");

  TParameter<double>* dataPOT = (TParameter<double>*)dataFile->Get("POTUsed");
  TParameter<double>* mcPOT = (TParameter<double>*)inFile->Get("POTUsed");

  double scale = dataPOT->GetVal()/mcPOT->GetVal();

  cout << scale << endl;

  TCanvas* c1 = DrawStack(tag,Side,inFile,h_data,scale);
  c1->Print(outDir+"NPlanes_Tgt1_"+tag+"_"+Side+".pdf");
  c1->Print(outDir+"NPlanes_Tgt1_"+tag+"_"+Side+".png");
  delete c1;

  cout << "HEY YOU DID IT!!!" << endl;
  return;
}
