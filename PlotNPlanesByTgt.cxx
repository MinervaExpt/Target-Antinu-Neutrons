//File: PlotNPlanesByTgt.cxx
//Info: Plots the number of planes from the target for interstitial plastic both upstream and downstream showing the tuning categories for the plastic.
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

map<int, vector<TString>> mats = {{1,{"Fe","Pb"}},
                                  {2,{"Fe","Pb"}},
                                  {3,{"C","Fe","Pb"}},
                                  {4,{"Pb"}},
                                  {5,{"Fe","Pb"}},
                                  {6,{""}}};

vector<TString> TgtTypes = {"C","Fe","Pb","Water","USPlastic","DSPlastic","Other"};

void DrawUSDSVtx(int tgtNum, TString tag, TFile* inFile, TFile* dataFile, double scale, TString outdir){

  TString mat="_"+tag;

  TString tgtName = "Tgt"+to_string(tgtNum);
  if (tgtNum==6){
    tgtName = "WaterTgt";
    mat = "";
  }

  TString nameToSave = outdir+"/NPlanes_USDS_"+tgtName+mat;

  cout << "sig" << endl;

  MnvH1D* h_sig = (MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_selected_signal_reco"))->Clone();
  cout << "sig" << endl;
  h_sig->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_selected_signal_reco")));
  h_sig->SetLineColor(TColor::GetColor("#999933"));
  h_sig->SetFillColor(TColor::GetColor("#999933"));
  h_sig->Scale(scale);
  MnvH1D* mcSum = (MnvH1D*)h_sig->Clone();

  cout << "1PiC" << endl;
  MnvH1D* h_1PiC_Bkg = (MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_1chargePi"))->Clone();
  cout << "1PiC" << endl;
  h_1PiC_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_1chargePi")));
  h_1PiC_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_1PiC_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_1PiC_Bkg->Scale(scale);
  mcSum->Add(h_1PiC_Bkg);

  cout << "1Pi0" << endl;
  MnvH1D* h_1Pi0_Bkg = (MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_1neutPi"))->Clone();
  cout << "1Pi0" << endl;
  h_1Pi0_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_1neutPi")));
  h_1Pi0_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_1Pi0_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_1Pi0_Bkg->Scale(scale);
  mcSum->Add(h_1Pi0_Bkg);

  cout << "NPi" << endl;
  MnvH1D* h_NPi_Bkg = (MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_NPi"))->Clone();
  cout << "NPi" << endl;
  h_NPi_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_NPi")));
  h_NPi_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_NPi_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_NPi_Bkg->Scale(scale);
  mcSum->Add(h_NPi_Bkg);

  cout << "Other" << endl;
  MnvH1D* h_Other_Bkg = (MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_Other"))->Clone();
  cout << "Other" << endl;
  h_Other_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_Other")));
  cout << "Other" << endl;
  h_Other_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_Wrong_Nucleus")));
  cout << "Other" << endl;
  h_Other_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_Wrong_Nucleus")));
  cout << "Other" << endl;
  h_Other_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_DSPlastic")));
  cout << "Other" << endl;
  h_Other_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_DSPlastic")));
  cout << "Other" << endl;
  h_Other_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_USPlastic")));
  cout << "Other" << endl;
  h_Other_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_Plastic/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Plastic_background_USPlastic")));
  h_Other_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_Other_Bkg->Scale(scale);
  mcSum->Add(h_Other_Bkg);

  cout << "Mat" << endl;
  MnvH1D* h_Mat_Bkg = (tgtName != "WaterTgt") ? (MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType"+mat+"/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut"+mat+"_data"))->Clone() : (MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_Water/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Water_data"))->Clone();
  cout << "Mat" << endl;
  if (tgtName != "WaterTgt") h_Mat_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType"+mat+"/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut"+mat+"_data")));
  else h_Mat_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_Water/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Water_data")));
  h_Mat_Bkg->SetLineColor(TColor::GetColor("#909497"));
  h_Mat_Bkg->SetFillColor(TColor::GetColor("#909497"));
  h_Mat_Bkg->Scale(scale);
  mcSum->Add(h_Mat_Bkg);

  MnvH1D* h_All_Bkg = new MnvH1D("h_All_Bkg","",mcSum->GetNbinsX(),mcSum->GetXaxis()->GetXbins()->GetArray());
  h_All_Bkg->AddMissingErrorBandsAndFillWithCV(*mcSum);
  //LOOP HERE For All Materials, skip the one we're interested in.
  for (auto typeName : TgtTypes){
    TString compare = "_"+typeName;
    if (compare == mat || (tgtName == "WaterTgt" && typeName == "Water")) continue;
    cout << "typeName" << endl;
    cout << "All" << endl;
    h_All_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_"+typeName+"/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_"+typeName+"_data")));
    cout << "typeName" << endl;
    cout << "All" << endl;
    h_All_Bkg->Add((MnvH1D*)(inFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_"+typeName+"/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_"+typeName+"_data")));
  }
  h_All_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_All_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_All_Bkg->Scale(scale);
  mcSum->Add(h_All_Bkg);

  cout << "data" << endl;  
  MnvH1D* h_data = (MnvH1D*)(dataFile->Get("ByTgt_"+tgtName+mat+"/US_ByType_Other/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Other_data"))->Clone();
  cout << "data" << endl;
  h_data->Add((MnvH1D*)(dataFile->Get("ByTgt_"+tgtName+mat+"/DS_ByType_Other/vtxZ_ByTgt_"+tgtName+mat+"_PreRecoilCut_Other_data")));
  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone();
  h_data->AddMissingErrorBandsAndFillWithCV(*h_sig);

  THStack* h = new THStack();
  h->Add((TH1D*)h_All_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Mat_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Other_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_NPi_Bkg->GetCVHistoWithError().Clone());  
  h->Add((TH1D*)h_1Pi0_Bkg->GetCVHistoWithError().Clone());  
  h->Add((TH1D*)h_1PiC_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_sig->GetCVHistoWithError().Clone());

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

  h->SetMaximum((dataHist->GetMaximum())*1.25);

  dataHist->Draw("same");
  c1->Update();

  TLegend* leg = new TLegend(0.15,0.7,0.35,0.9);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry(h_sig,"Plastic Signal");
  leg->AddEntry(h_1PiC_Bkg,"Plastic #pi^{#pm}");
  leg->AddEntry(h_1Pi0_Bkg,"Plastic single #pi^{0}");
  leg->AddEntry(h_NPi_Bkg,"Plastic N#pi");
  leg->AddEntry(h_Other_Bkg,"Plastic Other BKG");
  if (tgtNum!=6) leg->AddEntry(h_Mat_Bkg,tag);
  else leg->AddEntry(h_Mat_Bkg,"Water");
  leg->AddEntry(h_All_Bkg,"Other Materials");

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
  ratio->Draw();

  mcRatio->SetLineColor(kRed);
  mcRatio->SetLineWidth(3);
  mcRatio->SetFillColorAlpha(kPink + 1, 0.4);
  mcRatio->Draw("E2 SAME");

  TH1D* straightLine = (TH1D*)mcRatio->Clone();
  straightLine->SetFillStyle(0);
  straightLine->Draw("HIST SAME");

  c1->Update();

  c1->Print(nameToSave+".pdf");
  c1->Print(nameToSave+".png");
  top->SetLogy();
  c1->Update();
  c1->Print(nameToSave+"_log.pdf");
  c1->Print(nameToSave+"_log.png");         
  
  delete mcSum;
  delete dataHist;
  delete h;
  delete ratio;
  delete straightLine;
  delete h_data;
  delete h_sig;
  delete h_1PiC_Bkg;
  delete h_1Pi0_Bkg;
  delete h_NPi_Bkg;
  delete h_Mat_Bkg;
  delete h_All_Bkg;
  delete h_Other_Bkg;
  delete c1;

  return;
}

bool PathExists(string path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

int main(int argc, char* argv[]) {

  gStyle->SetOptStat(0);

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc != 5) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string MCfileName=string(argv[1]);
  string DATAfileName=string(argv[2]);
  string outDir=string(argv[3]);
  int tgtNum = atoi(argv[4]);

  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory doesn't exist. Exiting" << endl;
    return 3;
  }

  string rootExt = ".root";
  string slash = "/";
  string token;
  string fileNameStub = MCfileName;
  size_t pos=0;

  //cout << sigNameStub << endl;
  while ((pos = fileNameStub.find(slash)) != string::npos){
    //cout << sigNameStub << endl;
    token = fileNameStub.substr(0, pos);
    //cout << token << endl;
    fileNameStub.erase(0, pos+slash.length());
  }
  //cout << sigNameStub << endl;
  if ((pos=fileNameStub.find(rootExt)) == string::npos){
    cout << "MC Input need be .root file." << endl;
    return 4;
  }

  cout << "Input MC file name parsed to: " << fileNameStub << endl;

  rootExt = ".root";
  slash = "/";
  token = "";
  fileNameStub = DATAfileName;
  pos=0;

  //cout << sigNameStub << endl;
  while ((pos = fileNameStub.find(slash)) != string::npos){
    //cout << sigNameStub << endl;
    token = fileNameStub.substr(0, pos);
    //cout << token << endl;
    fileNameStub.erase(0, pos+slash.length());
  }
  //cout << sigNameStub << endl;
  if ((pos=fileNameStub.find(rootExt)) == string::npos){
    cout << "DATA Input need be .root file." << endl;
    return 5;
  }

  cout << "Input Data file name parsed to: " << fileNameStub << endl;

  //cout << "Setting up MnvPlotter" << endl;
  //MnvPlotter* plotter = new MnvPlotter(kCCQEAntiNuStyle);

  TFile* mcFile = new TFile(MCfileName.c_str(),"READ");
  TFile* dataFile = new TFile(DATAfileName.c_str(),"READ");

  double mcPOT = ((TParameter<double>*)mcFile->Get("POTUsed"))->GetVal();
  double dataPOT = ((TParameter<double>*)dataFile->Get("POTUsed"))->GetVal();

  double scale = dataPOT/mcPOT;
  cout << "POT scale factor: " << scale << endl;
  
  for (auto mat : mats[tgtNum]){
    cout << "Drawing for " << tgtNum << ", " << mat << endl;
    DrawUSDSVtx (tgtNum, mat, mcFile, dataFile, scale, outDir);
  }

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  mcFile->Close();
  dataFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
