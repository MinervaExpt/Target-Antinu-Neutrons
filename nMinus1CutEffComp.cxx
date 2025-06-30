//File: nMinus1CutEffComp
//Info: Compares Specific cut efficiencies in data versus total MC with cuts N compared to N-1
//
//Usage: signalBKGStack <mc file N> <mc file N-1> <data file N> <data file N-1> <outputDir> <plot title>
//Author: David Last david.last@rochester.edu/lastd44@gmail.com

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
#include "TLatex.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

using namespace std;
using namespace PlotUtils;

double Chi2(MnvH1D* m1, MnvH1D* m2){
  if (m1->GetNbinsX() != m2->GetNbinsX()){
    cout << "Histos to compare don't have the same number of bins! Returning -999!!!" << endl;
    return -999.0;
  }
  if (m1->GetEntries() == 0 || m2->GetEntries() == 0){
    cout << "No entries in one of the input histograms. Assume this is faulty." << endl;
    return -99.0;
  }
  TH1D* h1 = (TH1D*)m1->GetCVHistoWithError().Clone();
  TH1D* h2 = (TH1D*)m2->GetCVHistoWithError().Clone();
  double chi2 = 0.0;
  for (int whichBin = 1; whichBin <= h1->GetNbinsX(); ++whichBin){
    double h1Content = h1->GetBinContent(whichBin);
    double h2Content = h2->GetBinContent(whichBin);
    double h1Err = h1->GetBinError(whichBin);
    double h2Err = h2->GetBinError(whichBin);
    double err = sqrt(h1Err*h1Err+h2Err*h2Err);
    double diff = h1Content-h2Content;
    //cout << "Content 1: " << h1Content << ", Content 2: " << h2Content <<", Error 1: " << h1Err << ", Error 2: " << h2Err << ", Error: " << err << ", Contribution: " << max(1e-10, fabs(diff*diff/(err*err))) << endl;
    if(h2Err > 1e-10) chi2 += (diff*diff)/(err*err);
  }

  delete h1;
  delete h2;
  return chi2;
}

void DrawEffComp(TString name, TFile* mcFile_N, TFile* mcFile_NM1, TFile* dataFile_N, TFile* dataFile_NM1, TString title, TString nameToSave, bool useSig=false){
  cout << "Handling: " << name << endl;

  mcFile_N->cd();

  MnvH1D* MC_N = nullptr;
  MnvH1D* MC_BKG_N = nullptr;
  
  if (useSig){
    string nameSTR = string(name.Data());
    nameSTR.erase(nameSTR.length()-5);
    TString tmpName = (TString)(nameSTR.c_str());
    MC_N = (MnvH1D*)(mcFile_N->Get(tmpName+"_selected_signal_reco")->Clone());
    MC_BKG_N = (MnvH1D*)(mcFile_N->Get(tmpName+"_background_1chargePi")->Clone());
    MC_BKG_N->Add((MnvH1D*)(mcFile_N->Get(tmpName+"_background_1neutPi")));
    MC_BKG_N->Add((MnvH1D*)(mcFile_N->Get(tmpName+"_background_NPi")));
    MC_BKG_N->Add((MnvH1D*)(mcFile_N->Get(tmpName+"_background_Other")));
    MC_BKG_N->Add((MnvH1D*)(mcFile_N->Get(tmpName+"_background_USPlastic")));
    MC_BKG_N->Add((MnvH1D*)(mcFile_N->Get(tmpName+"_background_DSPlastic")));
    MC_BKG_N->Add((MnvH1D*)(mcFile_N->Get(tmpName+"_background_Wrong_Nucleus")));
  }
  else{
    MC_N = new MnvH1D(*(MnvH1D*)(mcFile_N->Get(name)));
  }

  if (MC_N->GetEntries()==0){
    cout << "MC_N empty. Skipping." << endl;
    delete MC_N;
    delete MC_BKG_N;
    return;
  }

  double MCN_Int = MC_N->Integral(0,-1);

  mcFile_NM1->cd();

  MnvH1D* MC_NM1 = nullptr;
  MnvH1D* MC_BKG_NM1 = nullptr;
  
  if (useSig){
    string nameSTR = string(name.Data());
    nameSTR.erase(nameSTR.length()-5);
    TString tmpName = (TString)(nameSTR.c_str());
    MC_NM1 = (MnvH1D*)(mcFile_NM1->Get(tmpName+"_selected_signal_reco")->Clone());
    MC_BKG_NM1 = (MnvH1D*)(mcFile_NM1->Get(tmpName+"_background_1chargePi")->Clone());
    MC_BKG_NM1->Add((MnvH1D*)(mcFile_NM1->Get(tmpName+"_background_1neutPi")));
    MC_BKG_NM1->Add((MnvH1D*)(mcFile_NM1->Get(tmpName+"_background_NPi")));
    MC_BKG_NM1->Add((MnvH1D*)(mcFile_NM1->Get(tmpName+"_background_Other")));
    MC_BKG_NM1->Add((MnvH1D*)(mcFile_NM1->Get(tmpName+"_background_USPlastic")));
    MC_BKG_NM1->Add((MnvH1D*)(mcFile_NM1->Get(tmpName+"_background_DSPlastic")));
    MC_BKG_NM1->Add((MnvH1D*)(mcFile_NM1->Get(tmpName+"_background_Wrong_Nucleus")));
  }
  else{
    MC_NM1 = new MnvH1D(*(MnvH1D*)(mcFile_NM1->Get(name)));
  }


  if (MC_NM1->GetEntries()==0){
    cout << "MC_NM1 empty. Skipping." << endl;
    delete MC_N;
    delete MC_BKG_N;
    delete MC_NM1;
    delete MC_BKG_NM1;
    return;
  }

  double MCNM1_Int = MC_NM1->Integral(0,-1);
  
  MnvH1D* MC_eff = (MnvH1D*)(MC_N->Clone());
  MC_eff->Divide(MC_eff,MC_NM1);
  
  //Don't want flux error bands for this comparison just in case
  MC_eff->PopVertErrorBand("Flux");

  dataFile_N->cd();
  MnvH1D* data_N = new MnvH1D(*(MnvH1D*)(dataFile_N->Get(name)));
  if (data_N->GetEntries()==0){
    cout << "data_N empty. Skipping." << endl;
    delete MC_N;
    delete MC_BKG_N;
    delete MC_NM1;
    delete MC_BKG_NM1;
    delete MC_eff;
    delete data_N;
    return;
  }
  data_N->AddMissingErrorBandsAndFillWithCV(*MC_N);
  if (MC_BKG_N)data_N->Add(MC_BKG_N,-1.0);

  double dataN_Int = data_N->Integral(0,-1);

  dataFile_NM1->cd();
  MnvH1D* data_NM1 = new MnvH1D(*(MnvH1D*)(dataFile_NM1->Get(name)));
  if (data_NM1->GetEntries()==0){
    cout << "data_NM1 empty. Skipping." << endl;
    delete MC_N;
    delete MC_BKG_N;
    delete MC_NM1;
    delete MC_BKG_NM1;
    delete MC_eff;
    delete data_N;
    delete data_NM1;
    return;
  }
  data_NM1->AddMissingErrorBandsAndFillWithCV(*MC_NM1);
  if (MC_BKG_NM1)data_NM1->Add(MC_BKG_NM1,-1.0);
  
  double dataNM1_Int = data_NM1->Integral(0,-1);
  
  MnvH1D* data_eff = (MnvH1D*)(data_N->Clone());
  data_eff->Divide(data_eff,data_NM1);

  cout << "MC_N: " << MCN_Int << ", MC_NM1: " << MCNM1_Int << ", data_N: " << dataN_Int << ", data_NM1: " << dataNM1_Int << endl;
  
  double chi2 = Chi2(MC_eff,data_eff);
  
  cout << "Chi2: " << chi2 << endl;

  if (chi2 == -99 || chi2 == -999){
    cout << "Bad Chi2. Skipping" << endl;
    delete MC_N;
    delete MC_NM1;
    delete MC_eff;
    delete data_N;
    delete data_NM1;
    delete data_eff;
    return;
  }
  
  TH1D* mcHist = new TH1D(MC_eff->GetCVHistoWithError());
  mcHist->SetLineColor(kRed);
  TH1D* errHist = (TH1D*)mcHist->Clone();
  errHist->SetFillColorAlpha(kPink + 1, 0.4);
  
  TH1D* dataHist = new TH1D(data_eff->GetCVHistoWithError());
  dataHist->SetLineColor(kBlack);
  dataHist->SetLineWidth(3);

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  mcHist->SetMaximum((dataHist->GetMaximum())*1.25);

  mcHist->SetTitle(title);
  
  mcHist->Draw("hist");
  errHist->Draw("E2 SAME");
  c1->Update();

  dataHist->Draw("same");
  c1->Update();

  TLatex* latex = new TLatex( 0.2, 0.83, "MINER#nuA Work in Progress" );
  latex->SetTextFont(43);
  latex->SetTextSize(32);
  latex->SetNDC();
  latex->Draw();
  latex->SetTextColor(kRed);
  c1->Update();

  TLegend* leg = new TLegend(0.7,0.9,0.7,0.9);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry(mcHist,"MC");

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  MnvH1D* ratio = (MnvH1D*)data_eff->Clone();
  ratio->Divide(ratio,MC_eff);
  TH1D* ratioHist = new TH1D(ratio->GetCVHistoWithError());
  TString Xtitle = mcHist->GetXaxis()->GetTitle();
  ratioHist->GetXaxis()->SetTitle(Xtitle);

  TH1D* mcRatio = new TH1D(MC_eff->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin, 1);
  }

  ratioHist->SetLineColor(kBlack);
  ratioHist->SetLineWidth(3);
  ratioHist->SetTitle("");
  ratioHist->GetYaxis()->SetTitle("Data / MC");
  ratioHist->GetYaxis()->SetTitleSize(0.05*areaScale);
  ratioHist->GetYaxis()->SetTitleOffset(0.75/areaScale);
  ratioHist->GetYaxis()->SetLabelSize(ratioHist->GetYaxis()->GetLabelSize()*areaScale);

  ratioHist->GetXaxis()->SetLabelSize(ratioHist->GetXaxis()->GetLabelSize()*areaScale);
  ratioHist->GetXaxis()->SetTitleSize(0.04*areaScale);
  ratioHist->SetMinimum(0.5);
  ratioHist->SetMaximum(1.5);
  
  ratioHist->Draw();

  mcRatio->SetLineColor(kRed);
  //mcRatio->SetLineWidth(3);
  mcRatio->SetFillColorAlpha(kPink + 1, 0.4);
  mcRatio->Draw("E2 SAME");

  TH1D* straightLine = (TH1D*)mcRatio->Clone();
  straightLine->SetFillStyle(0);
  straightLine->Draw("HIST SAME");

  ratioHist->Draw("SAME");

  c1->Update();

  c1->Print(nameToSave+"/"+name+"_Efficiency_Comparison.pdf");
  c1->Print(nameToSave+"/"+name+"_Efficiency_Comparison.png");
  c1->Print(nameToSave+"/"+name+"_Efficiency_Comparison.C");

  delete MC_N;
  delete MC_NM1;
  delete MC_eff;
  delete data_N;
  delete data_NM1;
  delete data_eff;
  delete mcHist;
  delete errHist;
  delete dataHist;
  delete latex;
  delete leg;
  delete ratio;
  delete ratioHist;
  delete mcRatio;
  delete straightLine;
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
  if (argc < 7 || argc > 8) {
    cout << "Check usage..." << endl;
    return 2;
  }

  TString MC_N_fileName= (TString)(argv[1]);
  TString MC_NM1_fileName= (TString)(argv[2]);
  TString DATA_N_fileName= (TString)(argv[3]);
  TString DATA_NM1_fileName= (TString)(argv[4]);
  string outDir=string(argv[5]);
  TString title=(TString)(argv[6]);
  bool useSig = (argc==8) ? (atoi(argv[7])!=0) : false;
  
  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory doesn't exist. Exiting" << endl;
    return 3;
  }

  if (!MC_N_fileName.Contains(".root") || !MC_NM1_fileName.Contains(".root") || !DATA_N_fileName.Contains(".root") || !DATA_NM1_fileName.Contains(".root")){
    cout << "Each file must be a ROOT file" << endl;
    return 4;
  }

  cout << "MC N File Name: " << MC_N_fileName << endl;
  cout << "MC NM1 File Name: " << MC_NM1_fileName << endl;
  cout << "DATA N File Name: " << DATA_N_fileName << endl;
  cout << "DATA NM1 File Name: " << DATA_NM1_fileName << endl;
  
  TFile* mcFile_N = new TFile(MC_N_fileName,"READ");
  TFile* mcFile_NM1 = new TFile(MC_NM1_fileName,"READ");
  TFile* dataFile_N = new TFile(DATA_N_fileName,"READ");
  TFile* dataFile_NM1 = new TFile(DATA_NM1_fileName,"READ");

  TList* keyList = mcFile_N->GetListOfKeys();
  if (!keyList){
    cout << "List of keys failed to get." << endl;
    return 5;
  }

  TIter next(keyList);
  TKey* key;
  while ( key = (TKey*)next() ){
    TString name = (TString)(key->GetName());
    TString className = (TString)(key->GetClassName());
    if (!className.Contains("MnvH1D") || !name.Contains("_data") || name.Contains("vtx")) continue;
    DrawEffComp(name, mcFile_N, mcFile_NM1, dataFile_N, dataFile_NM1, title, (TString)(outDir.c_str()), useSig);
  }

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  mcFile_N->Close();
  mcFile_NM1->Close();
  dataFile_N->Close();
  dataFile_NM1->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
