//File: BKGFitting.cxx
//Info: This script is intended to fit recoil/pTmu plots using TMinuit primarily for the neutron selected sample.
//
//Usage: BKGFitting <mc_file> <data_file> <outdir> <recoilE/pTmu/vtxZ> <doSyst (only 0 means no)> optional: <mainTagName> <lowFitBinNum> <hiFitBinNum> <do fits in bins of muon momentum (only 0 means no)> TODO: Save the information beyond just printing it out
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

//TODO: Same SegFault Business From My Plotting Code... I'm assuming I just need to delete things carefully that I'm not yet.

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
#include "TFractionFitter.h"
#include "Math/IFunction.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TMinuitMinimizer.h"

//Analysis includes
#include "fits/ScaleFactors.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

using namespace std;
using namespace PlotUtils;

map<int, TString> colors = {{0,"#117733"},{1,"#CC6677"},{2,"#88CCEE"}};

bool PathExists(string path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

//Borrowed Directly from Andrew.                                                
void printCorrMatrix(const ROOT::Math::Minimizer& minim, const int nPars)
{
  vector<double> covMatrix(nPars * nPars, 0);
  minim.GetCovMatrix(covMatrix.data());
  const double* errors = minim.Errors();

  cout << "Printing Covariance Matrix" << endl;
  for(int xPar = 0; xPar < nPars; ++xPar){
    cout << "[";
    for(int yPar = 0; yPar < nPars-1; ++yPar){
      cout << fixed << setprecision(7) << setw(15) << covMatrix[xPar * nPars + yPar] << ", ";
    }
    cout << fixed << setprecision(7) << setw(15) << covMatrix[xPar * nPars + nPars-1] << "]\n";
  }

  cout << "Printing Correlation Matrix" << endl;
  for(int xPar = 0; xPar < nPars; ++xPar){
    cout << "[";
    for(int yPar = 0; yPar < nPars-1; ++yPar){
      cout << fixed << setprecision(2) << setw(5) << covMatrix[xPar * nPars + yPar]/errors[xPar]/errors[yPar] << ", ";
    }
    cout << fixed << setprecision(2) << setw(5) << covMatrix[xPar * nPars + nPars-1]/errors[xPar]/errors[nPars-1] << "]\n";
  }
}

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

//TODO: Instead of Plotting Save Scale Factor Histograms for analysis variables.
map<TString,MnvH1D*> FitScaleFactorsAndDraw(MnvH1D* dataHist, map<TString, MnvH1D*> fitHistsAndNames, map<TString, MnvH1D*> unfitHistsAndNames, TString varName, TString outDir, int lowBin, int hiBin, bool doSyst, bool sigFit){
  map<TString,MnvH1D*> scaleHists = {};

  TString name = varName+"_low_"+(TString)(to_string(lowBin))+"_hi_"+(TString)(to_string(hiBin));

  if (PathExists((string)(outDir+name+"_postFit.pdf"))){
    cout << "Already performed fits over this range for this histo." << endl;
    cout << "If you are doing this because of updated histos, it is in your best interest to save this elsewhere or remove the old plots." << endl;
    return scaleHists;
  }

  if (!PathExists((string)(outDir+varName+"_preFit.pdf"))){
    DrawFromMnvH1Ds(dataHist,fitHistsAndNames,unfitHistsAndNames,sigFit,outDir+varName+"_preFit");
  }

  TH1D* hData = (TH1D*)dataHist->GetCVHistoWithStatError().Clone();
  vector<TH1D*> fitHists = {};
  vector<TH1D*> unfitHists = {};

  for (auto hists:fitHistsAndNames){
    fitHists.push_back((TH1D*)hists.second->GetCVHistoWithStatError().Clone());
    scaleHists[hists.first] = new MnvH1D(name+hists.first,"",hists.second->GetNbinsX(),hists.second->GetXaxis()->GetXbins()->GetArray());
  }

  for (auto hists:unfitHistsAndNames){
    unfitHists.push_back((TH1D*)hists.second->GetCVHistoWithStatError().Clone());
  }

  fit::ScaleFactors func(fitHists,unfitHists,hData,lowBin,hiBin);

  auto* mini = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);

  int nextPar = 0;
  for (auto hist:fitHistsAndNames){
    string var = hist.first.Data();
    mini->SetVariable(nextPar,var,1.0,1.0);
    nextPar++;
  }

  if (nextPar != func.NDim()){
    cout << "The number of parameters was unexpected for some reason..." << endl;
    return scaleHists;
  }

  mini->SetFunction(func);

  if (!mini->Minimize()){
    cout << "FIT FAILED" << endl;
    cout << "Printing Results." << endl;
    mini->PrintResults();
    printCorrMatrix(*mini, func.NDim());
  }
  else{
    cout << "FIT SUCCEEDED" << endl;
    cout << "Printing Results." << endl;
    mini->PrintResults();
    printCorrMatrix(*mini, func.NDim());
  }
  
  const double* scaleResults = mini->X();
  const double* scaleErrors = mini->Errors();
  map<TString, double> scaleByName;
  nextPar=0;
  for (auto hist:fitHistsAndNames){
    scaleByName[hist.first]=scaleResults[nextPar];
    for (int iBin=0; iBin <= scaleHists[hist.first]->GetNbinsX()+1; ++iBin){
      scaleHists[hist.first]->SetBinContent(iBin,scaleResults[nextPar]);
      scaleHists[hist.first]->SetBinError(iBin,scaleErrors[nextPar]);
    }
    scaleHists[hist.first]->AddMissingErrorBandsAndFillWithCV(*hist.second);
    ++nextPar;
  }

  /*
  cout << "Checking the grabbing of scale factors." << endl;
  for (auto scale: scaleByName){
    cout << scale.first << " scale: " << scale.second << endl;
  }
  */

  //cout << "Trying to draw pre-scaling." << endl;

  cout << "Scaling the entire MnvH1D first." << endl;
  for (auto hist:fitHistsAndNames){
    hist.second->Scale(scaleByName[hist.first]);
  }

  if (doSyst){
    cout << "Looping over systematic universes." << endl;
    
    //Info here largely lifted from Andrew's fitting script.
    map<TString, TH1D*> scaleFactorHists;
    int nUnivBins = 180;
    
    nextPar=0;
    for (auto hist:fitHistsAndNames){
      scaleFactorHists[hist.first] = new TH1D(name+"_ScaleFactorHists_par"+(TString)(to_string(nextPar)),"Scale Factor away from CV for "+hist.first+" vs. Universe;Universe No.;Scale Factor",nUnivBins,0,nUnivBins);
      ++nextPar;
    }
    
    double quelUniv=0.5;
    
    const auto errorBandNames = fitHistsAndNames.begin()->second->GetErrorBandNames();
    for (const auto& bandName: errorBandNames){
      cout << bandName << endl;
      const auto univs = fitHistsAndNames.begin()->second->GetVertErrorBand(bandName)->GetHists();
      for (size_t whichUniv=0; whichUniv < univs.size(); ++ whichUniv){
	cout << "" << endl;
	cout << "Fitting Universe " << whichUniv << " in " << bandName << " error band." << endl;
	
	vector<TH1D*> fitHistsUniv = {};
	vector<TH1D*> unfitHistsUniv = {};
	
	for (auto hists:fitHistsAndNames){
	  fitHistsUniv.push_back((TH1D*)hists.second->GetVertErrorBand(bandName)->GetHist(whichUniv)->Clone());
	}
	
	for (auto hists:unfitHistsAndNames){
	  unfitHistsUniv.push_back((TH1D*)hists.second->GetVertErrorBand(bandName)->GetHist(whichUniv)->Clone());
	}
	
	fit::ScaleFactors funcUniv(fitHistsUniv,unfitHistsUniv,hData,lowBin,hiBin);
	
	auto* miniUniv = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);
	
	int parNext = 0;
	for (auto hist:fitHistsAndNames){
	  string var = hist.first.Data();
	  miniUniv->SetVariable(parNext,var,1.0,1.0);
	  parNext++;
	}
	
	if (parNext != funcUniv.NDim()){
	  cout << "The number of parameters was unexpected for some reason..." << endl;
	  return scaleHists;
	}
	
	miniUniv->SetFunction(funcUniv);
	
	if (!miniUniv->Minimize()){
	  cout << "FIT FAILED" << endl;
	  cout << "Printing Results." << endl;
	  miniUniv->PrintResults();
	  printCorrMatrix(*miniUniv, funcUniv.NDim());
	}
	else{
	  cout << "FIT SUCCEEDED" << endl;
	  cout << "Printing Results." << endl;
	  miniUniv->PrintResults();
	  printCorrMatrix(*miniUniv, funcUniv.NDim());
	}
	
	const double* scaleResultsUniv = miniUniv->X();
	parNext=0;
	for (auto hist:fitHistsAndNames){
	  hist.second->GetVertErrorBand(bandName)->GetHist(whichUniv)->Scale(scaleResultsUniv[parNext]);	
	  scaleHists[hist.first]->GetVertErrorBand(bandName)->GetHist(whichUniv)->Scale(scaleResultsUniv[parNext]);
	  scaleFactorHists[hist.first]->Fill(quelUniv,scaleResultsUniv[parNext]);
	  ++parNext;
	}
	
	quelUniv += 1.0;
	for (auto hist:fitHistsUniv) delete hist;
	fitHistsUniv.clear();
	for (auto hist:unfitHistsUniv) delete hist;
	unfitHistsUniv.clear();
      }
    }

    TCanvas* c1 = new TCanvas("c1","c1",1200,800);
    c1->cd();
    for (auto hist:scaleFactorHists){
      hist.second->GetYaxis()->SetRangeUser(-3.0,3.0);
      hist.second->Draw("hist");
      c1->Print(outDir+hist.second->GetName()+".pdf");
      c1->Print(outDir+hist.second->GetName()+".png");
    }
    delete c1;
  }

  DrawFromMnvH1Ds(dataHist,fitHistsAndNames,unfitHistsAndNames,sigFit,outDir+name+"_postFit");
  
  cout << "Drawing of CV scale completed." << endl;

  delete hData;
  for (auto hist:fitHists) delete hist;
  fitHists.clear();
  for (auto hist:unfitHists) delete hist;
  unfitHists.clear();

  return scaleHists;
}

int main(int argc, char* argv[]) {

  gStyle->SetOptStat(0);

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc < 6 || argc > 10) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string MCfileName = string(argv[1]);
  string DATAfileName = string(argv[2]);
  string outDir = string(argv[3]);
  TString varName= argv[4];
  bool doSyst = (bool)(atoi(argv[5]));
  int fitMuonBins = 0;

  int lowBin = 1;//Will be truncated later. Lowest allowed value for pT or recoil fits.
  int hiBin = 50;//Will be truncated later. Highest allowed value for pT or recoil fits.
  TString mainTag = "";
  if (argc > 7) lowBin = max(atoi(argv[7]),1);//Not allowed lower than 1
  if (argc > 8) hiBin = min(50, atoi(argv[8]));//Not allowed higher than 50
  if (argc > 9) fitMuonBins = atoi(argv[9]);

  string rootExt = ".root";
  string slash = "/";
  string token;
  string fileNameStub = MCfileName;
  size_t pos=0;

  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory: " << outDir << " doesn't exist. Exiting" << endl;
    return 3;
  }

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

  if (varName == "pTmu"){
    if (lowBin >= 14){
      cout << "Not a valid fitting range for pTmu. Maximum bin is 14." << endl;
      return 100;
    }
    hiBin = min(14,hiBin);
    mainTag = "_RecoilSB";
  }
  else if (varName == "recoilE"){
    if (argc < 8) lowBin = 26;
    mainTag = "_PreRecoilCut";
  }
  else {
    cout << "Not a valid variable name to fit." << endl;
    return 111;
  }
  if (argc > 6) mainTag = argv[6];

  //NEED TO CODE IN SOMETHING THAT HANDLES THE VERTEX PLOTS IN THE TARGETS.

  //cout << "Setting up MnvPlotter" << endl;
  //MnvPlotter* plotter = new MnvPlotter(kCCQEAntiNuStyle);

  TFile* outFile = new TFile("fitResults.root","RECREATE");
  TFile* mcFile = new TFile(MCfileName.c_str(),"READ");
  TFile* dataFile = new TFile(DATAfileName.c_str(),"READ");

  TParameter<double>* mcPOT = (TParameter<double>*)mcFile->Get("POTUsed");
  TParameter<double>* dataPOT = (TParameter<double>*)dataFile->Get("POTUsed");

  vector<TString> tags = {mainTag};
  if (fitMuonBins){
    tags.push_back((TString)("_bin_lost"));
    for (int iBin=0; iBin < 14; ++iBin){
      TString tag = "_bin_"+to_string(iBin);
      //cout << tag << endl;
      tags.push_back(tag);
    }
  }

  double POTscale = dataPOT->GetVal()/mcPOT->GetVal();
  //cout << "POT scale factor: " << scale << endl;
  
  //TString varName = "recoilE";

  for (int iTag=0; iTag < tags.size(); ++iTag){

    TString tag = tags.at(iTag);
    TString name = varName+tag;

    cout << "Performing Fitting and Scaling for: " << name << endl;

    MnvH1D* dataHist = (MnvH1D*)(dataFile->Get(name+"_data"))->Clone();
    MnvH1D* sigHist = (MnvH1D*)(mcFile->Get(name+"_selected_signal_reco"))->Clone();
    sigHist->Scale(POTscale);

    //
    MnvH1D* chargePiHist = (MnvH1D*)(mcFile->Get(name+"_background_1chargePi"))->Clone();
    chargePiHist->Scale(POTscale);
    MnvH1D* neutPiHist = (MnvH1D*)(mcFile->Get(name+"_background_1neutPi"))->Clone();
    neutPiHist->Scale(POTscale);
    MnvH1D* NPiHist = (MnvH1D*)(mcFile->Get(name+"_background_NPi"))->Clone();
    NPiHist->Scale(POTscale);
    MnvH1D* otherHist = (MnvH1D*)(mcFile->Get(name+"_background_Other"))->Clone();
    otherHist->Scale(POTscale);

    MnvH1D* bkgNNeutPiHist = neutPiHist->Clone();
    bkgNNeutPiHist->Add(NPiHist);

    MnvH1D* bkg1PiHist = neutPiHist->Clone();
    bkg1PiHist->Add(chargePiHist);

    //
    MnvH1D* QEHist = (MnvH1D*)(mcFile->Get(name+"_bkg_IntType_QE"))->Clone();
    QEHist->Scale(POTscale);
    MnvH1D* RESHist = (MnvH1D*)(mcFile->Get(name+"_bkg_IntType_RES"))->Clone();
    RESHist->Scale(POTscale);
    MnvH1D* DISHist = (MnvH1D*)(mcFile->Get(name+"_bkg_IntType_DIS"))->Clone();
    DISHist->Scale(POTscale);
    MnvH1D* MECHist = (MnvH1D*)(mcFile->Get(name+"_bkg_IntType_2p2h"))->Clone();
    MECHist->Scale(POTscale);
    MnvH1D* OtherIntTypeHist = (MnvH1D*)(mcFile->Get(name+"_bkg_IntType_Other"))->Clone();
    OtherIntTypeHist->Scale(POTscale);

    MnvH1D* bkgNonRESHist = DISHist->Clone();
    bkgNonRESHist->Add(QEHist);
    bkgNonRESHist->Add(MECHist);
    bkgNonRESHist->Add(OtherIntTypeHist);

    MnvH1D* bkgTotHist = bkgNonRESHist->Clone();
    bkgTotHist->Add(RESHist);

    map<TString, MnvH1D*> fitHists1A, unfitHists1A;
    /*
    map<TString, MnvH1D*> fitHists2A, unfitHists2A;
    map<TString, MnvH1D*> fitHists3A, unfitHists3A;
    map<TString, MnvH1D*> fitHists4A, unfitHists4A;
    map<TString, MnvH1D*> fitHists5A, unfitHists5A;
    map<TString, MnvH1D*> fitHists6A, unfitHists6A;
    */

    map<TString, MnvH1D*> fitHists1B, unfitHists1B;
    /*
    map<TString, MnvH1D*> fitHists2B, unfitHists2B;
    map<TString, MnvH1D*> fitHists3B, unfitHists3B;
    map<TString, MnvH1D*> fitHists4B, unfitHists4B;
    map<TString, MnvH1D*> fitHists5B, unfitHists5B;
    map<TString, MnvH1D*> fitHists6B, unfitHists6B;
    */

    fitHists1A["BKG"]=(MnvH1D*)bkgTotHist->Clone();
    fitHists1A["Signal"]=(MnvH1D*)sigHist->Clone();

    fitHists1B["BKG"]=(MnvH1D*)bkgTotHist->Clone();
    unfitHists1B["Signal"]=(MnvH1D*)sigHist->Clone();

    /*
    fitHists2A["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    fitHists2A["single #pi^{0}"]=(MnvH1D*)neutPiHist->Clone();
    fitHists2A["N#pi"]=(MnvH1D*)NPiHist->Clone();
    fitHists2A["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists2A["Other"]=(MnvH1D*)otherHist->Clone();

    fitHists2B["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    fitHists2B["single #pi^{0}"]=(MnvH1D*)neutPiHist->Clone();
    fitHists2B["N#pi"]=(MnvH1D*)NPiHist->Clone();
    unfitHists2B["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists2B["Other"]=(MnvH1D*)otherHist->Clone();

    fitHists3A["single #pi"]=(MnvH1D*)bkg1PiHist->Clone();
    fitHists3A["N#pi"]=(MnvH1D*)NPiHist->Clone();
    fitHists3A["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists3A["Other"]=(MnvH1D*)otherHist->Clone();

    fitHists3B["single #pi"]=(MnvH1D*)bkg1PiHist->Clone();
    fitHists3B["N#pi"]=(MnvH1D*)NPiHist->Clone();
    unfitHists3B["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists3B["Other"]=(MnvH1D*)otherHist->Clone();

    fitHists4A["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    fitHists4A["N#pi & single #pi^{0}"]=(MnvH1D*)bkgNNeutPiHist->Clone();
    fitHists4A["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists4A["Other"]=(MnvH1D*)otherHist->Clone();

    fitHists4B["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    fitHists4B["N#pi & single #pi^{0}"]=(MnvH1D*)bkgNNeutPiHist->Clone();
    unfitHists4B["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists4B["Other"]=(MnvH1D*)otherHist->Clone();

    fitHists5A["RES"]=(MnvH1D*)RESHist->Clone();
    fitHists5A["non-RES"]=(MnvH1D*)bkgNonRESHist->Clone();
    fitHists5A["Signal"]=(MnvH1D*)sigHist->Clone();

    fitHists5B["RES"]=(MnvH1D*)RESHist->Clone();
    fitHists5B["non-RES"]=(MnvH1D*)bkgNonRESHist->Clone();
    unfitHists5B["Signal"]=(MnvH1D*)sigHist->Clone();

    fitHists6A["RES"]=(MnvH1D*)RESHist->Clone();
    fitHists6A["DIS"]=(MnvH1D*)DISHist->Clone();
    fitHists6A["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists6A["QE"]=(MnvH1D*)QEHist->Clone();
    unfitHists6A["2p2h"]=(MnvH1D*)MECHist->Clone();
    unfitHists6A["Other"]=(MnvH1D*)OtherIntTypeHist->Clone();

    fitHists6B["RES"]=(MnvH1D*)RESHist->Clone();
    fitHists6B["DIS"]=(MnvH1D*)DISHist->Clone();
    unfitHists6B["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists6B["QE"]=(MnvH1D*)QEHist->Clone();
    unfitHists6B["2p2h"]=(MnvH1D*)MECHist->Clone();
    unfitHists6B["Other"]=(MnvH1D*)OtherIntTypeHist->Clone();
    */

    cout << "Fitting 1A" << endl;
    map<TString,MnvH1D*> result = FitScaleFactorsAndDraw(dataHist, fitHists1A, unfitHists1A, name+"_fit1A", outDir, lowBin, hiBin, doSyst, true);
    map<TString,MnvH1D*> scaledHists1A = {};
    scaledHists1A["BKG"]=(MnvH1D*)bkgTotHist->Clone();
    scaledHists1A["Signal"]=(MnvH1D*)sigHist->Clone();
    for(auto hists:scaledHists1A){
      hists.second->Multiply(hists.second,result[hists.first]);
      result[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists1A,unfitHists1A,true,outDir+"TEST_NEW_"+name+"_fit1A_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result) delete hist.second;
    result.clear();

    cout << "Fitting 1B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists1B, unfitHists1B, name+"_fit1B", outDir, lowBin, hiBin, doSyst, false);
    map<TString,MnvH1D*> scaledHists1B = {};
    scaledHists1B["BKG"]=(MnvH1D*)bkgTotHist->Clone();
    for(auto hists:scaledHists1B){
      hists.second->Multiply(hists.second,result[hists.first]);
      result[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists1B,unfitHists1B,false,outDir+"TEST_NEW_"+name+"_fit1B_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result) delete hist.second;
    result.clear();

    /*
    cout << "Fitting 2A" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists2A, unfitHists2A, name+"_fit2A", outDir, lowBin, hiBin, doSyst, true);
    cout << "Result: " << result << endl;
    cout << "" << endl;

    cout << "Fitting 2B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists2B, unfitHists2B, name+"_fit2B", outDir, lowBin, hiBin, doSyst, false);
    cout << "Result: " << result << endl;
    cout << "" << endl;

    cout << "Fitting 3A" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists3A, unfitHists3A, name+"_fit3A", outDir, lowBin, hiBin, doSyst, true);
    cout << "Result: " << result << endl;
    cout << "" << endl;

    cout << "Fitting 3B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists3B, unfitHists3B, name+"_fit3B", outDir, lowBin, hiBin, doSyst, false);
    cout << "Result: " << result << endl;
    cout << "" << endl;

    cout << "Fitting 4A" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists4A, unfitHists4A, name+"_fit4A", outDir, lowBin, hiBin, doSyst, true);
    cout << "Result: " << result << endl;
    cout << "" << endl;

    cout << "Fitting 4B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists4B, unfitHists4B, name+"_fit4B", outDir, lowBin, hiBin, doSyst, false);
    cout << "Result: " << result << endl;
    cout << "" << endl;

    cout << "Fitting 5A" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists5A, unfitHists5A, name+"_fit5A", outDir, lowBin, hiBin, doSyst, true);
    cout << "Result: " << result << endl;
    cout << "" << endl;

    cout << "Fitting 5B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists5B, unfitHists5B, name+"_fit5B", outDir, lowBin, hiBin, doSyst, false);
    cout << "Result: " << result << endl;
    cout << "" << endl;

    cout << "Fitting 6A" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists6A, unfitHists6A, name+"_fit6A", outDir, lowBin, hiBin, doSyst, true);
    cout << "Result: " << result << endl;
    cout << "" << endl;

    cout << "Fitting 6B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists6B, unfitHists6B, name+"_fit6B", outDir, lowBin, hiBin, doSyst, false);
    cout << "Result: " << result << endl;
    cout << "" << endl;
    //if (result != 0) return result;
    */

    delete dataHist;
    delete sigHist;
    delete chargePiHist;
    delete neutPiHist;
    delete NPiHist;
    delete otherHist;
    delete bkgNNeutPiHist;
    delete bkg1PiHist;
    delete QEHist;
    delete RESHist;
    delete DISHist;
    delete MECHist;
    delete OtherIntTypeHist;
    delete bkgNonRESHist;
    delete bkgTotHist;
    for (auto hist:fitHists1A) delete hist.second;
    for (auto hist:unfitHists1A) delete hist.second;
    for (auto hist:fitHists1B) delete hist.second;
    for (auto hist:unfitHists1B) delete hist.second;
    for (auto hist:scaledHists1A) delete hist.second;
    for (auto hist:scaledHists1B) delete hist.second;
    /*
    for (auto hist:fitHists2A) delete hist.second;
    for (auto hist:unfitHists2A) delete hist.second;
    for (auto hist:fitHists2B) delete hist.second;
    for (auto hist:unfitHists2B) delete hist.second;
    for (auto hist:fitHists3A) delete hist.second;
    for (auto hist:unfitHists3A) delete hist.second;
    for (auto hist:fitHists3B) delete hist.second;
    for (auto hist:unfitHists3B) delete hist.second;
    for (auto hist:fitHists4A) delete hist.second;
    for (auto hist:unfitHists4A) delete hist.second;
    for (auto hist:fitHists4B) delete hist.second;
    for (auto hist:unfitHists4B) delete hist.second;
    for (auto hist:fitHists5A) delete hist.second;
    for (auto hist:unfitHists5A) delete hist.second;
    for (auto hist:fitHists5B) delete hist.second;
    for (auto hist:unfitHists5B) delete hist.second;
    for (auto hist:fitHists6A) delete hist.second;
    for (auto hist:unfitHists6A) delete hist.second;
    for (auto hist:fitHists6B) delete hist.second;
    for (auto hist:unfitHists6B) delete hist.second;
    */
  }

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  mcFile->Close();
  dataFile->Close();
  outFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
