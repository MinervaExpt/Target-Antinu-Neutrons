//File: BKGFitting.cxx
//Info: This script is intended to fit recoil/pTmu plots using TMinuit primarily for the neutron selected sample.
//
//Usage: BKGFitting <mc_file> <data_file> <outdir> <outFileTag> <recoilE/pTmu/vtxZ> <doSyst (only 0 means no)> <Tgts> optional: <mainTagName> <lowFitBinNum> <hiFitBinNum> <do fits in bins of muon momentum (only 0 means no)> TODO: Save the information beyond just printing it out
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
#include "TObjString.h"
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

map<TString,map<TString,MnvH1D*>> FitScaleFactorsAndDraw(MnvH1D* dataHist, map<TString, MnvH1D*> fitHistsAndNames, map<TString, MnvH1D*> unfitHistsAndNames, TString varName, TString outDir, TString fitName, int lowBin, int hiBin, bool doSyst, bool sigFit, map<TString,MnvH1D*> varsToSave, map<TString,vector<TString>> tagsToSave){
  map<TString,map<TString,MnvH1D*>> scaleHists = {};

  TString nameTag = fitName+"_low_"+(TString)(to_string(lowBin))+"_hi_"+(TString)(to_string(hiBin));
  TString name = varName+nameTag;

  if (PathExists((string)(outDir+name+"_postFit.pdf"))){
    cout << "Already performed fits over this range for this histo." << endl;
    cout << "If you are doing this because of updated histos, it is in your best interest to save this elsewhere or remove the old plots." << endl;
    return scaleHists;
  }

  if (!PathExists((string)(outDir+varName+fitName+"_preFit.pdf"))){
    DrawFromMnvH1Ds(dataHist,fitHistsAndNames,unfitHistsAndNames,sigFit,outDir+varName+fitName+"_preFit");
  }

  cout << "Getting Data" << endl;
  TH1D* hData = (TH1D*)dataHist->GetCVHistoWithStatError().Clone();
  vector<TH1D*> fitHists = {};
  vector<TH1D*> unfitHists = {};

  for (auto hists:fitHistsAndNames){
    TString dumpTag = "";//tag which lets you know which histos to scale later.
    for (auto tag:tagsToSave[hists.first]){
      dumpTag = dumpTag+"_t_"+tag;
    }
    fitHists.push_back((TH1D*)hists.second->GetCVHistoWithStatError().Clone());
    for (auto var:varsToSave){
      if (var.first == varName) continue;
      scaleHists[var.first][hists.first] = new MnvH1D(var.first+"_fit_"+name+dumpTag,"",var.second->GetNbinsX(),var.second->GetXaxis()->GetXbins()->GetArray());
    }
    scaleHists[varName][hists.first]= new MnvH1D(varName+"_fit_"+name+dumpTag,"",hists.second->GetNbinsX(),hists.second->GetXaxis()->GetXbins()->GetArray());
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
    cout << "HI" << endl;
    scaleByName[hist.first]=scaleResults[nextPar];
    for (auto var:scaleHists){
      cout << "HI 2" << endl;
      for (int iBin=0; iBin <= var.second[hist.first]->GetNbinsX()+1; ++iBin){
	cout << "HI 3" << endl;
	var.second[hist.first]->SetBinContent(iBin,scaleResults[nextPar]);
	var.second[hist.first]->SetBinError(iBin,scaleErrors[nextPar]);
      }
      cout << "HI 4" << endl;
      cout << hist.first << endl;
      var.second[hist.first]->AddMissingErrorBandsAndFillWithCV(*hist.second);
    }
    cout << "HI 5" << endl;
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
	  for (auto var:scaleHists) var.second[hist.first]->GetVertErrorBand(bandName)->GetHist(whichUniv)->Scale(scaleResultsUniv[parNext]);
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
  if (argc < 8 || argc > 12) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string MCfileName = string(argv[1]);
  string DATAfileName = string(argv[2]);
  string outDir = string(argv[3]);
  string outFileTag = string(argv[4]);
  TString varName= argv[5];
  bool doSyst = (bool)(atoi(argv[6]));
  bool Tgts = (bool)(atoi(argv[7]));
  int fitMuonBins = 0;
  
  vector<TString> namesToSave = {"pTmu","recoilE"};
  //vector<TString> namesToSave = {"pTmu","recoilE","NPlanes"};
  //vector<TString> namesToSave = {};

  int lowBin = 1;//Will be truncated later. Lowest allowed value for pT or recoil fits.
  int hiBin = 50;//Will be truncated later. Highest allowed value for pT or recoil fits.
  TString mainTag = "";
  if (argc > 9) lowBin = max(atoi(argv[9]),1);//Not allowed lower than 1
  if (argc > 10) hiBin = min(50, atoi(argv[10]));//Not allowed higher than 50
  if (argc > 11) fitMuonBins = atoi(argv[11]);

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
    if (argc < 9) lowBin = 26;
    mainTag = "_PreRecoilCut";
  }
  else {
    cout << "Not a valid variable name to fit." << endl;
    return 111;
  }
  if (argc > 8) mainTag = argv[8];

  //NEED TO CODE IN SOMETHING THAT HANDLES THE VERTEX PLOTS IN THE TARGETS.

  //cout << "Setting up MnvPlotter" << endl;
  //MnvPlotter* plotter = new MnvPlotter(kCCQEAntiNuStyle);

  TFile* outFile = new TFile((outDir+"fitResults_"+outFileTag+".root").c_str(),"RECREATE");
  TFile* mcFile = new TFile(MCfileName.c_str(),"READ");
  TFile* dataFile = new TFile(DATAfileName.c_str(),"READ");

  TParameter<double>* mcPOT = (TParameter<double>*)mcFile->Get("POTUsed");
  TParameter<double>* dataPOT = (TParameter<double>*)dataFile->Get("POTUsed");

  vector<TString> tags = {mainTag};

  if (Tgts){
    tags.push_back(mainTag+"_Tgt1");
    tags.push_back(mainTag+"_Tgt2");
    tags.push_back(mainTag+"_Tgt3");
    tags.push_back(mainTag+"_Tgt4");
    tags.push_back(mainTag+"_Tgt5");
  }

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
    TString grabName;
    if (tag.Contains("_Tgt")){
      TObjArray* tagArr = tag.Tokenize("Tgt");
      grabName="ByTgt_Tgt"+((TObjString*)(tagArr->At(tagArr->GetEntries()-1)))->String()+"/"+name;
    }
    else grabName = name;

    cout << "Performing Fitting and Scaling for: " << name << endl;
    cout << "Grabbing: " << grabName << endl;

    double binsNPlanes[22];
    for (unsigned int i=0; i < 22;++i) binsNPlanes[i]=1.0*i-10.5;

    map<TString,MnvH1D*> varsToSave = {};
    for (auto nameSave:namesToSave){
      if (nameSave == "NPlanes"){
	MnvH1D* hNPlanes = new MnvH1D(nameSave+tag,"",21,binsNPlanes);
	varsToSave[nameSave+tag] = hNPlanes->Clone();
      }
      else if (tag.Contains("_Tgt")){
	TObjArray* tagArr = tag.Tokenize("Tgt");
	TString grabNameSave = "ByTgt_Tgt"+((TObjString*)(tagArr->At(tagArr->GetEntries()-1)))->String()+"/"+nameSave+tag;
	varsToSave[nameSave+tag] = (MnvH1D*)(mcFile->Get(grabNameSave+"_selected_signal_reco"))->Clone();
      }
      else varsToSave[nameSave+tag] = (MnvH1D*)(mcFile->Get(nameSave+tag+"_selected_signal_reco"))->Clone();
    }

    MnvH1D* dataHist = (MnvH1D*)(dataFile->Get(grabName+"_data"))->Clone();
    MnvH1D* sigHist = (MnvH1D*)(mcFile->Get(grabName+"_selected_signal_reco"))->Clone();
    sigHist->Scale(POTscale);

    //
    MnvH1D* chargePiHist = (MnvH1D*)(mcFile->Get(grabName+"_background_1chargePi"))->Clone();
    chargePiHist->Scale(POTscale);
    MnvH1D* neutPiHist = (MnvH1D*)(mcFile->Get(grabName+"_background_1neutPi"))->Clone();
    neutPiHist->Scale(POTscale);
    MnvH1D* NPiHist = (MnvH1D*)(mcFile->Get(grabName+"_background_NPi"))->Clone();
    NPiHist->Scale(POTscale);
    MnvH1D* otherHist = (MnvH1D*)(mcFile->Get(grabName+"_background_Other"))->Clone();
    otherHist->Scale(POTscale);

    MnvH1D* bkgNNeutPiHist = neutPiHist->Clone();
    bkgNNeutPiHist->Add(NPiHist);

    MnvH1D* bkg1PiHist = neutPiHist->Clone();
    bkg1PiHist->Add(chargePiHist);

    //
    MnvH1D* QEHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_QE"))->Clone();
    QEHist->Scale(POTscale);
    MnvH1D* RESHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_RES"))->Clone();
    RESHist->Scale(POTscale);
    MnvH1D* DISHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_DIS"))->Clone();
    DISHist->Scale(POTscale);
    MnvH1D* MECHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_2p2h"))->Clone();
    MECHist->Scale(POTscale);
    MnvH1D* OtherIntTypeHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_Other"))->Clone();
    OtherIntTypeHist->Scale(POTscale);

    MnvH1D* bkgNonRESHist = DISHist->Clone();
    bkgNonRESHist->Add(QEHist);
    bkgNonRESHist->Add(MECHist);
    bkgNonRESHist->Add(OtherIntTypeHist);

    MnvH1D* bkgTotHist = bkgNonRESHist->Clone();
    bkgTotHist->Add(RESHist);

    map<TString, MnvH1D*> fitHists1A, unfitHists1A;
    map<TString, MnvH1D*> fitHists2A, unfitHists2A;
    map<TString, MnvH1D*> fitHists3A, unfitHists3A;
    map<TString, MnvH1D*> fitHists4A, unfitHists4A;
    map<TString, MnvH1D*> fitHists5A, unfitHists5A;
    map<TString, MnvH1D*> fitHists6A, unfitHists6A;

    map<TString, MnvH1D*> fitHists1B, unfitHists1B;
    map<TString, MnvH1D*> fitHists2B, unfitHists2B;
    map<TString, MnvH1D*> fitHists3B, unfitHists3B;
    map<TString, MnvH1D*> fitHists4B, unfitHists4B;
    map<TString, MnvH1D*> fitHists5B, unfitHists5B;
    map<TString, MnvH1D*> fitHists6B, unfitHists6B;

    map<TString, vector<TString>> nameKeys1A, nameKeys1B;
    map<TString, vector<TString>> nameKeys2A, nameKeys2B;
    map<TString, vector<TString>> nameKeys3A, nameKeys3B;
    map<TString, vector<TString>> nameKeys4A, nameKeys4B;
    map<TString, vector<TString>> nameKeys5A, nameKeys5B;
    map<TString, vector<TString>> nameKeys6A, nameKeys6B;

    fitHists1A["BKG"]=(MnvH1D*)bkgTotHist->Clone();
    fitHists1A["Signal"]=(MnvH1D*)sigHist->Clone();
    nameKeys1A["BKG"]={"bkg","background"};
    nameKeys1A["Signal"]={"sig","signal"};

    fitHists1B["BKG"]=(MnvH1D*)bkgTotHist->Clone();
    unfitHists1B["Signal"]=(MnvH1D*)sigHist->Clone();
    nameKeys1B["BKG"]=nameKeys1A["BKG"];

    fitHists2A["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    fitHists2A["single #pi^{0}"]=(MnvH1D*)neutPiHist->Clone();
    fitHists2A["N#pi"]=(MnvH1D*)NPiHist->Clone();
    fitHists2A["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists2A["Other"]=(MnvH1D*)otherHist->Clone();
    nameKeys2A["single #pi^{#pm}"]={"background_1chargePi"};
    nameKeys2A["single #pi^{0}"]={"background_1neutPi"};
    nameKeys2A["N#pi"]={"background_NPi"};
    nameKeys2A["Signal"]=nameKeys1A["Signal"];

    fitHists2B["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    fitHists2B["single #pi^{0}"]=(MnvH1D*)neutPiHist->Clone();
    fitHists2B["N#pi"]=(MnvH1D*)NPiHist->Clone();
    unfitHists2B["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists2B["Other"]=(MnvH1D*)otherHist->Clone();
    nameKeys2B["single #pi^{#pm}"]=nameKeys2A["single #pi^{#pm}"];
    nameKeys2B["single #pi^{0}"]=nameKeys2A["single #pi^{0}"];
    nameKeys2B["N#pi"]=nameKeys2A["N#pi"];

    fitHists3A["single #pi"]=(MnvH1D*)bkg1PiHist->Clone();
    fitHists3A["N#pi"]=(MnvH1D*)NPiHist->Clone();
    fitHists3A["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists3A["Other"]=(MnvH1D*)otherHist->Clone();
    nameKeys3A["Signal"]=nameKeys1A["Signal"];
    nameKeys3A["single #pi"]={"background_1chargePi","1neutPi"};
    nameKeys3A["N#pi"]={"background_NPi"};

    fitHists3B["single #pi"]=(MnvH1D*)bkg1PiHist->Clone();
    fitHists3B["N#pi"]=(MnvH1D*)NPiHist->Clone();
    unfitHists3B["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists3B["Other"]=(MnvH1D*)otherHist->Clone();
    nameKeys3B["single #pi"]=nameKeys3A["single #pi"];
    nameKeys3B["N#pi"]=nameKeys3A["N#pi"];

    fitHists4A["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    fitHists4A["N#pi & single #pi^{0}"]=(MnvH1D*)bkgNNeutPiHist->Clone();
    fitHists4A["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists4A["Other"]=(MnvH1D*)otherHist->Clone();
    nameKeys4A["Signal"]=nameKeys1A["Signal"];
    nameKeys4A["single #pi^{#pm}"]={"background_1chargePi"};
    nameKeys4A["N#pi & single #pi^{0}"]={"background_NPi","1neutPi"};

    fitHists4B["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    fitHists4B["N#pi & single #pi^{0}"]=(MnvH1D*)bkgNNeutPiHist->Clone();
    unfitHists4B["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists4B["Other"]=(MnvH1D*)otherHist->Clone();
    nameKeys4B["single #pi^{#pm}"]=nameKeys4A["single #pi^{#pm}"];
    nameKeys4B["N#pi & single #pi^{0}"]=nameKeys4A["N#pi & single #pi^{0}"];

    fitHists5A["RES"]=(MnvH1D*)RESHist->Clone();
    fitHists5A["nonRES"]=(MnvH1D*)bkgNonRESHist->Clone();
    fitHists5A["Signal"]=(MnvH1D*)sigHist->Clone();
    nameKeys5A["Signal"]=nameKeys1A["Signal"];
    nameKeys5A["RES"]={"bkg_IntType_RES"};
    nameKeys5A["nonRES"]={"IntType_DIS","2p2h","Other","Wrong_Nucleus","USPlastic","DSPlastic"};

    fitHists5B["RES"]=(MnvH1D*)RESHist->Clone();
    fitHists5B["nonRES"]=(MnvH1D*)bkgNonRESHist->Clone();
    unfitHists5B["Signal"]=(MnvH1D*)sigHist->Clone();
    nameKeys5B["RES"]=nameKeys5A["RES"];
    nameKeys5B["nonRES"]=nameKeys5A["nonRES"];

    fitHists6A["RES"]=(MnvH1D*)RESHist->Clone();
    fitHists6A["DIS"]=(MnvH1D*)DISHist->Clone();
    fitHists6A["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists6A["QE"]=(MnvH1D*)QEHist->Clone();
    unfitHists6A["2p2h"]=(MnvH1D*)MECHist->Clone();
    unfitHists6A["Other"]=(MnvH1D*)OtherIntTypeHist->Clone();
    nameKeys6A["Signal"]=nameKeys1A["Signal"];
    nameKeys6A["RES"]={"bkg_IntType_RES"};
    nameKeys6A["DIS"]={"bkg_IntType_DIS"};

    fitHists6B["RES"]=(MnvH1D*)RESHist->Clone();
    fitHists6B["DIS"]=(MnvH1D*)DISHist->Clone();
    unfitHists6B["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists6B["QE"]=(MnvH1D*)QEHist->Clone();
    unfitHists6B["2p2h"]=(MnvH1D*)MECHist->Clone();
    unfitHists6B["Other"]=(MnvH1D*)OtherIntTypeHist->Clone();
    nameKeys6B["RES"]=nameKeys6A["RES"];
    nameKeys6B["DIS"]=nameKeys6A["DIS"];

    cout << "Fitting 1A" << endl;
    map<TString,map<TString,MnvH1D*>> result = FitScaleFactorsAndDraw(dataHist, fitHists1A, unfitHists1A, name, outDir, "_fit1A", lowBin, hiBin, doSyst, true, varsToSave, nameKeys1A);
    map<TString,MnvH1D*> scaledHists1A = {};
    scaledHists1A["BKG"]=(MnvH1D*)bkgTotHist->Clone();
    scaledHists1A["Signal"]=(MnvH1D*)sigHist->Clone();
    for(auto hists:scaledHists1A){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists1A,unfitHists1A,true,outDir+"TEST_NEW_"+name+"_fit1A_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 1B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists1B, unfitHists1B, name, outDir, "_fit1B", lowBin, hiBin, doSyst, false, varsToSave, nameKeys1B);
    map<TString,MnvH1D*> scaledHists1B = {};
    scaledHists1B["BKG"]=(MnvH1D*)bkgTotHist->Clone();
    for(auto hists:scaledHists1B){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists1B,unfitHists1B,false,outDir+"TEST_NEW_"+name+"_fit1B_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 2A" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists2A, unfitHists2A, name, outDir, "_fit2A", lowBin, hiBin, doSyst, true, varsToSave, nameKeys2A);
    map<TString,MnvH1D*> scaledHists2A = {};
    scaledHists2A["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    scaledHists2A["single #pi^{0}"]=(MnvH1D*)neutPiHist->Clone();
    scaledHists2A["N#pi"]=(MnvH1D*)NPiHist->Clone();
    scaledHists2A["Signal"]=(MnvH1D*)sigHist->Clone();
    for(auto hists:scaledHists2A){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists2A,unfitHists2A,true,outDir+"TEST_NEW_"+name+"_fit2A_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 2B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists2B, unfitHists2B, name, outDir, "_fit2B", lowBin, hiBin, doSyst, false, varsToSave, nameKeys2B);
    map<TString,MnvH1D*> scaledHists2B = {};
    scaledHists2B["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    scaledHists2B["single #pi^{0}"]=(MnvH1D*)neutPiHist->Clone();
    scaledHists2B["N#pi"]=(MnvH1D*)NPiHist->Clone();
    for(auto hists:scaledHists2B){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists2B,unfitHists2B,false,outDir+"TEST_NEW_"+name+"_fit2B_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 3A" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists3A, unfitHists3A, name, outDir, "_fit3A", lowBin, hiBin, doSyst, true, varsToSave, nameKeys3A);
    map<TString,MnvH1D*> scaledHists3A = {};
    scaledHists3A["single #pi"]=(MnvH1D*)bkg1PiHist->Clone();
    scaledHists3A["N#pi"]=(MnvH1D*)NPiHist->Clone();
    scaledHists3A["Signal"]=(MnvH1D*)sigHist->Clone();
    for(auto hists:scaledHists3A){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists3A,unfitHists3A,true,outDir+"TEST_NEW_"+name+"_fit3A_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 3B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists3B, unfitHists3B, name, outDir, "_fit3B", lowBin, hiBin, doSyst, false, varsToSave, nameKeys3B);
    map<TString,MnvH1D*> scaledHists3B = {};
    scaledHists3B["single #pi"]=(MnvH1D*)bkg1PiHist->Clone();
    scaledHists3B["N#pi"]=(MnvH1D*)NPiHist->Clone();
    for(auto hists:scaledHists3B){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists3B,unfitHists3B,false,outDir+"TEST_NEW_"+name+"_fit3B_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 4A" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists4A, unfitHists4A, name, outDir, "_fit4A", lowBin, hiBin, doSyst, true, varsToSave, nameKeys4A);
    map<TString,MnvH1D*> scaledHists4A = {};
    scaledHists4A["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    scaledHists4A["N#pi & single #pi^{0}"]=(MnvH1D*)bkgNNeutPiHist->Clone();
    scaledHists4A["Signal"]=(MnvH1D*)sigHist->Clone();
    for(auto hists:scaledHists4A){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists4A,unfitHists4A,true,outDir+"TEST_NEW_"+name+"_fit4A_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 4B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists4B, unfitHists4B, name, outDir, "_fit4B", lowBin, hiBin, doSyst, false, varsToSave, nameKeys4B);
    map<TString,MnvH1D*> scaledHists4B = {};
    scaledHists4B["single #pi^{#pm}"]=(MnvH1D*)chargePiHist->Clone();
    scaledHists4B["N#pi & single #pi^{0}"]=(MnvH1D*)bkgNNeutPiHist->Clone();
    for(auto hists:scaledHists4B){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists4B,unfitHists4B,false,outDir+"TEST_NEW_"+name+"_fit4B_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 5A" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists5A, unfitHists5A, name, outDir, "_fit5A", lowBin, hiBin, doSyst, true, varsToSave, nameKeys5A);
    map<TString,MnvH1D*> scaledHists5A = {};
    scaledHists5A["RES"]=(MnvH1D*)RESHist->Clone();
    scaledHists5A["nonRES"]=(MnvH1D*)bkgNonRESHist->Clone();
    scaledHists5A["Signal"]=(MnvH1D*)sigHist->Clone();
    for(auto hists:scaledHists5A){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists5A,unfitHists5A,true,outDir+"TEST_NEW_"+name+"_fit5A_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 5B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists5B, unfitHists5B, name, outDir, "_fit5B", lowBin, hiBin, doSyst, false, varsToSave, nameKeys5B);
    map<TString,MnvH1D*> scaledHists5B = {};
    scaledHists5B["RES"]=(MnvH1D*)RESHist->Clone();
    scaledHists5B["nonRES"]=(MnvH1D*)bkgNonRESHist->Clone();
    for(auto hists:scaledHists5B){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists5B,unfitHists5B,false,outDir+"TEST_NEW_"+name+"_fit5B_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 6A" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists6A, unfitHists6A, name, outDir, "_fit6A", lowBin, hiBin, doSyst, true, varsToSave, nameKeys6A);
    map<TString,MnvH1D*> scaledHists6A = {};
    scaledHists6A["RES"]=(MnvH1D*)RESHist->Clone();
    scaledHists6A["DIS"]=(MnvH1D*)DISHist->Clone();
    scaledHists6A["Signal"]=(MnvH1D*)sigHist->Clone();
    for(auto hists:scaledHists6A){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists6A,unfitHists6A,true,outDir+"TEST_NEW_"+name+"_fit6A_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting 6B" << endl;
    result = FitScaleFactorsAndDraw(dataHist, fitHists6B, unfitHists6B, name, outDir, "_fit6B", lowBin, hiBin, doSyst, false, varsToSave, nameKeys6B);
    map<TString,MnvH1D*> scaledHists6B = {};
    scaledHists6B["RES"]=(MnvH1D*)RESHist->Clone();
    scaledHists6B["DIS"]=(MnvH1D*)DISHist->Clone();
    for(auto hists:scaledHists6B){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    DrawFromMnvH1Ds(dataHist,scaledHists6B,unfitHists6B,false,outDir+"TEST_NEW_"+name+"_fit6B_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

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
    for (auto hist:scaledHists1A) delete hist.second;
    for (auto hist:scaledHists1B) delete hist.second;
    for (auto hist:scaledHists2A) delete hist.second;
    for (auto hist:scaledHists2B) delete hist.second;
    for (auto hist:scaledHists3A) delete hist.second;
    for (auto hist:scaledHists3B) delete hist.second;
    for (auto hist:scaledHists4A) delete hist.second;
    for (auto hist:scaledHists4B) delete hist.second;
    for (auto hist:scaledHists5A) delete hist.second;
    for (auto hist:scaledHists5B) delete hist.second;
    for (auto hist:scaledHists6A) delete hist.second;
    for (auto hist:scaledHists6B) delete hist.second;
    for (auto hist:varsToSave) delete hist.second;
    varsToSave.clear();
  }

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  mcFile->Close();
  dataFile->Close();
  outFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
