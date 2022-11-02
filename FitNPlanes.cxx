//File: FitNPlanes.cxx
//Info: Tries to fit one given target material for its US plastic contamination
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
#include "Math/IFunction.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TMinuitMinimizer.h"

//Analysis includes
#include "fits/ScaleFactors.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"

using namespace std;
using namespace PlotUtils;

//vector<TString> TgtTypes = {"C","Fe","Pb","Water","USPlastic","DSPlastic","Other"};
vector<TString> TgtTypes = {"C","Fe","Pb","Water","Other"};//,"USPlastic","DSPlastic","Other"};

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

void DrawFromMnvH1Ds(MnvH1D* h_data, map<TString, MnvH1D*> hFit, map<TString, MnvH1D*> hUnfit, bool sigFit, TString nameToSave, TString whichPlastic){

  if (hFit.size()==0){
    cout << "This script should not be used to plot histograms not involved with fitting at this time." << endl;
    return;
  }

  MnvH1D* mcSum;
  MnvH1D* unfitSum;
  MnvH1D* h_sig;

  bool allFit = ((sigFit && hUnfit.size()==0) || (!sigFit && hUnfit.size()==1));
  
  //cout << "Getting signal hist." << endl;

  if (sigFit) h_sig = hFit[whichPlastic];
  else h_sig = hUnfit[whichPlastic];

  mcSum = h_sig->Clone();

  //cout << "Summing hists." << endl;

  for (auto hist:hFit){
    if (hist.first != whichPlastic) mcSum->Add(hist.second);
  }

  int iHist = 0;
  for (auto hist:hUnfit){
    if (hist.first != whichPlastic){
      mcSum->Add(hist.second);
      if (iHist!=0)unfitSum->Add(hist.second);
      else unfitSum = hist.second->Clone();
      iHist++;
    }
  }

  //cout << "Coloring hists." << endl;  
  h_sig->SetLineColor(TColor::GetColor("#999933"));
  h_sig->SetFillColor(TColor::GetColor("#999933"));

  iHist = 0;
  for (auto hist:hUnfit){
    if (hist.first != whichPlastic){
      hist.second->SetLineColor(TColor::GetColor(colors[iHist%3]));
      hist.second->SetFillColor(TColor::GetColor(colors[iHist%3]));
      iHist++;
    }
  }

  //cout << "Stacking hists." << endl;

  THStack* h = new THStack();
  //if(!allFit) h->Add((TH1D*)unfitSum->GetCVHistoWithError().Clone());
  for (auto hist:hUnfit){
    if(hist.first != whichPlastic) h->Add((TH1D*)hist.second->GetCVHistoWithError().Clone());
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
  leg->AddEntry(h_sig,whichPlastic);
  for (auto hist = hUnfit.rbegin(); hist != hUnfit.rend(); ++hist){
    if (hist->first != whichPlastic) leg->AddEntry(hist->second,hist->first);
  }
  //if (!allFit) leg->AddEntry(unfitSum,"Not fit BKGs");

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

map<TString,map<TString,MnvH1D*>> FitScaleFactorsAndDraw(MnvH1D* dataHist, map<TString, MnvH1D*> fitHistsAndNames, map<TString, MnvH1D*> unfitHistsAndNames, TString varName, TString outDir, TString fitName, int lowBin, int hiBin, bool doSyst, bool sigFit, map<TString,MnvH1D*> varsToSave, TString whichPlastic){
  map<TString,map<TString,MnvH1D*>> scaleHists = {};

  TString nameTag = fitName+"_low_"+(TString)(to_string(lowBin))+"_hi_"+(TString)(to_string(hiBin));
  TString name = varName+nameTag;

  if (PathExists((string)(outDir+name+"_postFit.pdf"))){
    cout << "Already performed fits over this range for this histo." << endl;
    cout << "If you are doing this because of updated histos, it is in your best interest to save this elsewhere or remove the old plots." << endl;
    return scaleHists;
  }

  if (!PathExists((string)(outDir+varName+fitName+"_preFit.pdf"))){
    DrawFromMnvH1Ds(dataHist,fitHistsAndNames,unfitHistsAndNames,sigFit,outDir+varName+fitName+"_preFit", whichPlastic);
  }

  cout << "Getting Data" << endl;
  TH1D* hData = (TH1D*)dataHist->GetCVHistoWithStatError().Clone();
  vector<TH1D*> fitHists = {};
  vector<TH1D*> unfitHists = {};

  for (auto hists:fitHistsAndNames){
    fitHists.push_back((TH1D*)hists.second->GetCVHistoWithStatError().Clone());
    for (auto var:varsToSave){
      if (var.first == varName) continue;
      scaleHists[var.first][hists.first] = new MnvH1D(var.first+"_fit_"+name+"_"+hists.first,"",var.second->GetNbinsX(),var.second->GetXaxis()->GetXbins()->GetArray());
    }
    scaleHists[varName][hists.first]= new MnvH1D(varName+"_fit_"+name+"_"+hists.first,"",hists.second->GetNbinsX(),hists.second->GetXaxis()->GetXbins()->GetArray());
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
    for (auto var:scaleHists){
      for (int iBin=0; iBin <= var.second[hist.first]->GetNbinsX()+1; ++iBin){
	var.second[hist.first]->SetBinContent(iBin,scaleResults[nextPar]);
	var.second[hist.first]->SetBinError(iBin,scaleErrors[nextPar]);
      }
      var.second[hist.first]->AddMissingErrorBandsAndFillWithCV(*hist.second);
    }
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

  DrawFromMnvH1Ds(dataHist,fitHistsAndNames,unfitHistsAndNames,sigFit,outDir+name+"_postFit",whichPlastic);
  
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
  if (argc != 8){
    cout << "Check usage..." << endl;
    return 2;
  }

  string MCfileName = string(argv[1]);
  string DATAfileName = string(argv[2]);
  string outDir = string(argv[3]);
  string outFileTag = string(argv[4]);
  //int tgtNum = atoi(argv[4]);
  //TString tgtName = (tgtNum != 6) ? "Tgt"+to_string(tgtNum) : "WaterTgt";
  TString material = argv[5];
  TString mat = "_"+material;
  bool doSyst = (bool)(atoi(argv[6]));
  TString USDS = argv[7];
  
  int lowBin = 0;//This is US default
  int hiBin = -1;//This is US default

  if (USDS == "US"){
    lowBin = 5;
    hiBin = 8;
  }
  else if (USDS == "DS"){
    lowBin = 2;
    hiBin = 5;
  }
  else{
    cout << "Argument for upstream v. downstream must be US or DS exactly." << endl;
    return 626;//Cause inter"Stitch"ial plastic...
  }

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

  vector<TString> namesToSave = {"pTmu","recoilE"};

  cout << "Input Data file name parsed to: " << fileNameStub << endl;
  TFile* outFile = new TFile((outDir+"fitNPlanes"+(string)(USDS)+"Results_"+outFileTag+".root").c_str(),"RECREATE");
  TFile* mcFile = new TFile(MCfileName.c_str(),"READ");
  TFile* dataFile = new TFile(DATAfileName.c_str(),"READ");

  TParameter<double>* mcPOT = (TParameter<double>*)mcFile->Get("POTUsed");
  TParameter<double>* dataPOT = (TParameter<double>*)dataFile->Get("POTUsed");

  double scale = dataPOT->GetVal()/mcPOT->GetVal();
  MnvH1D* h_sig = (MnvH1D*)(mcFile->Get(USDS+"_ByType_Plastic/vtxZ_ByTgt_"+material+"_PreRecoilCut_Plastic_selected_signal_reco"))->Clone();
  h_sig->Scale(scale);

  MnvH1D* h_1PiC_Bkg = (MnvH1D*)(mcFile->Get(USDS+"_ByType_Plastic/vtxZ_ByTgt_"+material+"_PreRecoilCut_Plastic_background_1chargePi"))->Clone();
  h_1PiC_Bkg->Scale(scale);

  MnvH1D* h_1Pi0_Bkg = (MnvH1D*)(mcFile->Get(USDS+"_ByType_Plastic/vtxZ_ByTgt_"+material+"_PreRecoilCut_Plastic_background_1neutPi"))->Clone();
  h_1Pi0_Bkg->Scale(scale);

  MnvH1D* h_NPi_Bkg = (MnvH1D*)(mcFile->Get(USDS+"_ByType_Plastic/vtxZ_ByTgt_"+material+"_PreRecoilCut_Plastic_background_NPi"))->Clone();
  h_NPi_Bkg->Scale(scale);

  MnvH1D* h_Other_Bkg = (MnvH1D*)(mcFile->Get(USDS+"_ByType_Plastic/vtxZ_ByTgt_"+material+"_PreRecoilCut_Plastic_background_Other"))->Clone();
  h_Other_Bkg->Add((MnvH1D*)(mcFile->Get(USDS+"_ByType_Plastic/vtxZ_ByTgt_"+material+"_PreRecoilCut_Plastic_background_Wrong_Nucleus")));
  h_Other_Bkg->Add((MnvH1D*)(mcFile->Get(USDS+"_ByType_Plastic/vtxZ_ByTgt_"+material+"_PreRecoilCut_Plastic_background_DSPlastic")));
  h_Other_Bkg->Add((MnvH1D*)(mcFile->Get(USDS+"_ByType_Plastic/vtxZ_ByTgt_"+material+"_PreRecoilCut_Plastic_background_USPlastic")));
  cout << "Other" << endl;
  h_Other_Bkg->Scale(scale);

  MnvH1D* h_Plastic = h_sig->Clone();
  h_Plastic->Add(h_1PiC_Bkg);
  h_Plastic->Add(h_1Pi0_Bkg);
  h_Plastic->Add(h_NPi_Bkg);
  h_Plastic->Add(h_Other_Bkg);

  MnvH1D* h_Mat_Bkg = (MnvH1D*)(mcFile->Get(USDS+"_ByType"+mat+"/vtxZ_ByTgt_"+material+"_PreRecoilCut"+mat+"_data"))->Clone();
  cout << "Mat" << endl;
  h_Mat_Bkg->Scale(scale);

  MnvH1D* h_All_Bkg = new MnvH1D("h_All_Bkg","",h_sig->GetNbinsX(),h_sig->GetXaxis()->GetXbins()->GetArray());
  h_All_Bkg->AddMissingErrorBandsAndFillWithCV(*h_sig);
  //LOOP HERE For All Materials, skip the one we're interested in.
  for (auto typeName : TgtTypes){
    if (typeName == material) continue;
    h_All_Bkg->Add((MnvH1D*)(mcFile->Get(USDS+"_ByType_"+typeName+"/vtxZ_ByTgt_"+material+"_PreRecoilCut_"+typeName+"_data")));
  }
  h_All_Bkg->Scale(scale);

  MnvH1D* h_data = (MnvH1D*)(dataFile->Get(USDS+"_ByType_Other/vtxZ_ByTgt_"+material+"_PreRecoilCut_Other_data"))->Clone();

  TString whichPlastic = "";
  if (USDS == "US") whichPlastic = "USPlastic";
  else if (USDS == "DS") whichPlastic = "DSPlastic";

  map<TString,MnvH1D*> varsToSave = {};
  for (auto nameSave:namesToSave){
    if (nameSave == "NPlanes") continue;
    TString getName = nameSave+"_Inner"+whichPlastic+"_PreRecoilCut"+mat+"_selected_signal_reco";
    cout << getName << endl;
    varsToSave[nameSave] = (MnvH1D*)(mcFile->Get(getName))->Clone();
  }

  map<TString, MnvH1D*> fitHists, unfitHists;

  fitHists[whichPlastic] = (MnvH1D*)h_Plastic->Clone();
  unfitHists[material] = (MnvH1D*)h_Mat_Bkg->Clone();
  //"if" the above line constrained by tgtNum != 6...
  //else unfitHists["Water"] = (MnvH1D*)h_Mat_Bkg->Clone();
  unfitHists["Other Materials"] = (MnvH1D*)h_All_Bkg->Clone();
  
  cout << "Fitting Vtx" << endl;
  map<TString,map<TString,MnvH1D*>> result = FitScaleFactorsAndDraw(h_data, fitHists, unfitHists, "NPlanes", outDir, "NPlanes_fit", lowBin, hiBin, doSyst, true, varsToSave, whichPlastic);

  map<TString,MnvH1D*> scaledHists = {};
  scaledHists[whichPlastic]=(MnvH1D*)h_Plastic->Clone();

  for(auto hists:scaledHists){
    hists.second->Multiply(hists.second,result["NPlanes"][hists.first]);
    for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
  }
  DrawFromMnvH1Ds(h_data,scaledHists,unfitHists,true,outDir+"TEST_NEW_NPlanes_postFit",whichPlastic);

  outFile->Write();
  
  for (auto hist:result){
    for (auto var:hist.second) delete var.second;
  }
  result.clear();
  
  delete h_data;
  delete h_sig;
  delete h_1PiC_Bkg;
  delete h_1Pi0_Bkg;
  delete h_NPi_Bkg;
  delete h_Other_Bkg;
  delete h_Mat_Bkg;
  delete h_All_Bkg;
  for (auto hist:fitHists) delete hist.second;
  for (auto hist:unfitHists) delete hist.second;
  for (auto hist:scaledHists) delete hist.second;
  for (auto hist:varsToSave) delete hist.second;
  varsToSave.clear();

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  mcFile->Close();
  dataFile->Close();
  outFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
