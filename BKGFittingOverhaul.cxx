//File: BKGFittingOverhaul.cxx
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
#include "fits/ScaleFactor.h"
#include "fits/Line.h"
#include "fits/NonFit.h"
#include "fits/FitMgr.h"
#include "fits/Piecewise.h"

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
  if(!allFit) h->Add((TH1D*)(unfitSum->GetBinNormalizedCopy().GetCVHistoWithError().Clone()));
  for (auto hist:hFit){
    if(hist.first != "Signal") h->Add((TH1D*)(hist.second->GetBinNormalizedCopy().GetCVHistoWithError().Clone()));
  }
  h->Add((TH1D*)(h_sig->GetBinNormalizedCopy().GetCVHistoWithError().Clone()));

  MnvH1D* hDataClone = new MnvH1D(h_data->GetBinNormalizedCopy());
  TH1D* dataHist = (TH1D*)hDataClone->GetCVHistoWithError().Clone();
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

map<TString,MnvH1D*> PerformFit(map<TString, map<TString, vector<MnvH1D*>>> mcHistsAndFitTypes, map<TString, vector<tuple<TString, int, int>>> piecewiseBins, vector<MnvH1D*> dataHists, TString varName, TString outDir, TString fitName, int lowBin, int hiBin, bool doSyst, bool sigFit, map<TString, vector<TString>> tagsToSave)
{
  cout << "low bin: " << lowBin << endl;
  cout << "hi bin: " << hiBin << endl;
  TString nameTag = fitName+"_low_"+(TString)(to_string(lowBin))+"_hi_"+(TString)(to_string(hiBin));
  TString name = varName+nameTag;

  map<TString, MnvH1D*> scaleHists;
  MnvH1D* refHist;
  bool setRef = false;
  //return scaleHists;

  //TODO: Fill these to get initial hist plots.
  map<TString, MnvH1D*> tmpMap;
  vector<map<TString, MnvH1D*>> fitHistsAndNames(dataHists.size(),map<TString,MnvH1D*>(tmpMap));
  vector<map<TString, MnvH1D*>> unfitHistsAndNames(dataHists.size(),map<TString,MnvH1D*>(tmpMap));
  for (auto hists:mcHistsAndFitTypes){
    for (auto hist:hists.second){
      for (unsigned int iRegion=0; iRegion<hist.second.size(); ++iRegion){
	if (hists.first=="NonFit"){
	  unfitHistsAndNames.at(iRegion)[hist.first]=hist.second.at(iRegion);
	}
	else{
	  fitHistsAndNames.at(iRegion)[hist.first]=hist.second.at(iRegion);
	}
      }
    }
  }

  for (unsigned int iRegion=0; iRegion<dataHists.size(); ++iRegion){
    if (!PathExists((string)(outDir+varName+fitName+"_Region_"+(TString)(to_string(iRegion))+"_preFit.pdf"))){
      DrawFromMnvH1Ds(dataHists.at(iRegion),fitHistsAndNames.at(iRegion),unfitHistsAndNames.at(iRegion),sigFit,outDir+varName+fitName+"_Region_"+(TString)(to_string(iRegion))+"_preFit");
    }
  }

  cout << "Getting Data" << endl;
  vector<TH1D*> hDataVec = {};
  for (auto dataHist : dataHists) hDataVec.push_back((TH1D*)dataHist->GetCVHistoWithStatError().Clone());

  vector<fit::Fit*> fits;
  //  vector<TH1D*> unfitHists = {};

  for (auto hists:mcHistsAndFitTypes["ScaleFactor"]){
    TString dumpTag = "";//tag which lets you know which histos to scale later.
    for (auto tag:tagsToSave[hists.first]){
      dumpTag = dumpTag+"_t_"+tag;
    }
    vector<TH1D*> fitHists = {};
    for (auto hist:hists.second){
      if (!setRef){
	refHist = (MnvH1D*)(hist->Clone());
	setRef = true;
      }
      fitHists.push_back((TH1D*)hist->GetCVHistoWithStatError().Clone());
    }
    fits.push_back(new fit::ScaleFactor(fitHists,hists.first,lowBin,hiBin));

    if (hists.second.size()){
      scaleHists[hists.first]= new MnvH1D(varName+"_fit_"+name+dumpTag,"",hists.second.at(0)->GetNbinsX(),hists.second.at(0)->GetXaxis()->GetXbins()->GetArray());
    }
  }

  for (auto hists:mcHistsAndFitTypes["Line"]){
    TString dumpTag = "";//tag which lets you know which histos to scale later.
    for (auto tag:tagsToSave[hists.first]){
      dumpTag = dumpTag+"_t_"+tag;
    }
    vector<TH1D*> fitHists = {};
    for (auto hist:hists.second){
      if (!setRef){
	refHist = (MnvH1D*)(hist->Clone());
	setRef = true;
      }
      fitHists.push_back((TH1D*)hist->GetCVHistoWithStatError().Clone());
    }
    fits.push_back(new fit::Line(fitHists,hists.first,lowBin,hiBin));

    if (hists.second.size()){
      scaleHists[hists.first]= new MnvH1D(varName+"_fit_"+name+dumpTag,"",hists.second.at(0)->GetNbinsX(),hists.second.at(0)->GetXaxis()->GetXbins()->GetArray());
    }
  }

  for (auto hists:mcHistsAndFitTypes["Piecewise"]){
    TString dumpTag = "";//tag which lets you know which histos to scale later.
    for (auto tag:tagsToSave[hists.first]){
      dumpTag = dumpTag+"_t_"+tag;
    }
    vector<TH1D*> fitHists = {};
    for (auto hist:hists.second){
      if (!setRef){
	refHist = (MnvH1D*)(hist->Clone());
	setRef = true;
      }
      fitHists.push_back((TH1D*)hist->GetCVHistoWithStatError().Clone());
    }
    vector<fit::Fit*> fitPieces;
    for (auto piece:piecewiseBins[hists.first]){
      TString fitType = get<0>(piece);
      int pieceLoBin = get<1>(piece);
      int pieceHiBin = get<2>(piece);
      if(fitType == "ScaleFactor"){
	fitPieces.push_back(new fit::ScaleFactor(fitHists,hists.first,pieceLoBin,pieceHiBin));
      }
      else if (fitType == "Line"){	
	fitPieces.push_back(new fit::Line(fitHists,hists.first,pieceLoBin,pieceHiBin));
      }
      else{
	cout << "Unsupported piecewise fit piece passed. Aborting Fits." << endl;
	return scaleHists;
      }
    }
    sort(fitPieces.begin(), fitPieces.end(), fit::compFitStartBins);
    fits.push_back(new fit::Piecewise(fitPieces,fitHists,hists.first,lowBin,hiBin));

    if (hists.second.size()){
      scaleHists[hists.first]= new MnvH1D(varName+"_fit_"+name+dumpTag,"",hists.second.at(0)->GetNbinsX(),hists.second.at(0)->GetXaxis()->GetXbins()->GetArray());
    }
  }

  for (auto hists:mcHistsAndFitTypes["NonFit"]){
    vector<TH1D*> fitHists = {};
    for (auto hist:hists.second){
      if (!setRef){
	refHist = (MnvH1D*)(hist->Clone());
	setRef = true;
      }
      fitHists.push_back((TH1D*)hist->GetCVHistoWithStatError().Clone());
    }
    fits.push_back(new fit::NonFit(fitHists,hists.first,lowBin,hiBin));
  }

  fit::FitMgr func(fits,hDataVec);
  if (!func.CheckFit()){
    cout << "OOPS" << endl;
    return scaleHists;
  }

  auto* mini = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);

  int nextPar = 0;
  cout << "SETTING PARAMETERS" << endl;
  for (auto fit:fits){
    TString fitName = fit->FitName();
    if (!fitName.Contains("NonFit")){
      for (int iPar=0; iPar < fit->NDim(); ++iPar){
	TString var = fitName+"_"+to_string(iPar);
	mini->SetVariable(nextPar+iPar,(string)var.Data(),1.0,1.0);
      }
      nextPar += fit->NDim();
    }
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

  cout << "Moving ON" << endl;

  const double* scaleResults = mini->X();
  const double* scaleErrors = mini->Errors();
  map<TString, double> scaleByName;
  nextPar=0;
  for (auto fit:fits){
    TString fitName = fit->FitName();
    if (!fitName.Contains("NonFit")){
      for (int iBin=0; iBin <= scaleHists[fitName]->GetNbinsX()+1; ++iBin){
	int whichBin = iBin;
	int fitLowBin = fit->GetFirstFitBin();
	int fitHiBin = fit->GetLastFitBin();
	if (iBin < fitLowBin) whichBin = fitLowBin;
	else if (iBin > fitHiBin) whichBin = fitHiBin;
	double fitVal = fit->GetFitVal(scaleResults,nextPar,whichBin);
	double fitErr = fit->GetFitErr(scaleResults,scaleErrors,nextPar,whichBin);
	scaleHists[fitName]->SetBinContent(iBin,fitVal);
	scaleHists[fitName]->SetBinError(iBin,fitErr);
      }
      cout << fitName << endl;
      scaleHists[fitName]->AddMissingErrorBandsAndFillWithCV(*refHist);
      nextPar += fit->NDim();
    }
  }

  if (doSyst){
    const auto errorBandNames = refHist->GetErrorBandNames();
    for (const auto& bandName: errorBandNames){
      cout << bandName << endl;
      const auto univs = refHist->GetVertErrorBand(bandName)->GetHists();
      for (size_t whichUniv=0; whichUniv < univs.size(); ++ whichUniv){
	cout << "" << endl;
	cout << "Fitting Universe " << whichUniv << " in " << bandName << " error band." << endl;
	
	vector<fit::Fit*> fitsUniv;
	//  vector<TH1D*> unfitHists = {};

	for (auto hists:mcHistsAndFitTypes["ScaleFactor"]){
	  vector<TH1D*> fitHists = {};
	  for (auto hist:hists.second){
	    fitHists.push_back((TH1D*)hist->GetVertErrorBand(bandName)->GetHist(whichUniv)->Clone());
	  }
	  fitsUniv.push_back(new fit::ScaleFactor(fitHists,hists.first,lowBin,hiBin));
	}
	
	for (auto hists:mcHistsAndFitTypes["Line"]){
	  vector<TH1D*> fitHists = {};
	  for (auto hist:hists.second){
	    fitHists.push_back((TH1D*)hist->GetVertErrorBand(bandName)->GetHist(whichUniv)->Clone());
	  }
	  fitsUniv.push_back(new fit::Line(fitHists,hists.first,lowBin,hiBin));
	}

	for (auto hists:mcHistsAndFitTypes["Piecewise"]){
	  vector<TH1D*> fitHists = {};
	  for (auto hist:hists.second){
	    fitHists.push_back((TH1D*)hist->GetVertErrorBand(bandName)->GetHist(whichUniv)->Clone());
	  }
	  vector<fit::Fit*> fitPieces;
	  for (auto piece:piecewiseBins[hists.first]){
	    TString fitType = get<0>(piece);
	    int pieceLoBin = get<1>(piece);
	    int pieceHiBin = get<2>(piece);
	    if(fitType == "ScaleFactor"){
	      fitPieces.push_back(new fit::ScaleFactor(fitHists,hists.first,pieceLoBin,pieceHiBin));
	    }
	    else if (fitType == "Line"){	
	      fitPieces.push_back(new fit::Line(fitHists,hists.first,pieceLoBin,pieceHiBin));
	    }
	    else{
	      cout << "Unsupported piecewise fit piece passed inside universe loop. Aborting Fits." << endl;
	      return scaleHists;
	    }
	  }
	  sort(fitPieces.begin(), fitPieces.end(), fit::compFitStartBins);
	  fitsUniv.push_back(new fit::Piecewise(fitPieces,fitHists,hists.first,lowBin,hiBin));
	}
	
	for (auto hists:mcHistsAndFitTypes["NonFit"]){
	  vector<TH1D*> fitHists = {};
	  for (auto hist:hists.second){
	    fitHists.push_back((TH1D*)hist->GetVertErrorBand(bandName)->GetHist(whichUniv)->Clone());
	  }
	  fitsUniv.push_back(new fit::NonFit(fitHists,hists.first,lowBin,hiBin));
	}
	
	fit::FitMgr funcUniv(fitsUniv,hDataVec);
	
	auto* miniUniv = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);
	
	int parNext = 0;
	for (auto fit:fitsUniv){
	  TString fitName = fit->FitName();
	  if (!fitName.Contains("NonFit")){
	    for (int iPar=0; iPar < fit->NDim(); ++iPar){
	      TString var = fitName+"_"+to_string(iPar);
	      miniUniv->SetVariable(parNext+iPar,(string)var.Data(),1.0,1.0);
	    }
	    parNext += fit->NDim();
	  }
	}
	
	if (parNext != funcUniv.NDim()){
	  cout << "The number of parameters was unexpected for some reason..." << endl;
	  return scaleHists;
	}

	cout << "Setting Function." << endl;
	
	miniUniv->SetFunction(funcUniv);
	
	cout << "MINIMIZING Function." << endl;

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
	for (auto fit:fitsUniv){
	  TString fitName = fit->FitName();
	  if (!fitName.Contains("NonFit")){
	    for (int iBin=0; iBin <= scaleHists[fitName]->GetVertErrorBand(bandName)->GetHist(whichUniv)->GetNbinsX()+1; ++iBin){
	      int whichBin = iBin;
	      int fitLowBin = fit->GetFirstFitBin();
	      int fitHiBin = fit->GetLastFitBin();
	      if (iBin < fitLowBin) whichBin = fitLowBin;
	      else if (iBin > fitHiBin) whichBin = fitHiBin;
	      double fitValUniv = fit->GetFitVal(scaleResultsUniv,parNext,whichBin);
	      scaleHists[fitName]->GetVertErrorBand(bandName)->GetHist(whichUniv)->SetBinContent(iBin,fitValUniv);
	    }
	    //scaleFactorHists[hist.first]->Fill(quelUniv,scaleResultsUniv[parNext]);
	    parNext+= fit->NDim();
	  }
	}
      }
    }
  }

  return scaleHists;
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

  vector<fit::Fit*> fits;
  for (auto hist:fitHists){
    std::vector<TH1D*> hists(1,hist);
    fits.push_back(new fit::ScaleFactor(hists,lowBin,hiBin));
  }

  for (auto hist:unfitHists){
    std::vector<TH1D*> hists(1,hist);
    fits.push_back(new fit::NonFit(hists,lowBin,hiBin));
  }

  std::vector<TH1D*> dataHists(1,hData);

  fit::FitMgr func(fits,dataHists);

  auto* mini = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);

  int nextPar = 0;
  int nextFit = 0;
  for (auto hist:fitHistsAndNames){
    for (int iPar=0; iPar < fits.at(nextFit)->NDim(); ++iPar){
      TString var = hist.first+"_"+to_string(iPar);
      mini->SetVariable(nextPar+iPar,(string)var.Data(),1.0,1.0);
    }
    nextPar += fits.at(nextFit)->NDim();
    ++nextFit;
  }

  if (nextPar != func.NDim()){
    cout << "The number of parameters was unexpected for some reason..." << endl;
    return scaleHists;
  }

  mini->SetFunction(func);

  /*
  cout << "CHECKING SOMETHING" << endl;
  cout << func.NDim() << endl;
  */

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
  nextFit=0;
  nextPar=0;
  for (auto hist:fitHistsAndNames){
    cout << "HI" << endl;
    scaleByName[hist.first]=scaleResults[nextPar]; //TODO: Modify this such that it affects later details correctly.
    for (auto var:scaleHists){
      cout << "HI 2" << endl;
      for (int iBin=0; iBin <= var.second[hist.first]->GetNbinsX()+1; ++iBin){
	cout << "HI 3" << endl;
	int whichBin = iBin;
	if (iBin < lowBin) whichBin = lowBin;
	else if (iBin > hiBin) whichBin = hiBin;
	double fitVal = fits.at(nextFit)->GetFitVal(scaleResults,nextPar,whichBin);
	double fitErr = fits.at(nextFit)->GetFitErr(scaleResults,scaleErrors,nextPar,whichBin);
	var.second[hist.first]->SetBinContent(iBin,fitVal);
	var.second[hist.first]->SetBinError(iBin,fitErr);
      }
      cout << "HI 4" << endl;
      cout << hist.first << endl;
      var.second[hist.first]->AddMissingErrorBandsAndFillWithCV(*hist.second);
    }
    cout << "HI 5" << endl;
    nextPar += fits.at(nextFit)->NDim();
    ++nextFit;
  }

  /*
  cout << "Checking the grabbing of scale factors." << endl;
  for (auto scale: scaleByName){
    cout << scale.first << " scale: " << scale.second << endl;
  }
  */

  //cout << "Trying to draw pre-scaling." << endl;

  //cout << "Scaling the entire MnvH1D first." << endl;
  /* Removed. Think I don't need to scale before fitting universes. Though I do lose the check I had before where I could see the effect without any error modificaitons, but I don't need that since the full scaling is more correct.
  for (auto hist:fitHistsAndNames){
    hist.second->Multiply(scaleHists[hist.first]);
  }
  */

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

	vector<fit::Fit*> fitsUniv;
	for (auto hist:fitHistsUniv){
	  std::vector<TH1D*> hists(1,hist);
	  fitsUniv.push_back(new fit::ScaleFactor(hists,lowBin,hiBin));
	}

	for (auto hist:unfitHistsUniv){
	  std::vector<TH1D*> hists(1,hist);
	  fitsUniv.push_back(new fit::NonFit(hists,lowBin,hiBin));
	}

	fit::FitMgr funcUniv(fitsUniv,dataHists);
	
	auto* miniUniv = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);
	
	int parNext = 0;
	int fitNext = 0;
	for (auto hist:fitHistsAndNames){
	  for (int iPar=0; iPar < fitsUniv.at(fitNext)->NDim(); ++iPar){
	    TString var = hist.first+"_"+to_string(iPar);
	    miniUniv->SetVariable(parNext+iPar,(string)var.Data(),1.0,1.0);
	  }
	  parNext += fitsUniv.at(fitNext)->NDim();
	  ++fitNext;
	}
	
	if (parNext != funcUniv.NDim()){
	  cout << "The number of parameters was unexpected for some reason..." << endl;
	  return scaleHists;
	}

	cout << "Setting Function." << endl;
	
	miniUniv->SetFunction(funcUniv);
	
	cout << "MINIMIZING Function." << endl;

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
	fitNext=0;
	for (auto hist:fitHistsAndNames){
	  //hist.second->GetVertErrorBand(bandName)->GetHist(whichUniv)->Scale(scaleResultsUniv[parNext]);	
	  for (auto var:scaleHists){
	    for (int iBin=0; iBin <= var.second[hist.first]->GetVertErrorBand(bandName)->GetHist(whichUniv)->GetNbinsX()+1; ++iBin){
	      int whichBin = iBin;
	      if (iBin < lowBin) whichBin = lowBin;
	      else if (iBin > hiBin) whichBin = hiBin;
	      double fitValUniv = fitsUniv.at(fitNext)->GetFitVal(scaleResultsUniv,parNext,whichBin);
	      var.second[hist.first]->GetVertErrorBand(bandName)->GetHist(whichUniv)->SetBinContent(iBin,fitValUniv);
	    }
	    //var.second[hist.first]->GetVertErrorBand(bandName)->GetHist(whichUniv)->Scale(scaleResultsUniv[parNext]);
	  }
	  //scaleFactorHists[hist.first]->Fill(quelUniv,scaleResultsUniv[parNext]);
	  parNext+= fitsUniv.at(fitNext)->NDim();
	  ++fitNext;
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
  
  //vector<TString> namesToSave = {"pTmu","recoilE"};
  vector<TString> namesToSave = {varName};//, "NPlanes"};
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

  if (varName.Contains("pTmu")){
    if (lowBin >= 9){
      cout << "Not a valid fitting range for pTmu. Maximum bin is 14." << endl;
      return 100;
    }
    hiBin = min(9,hiBin);
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
    TString grabName, grabName_noTag;
    if (tag.Contains("_Tgt")){
      TObjArray* tagArr = tag.Tokenize("Tgt");
      grabName="ByTgt_Tgt"+((TObjString*)(tagArr->At(tagArr->GetEntries()-1)))->String()+"/"+name;
      grabName_noTag = "ByTgt_Tgt"+((TObjString*)(tagArr->At(tagArr->GetEntries()-1)))->String()+"/"+varName;
    }
    else{
      grabName = name;
      grabName_noTag = varName;
    }

    cout << "Performing Fitting and Scaling for: " << name << endl;
    cout << "Grabbing: " << grabName << endl;
    cout << "TESTING NO TAG: " << grabName_noTag << endl;


    //TODO:
    // 
    double coarseVtxBinsUS[11];//Modified to be a number of planes instead of vertex Z values.                                             
    for (int iBin=0; iBin<11;++iBin) coarseVtxBinsUS[iBin] = (double)iBin-10.5;

    double coarseVtxBinsDS[11];//Modified to be a number of planes instead of vertex Z values.                                             
    for (int iBin=11; iBin<22;++iBin) coarseVtxBinsDS[iBin-11] = (double)iBin-10.5;

    double binsNPlanes[22];
    for (unsigned int i=0; i < 22;++i) binsNPlanes[i]=1.0*i-10.5;

    map<TString,MnvH1D*> varsToSave = {};
    for (auto nameSave:namesToSave){
      //TODO:
      // 
      if (nameSave == "NPlanes"){
	MnvH1D* hNPlanesUS = new MnvH1D("vtxZ_ByTgt_US"+tag,"",10,coarseVtxBinsUS);
	varsToSave[nameSave+"US"+tag] = hNPlanesUS->Clone();
	MnvH1D* hNPlanesDS = new MnvH1D(nameSave+"DS"+tag,"",10,coarseVtxBinsDS);
	varsToSave[nameSave+"DS"+tag] = hNPlanesDS->Clone();
      }
      //
      else if (tag.Contains("_Tgt")){
	TObjArray* tagArr = tag.Tokenize("Tgt");
	TString grabNameSave = "ByTgt_Tgt"+((TObjString*)(tagArr->At(tagArr->GetEntries()-1)))->String()+"/"+nameSave+tag;
	varsToSave[nameSave+tag] = (MnvH1D*)(mcFile->Get(grabNameSave+"_selected_signal_reco"))->Clone();
      }
      else varsToSave[nameSave+tag] = (MnvH1D*)(mcFile->Get(nameSave+tag+"_selected_signal_reco"))->Clone();
    }

    MnvH1D* dataHist = (MnvH1D*)(dataFile->Get(grabName+"_data"))->Clone();
    MnvH1D* dataHist_noTag = (MnvH1D*)(dataFile->Get(grabName_noTag+"_data"))->Clone();
    vector<MnvH1D*> dataHists;
    dataHists.push_back((MnvH1D*)dataHist->Clone());
    dataHists.push_back((MnvH1D*)dataHist_noTag->Clone());

    MnvH1D* sigHist = (MnvH1D*)(mcFile->Get(grabName+"_selected_signal_reco"))->Clone();
    MnvH1D* sigHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_selected_signal_reco"))->Clone();
    sigHist->Scale(POTscale);
    sigHist_noTag->Scale(POTscale);

    //
    MnvH1D* chargePiHist = (MnvH1D*)(mcFile->Get(grabName+"_background_1chargePi"))->Clone();
    MnvH1D* chargePiHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_1chargePi"))->Clone();
    chargePiHist->Scale(POTscale);
    chargePiHist_noTag->Scale(POTscale);
    MnvH1D* neutPiHist = (MnvH1D*)(mcFile->Get(grabName+"_background_1neutPi"))->Clone();
    MnvH1D* neutPiHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_1neutPi"))->Clone();
    neutPiHist->Scale(POTscale);
    neutPiHist_noTag->Scale(POTscale);
    MnvH1D* NPiHist = (MnvH1D*)(mcFile->Get(grabName+"_background_NPi"))->Clone();
    MnvH1D* NPiHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_NPi"))->Clone();
    NPiHist->Scale(POTscale);
    NPiHist_noTag->Scale(POTscale);
    MnvH1D* otherHist = (MnvH1D*)(mcFile->Get(grabName+"_background_Other"))->Clone();
    MnvH1D* otherHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_Other"))->Clone();
    otherHist->Scale(POTscale);
    otherHist_noTag->Scale(POTscale);

    MnvH1D* wrongNuclHist = (MnvH1D*)(mcFile->Get(grabName+"_background_Wrong_Nucleus"))->Clone();
    MnvH1D* wrongNuclHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_Wrong_Nucleus"))->Clone();
    wrongNuclHist->Scale(POTscale);
    wrongNuclHist_noTag->Scale(POTscale);

    MnvH1D* USHist = (MnvH1D*)(mcFile->Get(grabName+"_background_USPlastic"))->Clone();
    MnvH1D* USHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_USPlastic"))->Clone();
    USHist->Scale(POTscale);
    USHist_noTag->Scale(POTscale);

    MnvH1D* DSHist = (MnvH1D*)(mcFile->Get(grabName+"_background_DSPlastic"))->Clone();
    MnvH1D* DSHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_DSPlastic"))->Clone();
    DSHist->Scale(POTscale);
    DSHist_noTag->Scale(POTscale);

    MnvH1D* bkgNNeutPiHist = neutPiHist->Clone();
    MnvH1D* bkgNNeutPiHist_noTag = neutPiHist_noTag->Clone();
    bkgNNeutPiHist->Add(NPiHist);
    bkgNNeutPiHist_noTag->Add(NPiHist_noTag);

    MnvH1D* bkg1PiHist = neutPiHist->Clone();
    MnvH1D* bkg1PiHist_noTag = neutPiHist_noTag->Clone();
    bkg1PiHist->Add(chargePiHist);
    bkg1PiHist_noTag->Add(chargePiHist_noTag);

    //
    MnvH1D* QEHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_QE"))->Clone();
    MnvH1D* QEHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_bkg_IntType_QE"))->Clone();
    QEHist->Scale(POTscale);
    QEHist_noTag->Scale(POTscale);
    MnvH1D* RESHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_RES"))->Clone();
    MnvH1D* RESHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_bkg_IntType_RES"))->Clone();
    RESHist->Scale(POTscale);
    RESHist_noTag->Scale(POTscale);
    MnvH1D* DISHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_DIS"))->Clone();
    MnvH1D* DISHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_bkg_IntType_DIS"))->Clone();
    DISHist->Scale(POTscale);
    DISHist_noTag->Scale(POTscale);
    MnvH1D* MECHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_2p2h"))->Clone();
    MnvH1D* MECHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_bkg_IntType_2p2h"))->Clone();
    MECHist->Scale(POTscale);
    MECHist_noTag->Scale(POTscale);
    MnvH1D* OtherIntTypeHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_Other"))->Clone();
    MnvH1D* OtherIntTypeHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_bkg_IntType_Other"))->Clone();
    OtherIntTypeHist->Scale(POTscale);
    OtherIntTypeHist_noTag->Scale(POTscale);

    MnvH1D* bkgNonRESHist = DISHist->Clone();
    bkgNonRESHist->Add(QEHist);
    bkgNonRESHist->Add(MECHist);
    bkgNonRESHist->Add(OtherIntTypeHist);
    MnvH1D* bkgNonRESHist_noTag = DISHist_noTag->Clone();
    bkgNonRESHist_noTag->Add(QEHist_noTag);
    bkgNonRESHist_noTag->Add(MECHist_noTag);
    bkgNonRESHist_noTag->Add(OtherIntTypeHist_noTag);

    MnvH1D* bkgTotHist = bkgNonRESHist->Clone();
    bkgTotHist->Add(RESHist);
    bkgTotHist->Add(wrongNuclHist);
    MnvH1D* bkgTotHist_noTag = bkgNonRESHist_noTag->Clone();
    bkgTotHist_noTag->Add(RESHist_noTag);
    bkgTotHist_noTag->Add(wrongNuclHist_noTag);

    map<TString, MnvH1D*> fitHists1A, unfitHists1A;
    map<TString, map<TString, vector<MnvH1D*>>> fitTEST1A;
    map<TString, vector<tuple<TString, int, int>>> fitPieces1A;

    map<TString, vector<TString>> nameKeys1A;
    map<TString, vector<TString>> nameKeysTEST1A;

    fitHists1A["BKG"]=(MnvH1D*)bkgTotHist->Clone();
    fitHists1A["Signal"]=(MnvH1D*)sigHist->Clone();
    unfitHists1A["USPlastic"]=(MnvH1D*)USHist->Clone();
    unfitHists1A["DSPlastic"]=(MnvH1D*)DSHist->Clone();
    //unfitHists1A["WrongNucleus"]=(MnvH1D*)wrongNuclHist->Clone();
    nameKeys1A["BKG"]={"1chargePi","1neutPi","NPi","Other","Wrong_Nucleus"};
    nameKeys1A["Signal"]={"sig","signal"};

    fitTEST1A["Piecewise"]["Signal"].push_back((MnvH1D*)sigHist->Clone());
    fitTEST1A["Piecewise"]["Signal"].push_back((MnvH1D*)sigHist_noTag->Clone());
    fitTEST1A["Piecewise"]["BKG"].push_back((MnvH1D*)bkgTotHist->Clone());
    fitTEST1A["Piecewise"]["BKG"].push_back((MnvH1D*)bkgTotHist_noTag->Clone());
    fitTEST1A["NonFit"]["USPlastic"].push_back((MnvH1D*)USHist->Clone());
    fitTEST1A["NonFit"]["USPlastic"].push_back((MnvH1D*)USHist_noTag->Clone());
    fitTEST1A["NonFit"]["DSPlastic"].push_back((MnvH1D*)DSHist->Clone());
    fitTEST1A["NonFit"]["DSPlastic"].push_back((MnvH1D*)DSHist_noTag->Clone());

    fitPieces1A["Signal"].push_back(make_tuple("Line",5,999));
    fitPieces1A["Signal"].push_back(make_tuple("ScaleFactor",3,4));
    fitPieces1A["Signal"].push_back(make_tuple("ScaleFactor",1,2));

    //fitPieces1A["BKG"].push_back(make_tuple("ScaleFactor",11,999));
    fitPieces1A["BKG"].push_back(make_tuple("Line",6,999));
    fitPieces1A["BKG"].push_back(make_tuple("Line",1,6));

    nameKeysTEST1A["BKG"]=nameKeys1A["BKG"];
    nameKeysTEST1A["Signal"]=nameKeys1A["Signal"];

    cout << "Fitting 1A" << endl;
    map<TString,map<TString,MnvH1D*>> result;// = FitScaleFactorsAndDraw(dataHist, fitHists1A, unfitHists1A, name, outDir, "_fit1A", lowBin, hiBin, doSyst, true, varsToSave, nameKeys1A);
    map<TString,MnvH1D*> scaledHists1A = {};
    //scaledHists1A["BKG"]=(MnvH1D*)bkgTotHist->Clone();
    //scaledHists1A["Signal"]=(MnvH1D*)sigHist->Clone();
    for(auto hists:scaledHists1A){
      hists.second->Multiply(hists.second,result[name][hists.first]);
      for (auto var:result)var.second[hists.first]->SetDirectory(outFile);
    }
    //DrawFromMnvH1Ds(dataHist,scaledHists1A,unfitHists1A,true,outDir+"TEST_NEW_"+name+"_fit1A_postFit");
    //cout << "Result Has Size: " << result << endl;
    cout << "" << endl;

    outFile->Write();

    for (auto hist:result){
      for (auto var:hist.second) delete var.second;
    }
    result.clear();

    cout << "Fitting NEW TEST" << endl;
    map<TString, MnvH1D*> TESTresult = PerformFit(fitTEST1A, fitPieces1A, dataHists, name, outDir, "_fitTEST", lowBin, hiBin, doSyst, true, nameKeysTEST1A);

    map<TString, MnvH1D*> tmpMap;
    vector<map<TString, MnvH1D*>> scaledHistsAndNamesTEST(dataHists.size(),map<TString,MnvH1D*>(tmpMap));
    vector<map<TString, MnvH1D*>> unfitHistsAndNamesTEST(dataHists.size(),map<TString,MnvH1D*>(tmpMap));
    for (auto hists:fitTEST1A){
      for (auto hist:hists.second){
	for (unsigned int iRegion=0; iRegion<hist.second.size(); ++iRegion){
	  if (hists.first=="NonFit"){
	    unfitHistsAndNamesTEST.at(iRegion)[hist.first]=(MnvH1D*)hist.second.at(iRegion)->Clone();
	  }
	  else{
	    scaledHistsAndNamesTEST.at(iRegion)[hist.first]=(MnvH1D*)hist.second.at(iRegion)->Clone();
	    scaledHistsAndNamesTEST.at(iRegion)[hist.first]->Multiply(scaledHistsAndNamesTEST.at(iRegion)[hist.first], TESTresult[hist.first]);
	  }
	}
      }
    }

    for (unsigned int iRegion=0; iRegion<dataHists.size(); ++iRegion){
      if (!PathExists((string)(outDir+name+"_fitTEST_Region_"+(TString)(to_string(iRegion))+"_postFit.pdf"))){
	DrawFromMnvH1Ds(dataHists.at(iRegion),scaledHistsAndNamesTEST.at(iRegion),unfitHistsAndNamesTEST.at(iRegion),true,outDir+name+"_fitTEST_Region_"+(TString)(to_string(iRegion))+"_postFit");
      }
    }

    for (auto hist:TESTresult) hist.second->SetDirectory(outFile);    
    outFile->Write();
    for (auto hist:TESTresult) delete hist.second;
    TESTresult.clear();
 
    delete dataHist;
    delete sigHist;
    delete chargePiHist;
    delete neutPiHist;
    delete NPiHist;
    delete otherHist;
    delete USHist;
    delete DSHist;
    delete wrongNuclHist;
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
    for (auto hist:scaledHists1A) delete hist.second;
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
