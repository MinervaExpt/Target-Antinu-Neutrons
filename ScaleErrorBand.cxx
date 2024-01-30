//File: ScaleHistos.cxx
//Info: This script takes an input MC file and a scale factors file to produce a copy with the right histos scaled
//
//Usage: ScaleHistos <inFile> <scaleFile> <fitName> <outFile> <scaleSig yes if anything other than 0> optional: <enforceTag>
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
#include "PlotUtils/MnvVertErrorBand.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvVertErrorBand2D.h"
#include "PlotUtils/MnvPlotter.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

using namespace std;
using namespace PlotUtils;

MnvH1D* ScaleErrorBand(MnvH1D* h1D, std::string name, double scale){
  MnvH1D* hOut = new MnvH1D(*h1D);
  bool hasBand = hOut->HasErrorBand(name);

  if (!hasBand){
    std::vector<TH1D*> univs;
    TH1D* hUp = (TH1D*)(hOut->GetCVHistoWithStatError().Clone());
    hUp->SetDirectory(nullptr);
    TH1D* hDown = (TH1D*)(hUp->Clone());
    hUp->Scale(1+scale);
    univs.push_back(hUp);
    hDown->Scale(1-scale);
    univs.push_back(hDown);
    hOut->AddVertErrorBand(name, univs);
  }
  else{
    TH1D* hCV = (TH1D*)(hOut->GetCVHistoWithStatError().Clone());
    const auto univs = hOut->GetVertErrorBand(name)->GetHists();
    for (auto hist:univs){
      hist->Add(hCV,-1.0);
      hist->Scale(scale);
      hist->Add(hCV,1.0);
    }
  }
  return hOut;
}

MnvH2D* ScaleErrorBand(MnvH2D* h2D, std::string name, double scale){
  MnvH2D* hOut = new MnvH2D(*h2D);
  bool hasBand = hOut->HasErrorBand(name);

  if (!hasBand){
    std::vector<TH2D*> univs;
    TH2D* hUp = (TH2D*)(hOut->GetCVHistoWithStatError().Clone());
    hUp->SetDirectory(nullptr);
    TH2D* hDown = (TH2D*)(hUp->Clone());
    hUp->Scale(1+scale);
    univs.push_back(hUp);
    hDown->Scale(1-scale);
    univs.push_back(hDown);
    hOut->AddVertErrorBand(name, univs);
  }
  else{
    TH2D* hCV = (TH2D*)(hOut->GetCVHistoWithStatError().Clone());
    const auto univs = hOut->GetVertErrorBand(name)->GetHists();
    for (auto hist:univs){
      hist->Add(hCV,-1.0);
      hist->Scale(scale);
      hist->Add(hCV,1.0);
    }
  }
  return hOut;
}

int main(int argc, char* argv[]) {

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Info needed needed
  //1) Input MC File
  //1A) Output File with same POT and scaling
  //2) Scale Files
  //3) Which Scale To Use
  //4) What Gets Scaled

  //1(/1A) and 2 both just arguments
  //3) depends on chose naming for the scaling... Probably best to run different fits separately with different file names and then somehow match the naming more easily in the fitting files.
  //4) If histogram names well-matched, then this will follow naturally.
  //This all requires better development of how the naming is saved from BKG fitting instead of the 
  //Worth combining the different scalings together?

  //Pass an input file name to this script now
  if (argc != 5) {
    cout << "Check usage..." << endl;
    return 2;
  }

  TString inFileName = argv[1];
  std::string bandName = argv[2];
  double scale = atof(argv[3]);
  TString outFileName = argv[4];

  TFile* inFile = new TFile(inFileName,"READ");
  TFile* outFile = new TFile(outFileName,"RECREATE");

  TList* keyList = inFile->GetListOfKeys();
  if (!keyList){
    cout << "issue with input file" << endl;
    return 4;
  }

  TIter next(keyList);
  TKey* key;
  while ( key = (TKey*)next() ){
    TString className = (TString)key->GetClassName();
    TString nameObj = (TString)key->GetName();
    /*
    if (className == "TDirectoryFile"){
      TDirectory* newOutDir = outFile->mkdir(nameObj);
      TDirectoryFile* dirInt = (TDirectoryFile*)inFile->Get(nameObj);
      TList*  keyIntList = dirInt->GetListOfKeys();
      if(!keyIntList){
        cout << "List of keys failed to get inside second directory" << endl;
        return 20;
      }
      TIter nextKeyInt(keyIntList);
      TKey* keyInt;
      while ( keyInt = (TKey*)nextKeyInt() ){
	TString classNameInt = (TString)keyInt->GetClassName();
        TString nameObjInt = (TString)keyInt->GetName();
        if (!(classNameInt.Contains("MnvH")) || nameObjInt.Contains("MYBins")) continue;
        else if (classNameInt.Contains("MnvH2")){
	  //Could eventually be used as a way to check nuisance variables from a fit... for now just a direct copy over.
	  bool scaled = false;
	  MnvH2D* h2D = (MnvH2D*)(inFile->Get(nameObj+"/"+nameObjInt))->Clone(nameObjInt);
	  for (auto varName: varTagsMap){
	    if (scaled){
	      continue;
	    }
	    for (auto tag: varName.second){
	      if(!scaleSig && tag.Contains("sig")){
		continue;
	      }
	      if(!nameObjInt.Contains(tag) || !nameObjInt.Contains(forceTag)){
		continue;
	      }
	      cout << "Scaling: " << nameObjInt << endl;
	      TString nameOfScale = varNameMap[varName.first+tag];
	      cout << "With Scale: " << nameOfScale << endl;
	      MnvH1D* hScale = (MnvH1D*)(scaleFile->Get(nameOfScale)->Clone());
	      hScale->AddMissingErrorBandsAndFillWithCV(*h2D);
	      MnvH2D* hScale2D = Make2DX(hScale,h2D);//Assume x axis, not worth trying to specify otherwise right now. Plan is that one could change the directory to match the tag of what's scaling, but just need to get something going first.
	      h2D->Multiply(h2D,hScale2D);
	      delete hScale;
	      delete hScale2D;
	      scaled = true;
	      break;
	    }
	  }
	  newOutDir->cd();
	  h2D->Write();
	  delete h2D;
	}
	else if (classNameInt.Contains("MnvH1")){
	  bool scaled = false;
	  MnvH1D* h1D = (MnvH1D*)(inFile->Get(nameObj+"/"+nameObjInt))->Clone(nameObjInt);
	  for (auto varName: varTagsMap){
	    if (!nameObjInt.Contains(varName.first) || scaled) continue;
	    for (auto tag: varName.second){
	      if(!scaleSig && tag.Contains("sig")) continue;
	      if(!nameObjInt.Contains(tag) || !nameObjInt.Contains(forceTag)) continue;
	      cout << "Scaling: " << nameObjInt << endl;
	      TString nameOfScale = varNameMap[varName.first+tag];
	      cout << "With Scale: " << nameOfScale << endl;
	      MnvH1D* hScale = (MnvH1D*)(scaleFile->Get(nameOfScale)->Clone());
	      hScale->AddMissingErrorBandsAndFillWithCV(*h1D);
	      h1D->Multiply(h1D,hScale);
	      delete hScale;
	      scaled = true;
	      break;
	    }
	  }
	  newOutDir->cd();
	  h1D->Write();
	  delete h1D;
	}
	else {
	  cout << "HUH Inside?" << endl;
	}
      }
      }*/

    if (!(className.Contains("MnvH") || className == "TParameter<double>") || nameObj.Contains("MYBins")) continue;
    else if (className == "TParameter<double>"){
      TParameter<double>* tPar = (TParameter<double>*)(inFile->Get(nameObj))->Clone(nameObj);
      outFile->cd();
      tPar->Write();
      delete tPar;
    }
    else if (className.Contains("MnvH2")){
      MnvH2D* h2D = (MnvH2D*)(inFile->Get(nameObj))->Clone(nameObj);
      MnvH2D* h2D2 = ScaleErrorBand(h2D,bandName,scale);
      outFile->cd();
      h2D2->Write();
      delete h2D2;
      delete h2D;
    }
    else if (className.Contains("MnvH1")){
      MnvH1D* h1D = (MnvH1D*)(inFile->Get(nameObj))->Clone(nameObj);
      MnvH1D* h1D2 = ScaleErrorBand(h1D,bandName,scale);
      outFile->cd();
      h1D2->Write();
      delete h1D2;
      delete h1D;
    }
    else{
      cout << "HUH?" << endl;
    }
  }

  outFile->cd();

  cout << "Closing outFile" << endl;
  outFile->Close();

  cout << "Closing input file" << endl;
  inFile->Close();

  cout << "HEY YOU DID IT" << endl;
  return 0;
}
