//File: ScaleMCToData.cxx
//Info: This script takes an input MC file and a scale factors file to produce a copy with the right histos scaled
//
//Usage: ScaleMCToData <mcFile> <dataFile> <outFile> : optional <fixed scale>
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
#include "PlotUtils/MnvPlotter.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

using namespace std;
using namespace PlotUtils;

MnvH1D* SwapSysUniverse(MnvH1D* input, TString errorBandName, int univ){
  cout << "Working with: " << input->GetName() << endl;
  MnvH1D* ret = nullptr;
  int nUniv = input->GetVertErrorBand(errorBandName.Data())->GetHists().size();
  if (univ >= nUniv || nUniv <= 0) return ret;
  TH1D* CV = (TH1D*)(input->GetCVHistoWithStatError().Clone());
  TH1D* errBand = (TH1D*)(input->GetVertErrorBand(errorBandName.Data())->GetHist(univ)->Clone());
  if (!errBand){
    delete CV;
    delete input;

    return ret;
  }
  errBand->Divide(CV,errBand);
  for (int iBin=0; iBin <= errBand->GetNbinsX()+1; ++iBin){
    errBand->SetBinError(iBin,0);
  }
  MnvH1D* rat = new MnvH1D(*errBand);
  rat->AddMissingErrorBandsAndFillWithCV(*input);
  ret = (MnvH1D*)(input->Clone());
  ret->Divide(ret,rat);
  if (nUniv == 1){
    TH1D* newErrHist = ret->GetVertErrorBand(errorBandName.Data())->GetHist(univ);
    for (int iBin=0; iBin <= newErrHist->GetNbinsX()+1; ++iBin){
      newErrHist->SetBinContent(iBin, CV->GetBinContent(iBin));
    }
  }
  delete CV;
  delete errBand;
  delete rat;
  delete input;
  
  return ret;
}

MnvH2D* SwapSysUniverse(MnvH2D* input, TString errorBandName, int univ){
  cout << "Working with: " << input->GetName() << endl;
  MnvH2D* ret = nullptr;
  int nUniv = input->GetVertErrorBand(errorBandName.Data())->GetHists().size();
  if (univ >= nUniv || nUniv <= 0) return ret;
  TH2D* CV = (TH2D*)(input->GetCVHistoWithStatError().Clone());
  TH2D* errBand = (TH2D*)(input->GetVertErrorBand(errorBandName.Data())->GetHist(univ)->Clone());
  if (!errBand){
    delete CV;
    delete input;

    return ret;
  }
  errBand->Divide(CV,errBand);
  for (int iBinX=0; iBinX <= errBand->GetNbinsX()+1; ++iBinX){
    for (int iBinY=0; iBinY <= errBand->GetNbinsY()+1; ++iBinY){    
      errBand->SetBinError(iBinX,iBinY,0);
    }
  }
  MnvH2D* rat = new MnvH2D(*errBand);
  rat->AddMissingErrorBandsAndFillWithCV(*input);
  ret = (MnvH2D*)(input->Clone());
  ret->Divide(ret,rat);
  if (nUniv == 1){
    TH2D* newErrHist = ret->GetVertErrorBand(errorBandName.Data())->GetHist(univ);
    for (int iBinX=0; iBinX <= newErrHist->GetNbinsX()+1; ++iBinX){
      for (int iBinY=0; iBinY <= newErrHist->GetNbinsY()+1; ++iBinY){
	newErrHist->SetBinContent(iBinX, iBinY, CV->GetBinContent(iBinX, iBinY));
      }
    }
  }
  delete CV;
  delete errBand;
  delete rat;
  delete input;
  
  return ret;
}

TString MashNames(TString tag, vector<TString> pieces){
  TString name;
  for(unsigned int i=0; i<pieces.size()-1; ++i){
    name = name+pieces.at(i)+tag;
  }
  name = name+pieces.at(pieces.size()-1);
  return name;
}

vector<TString> BreakName(TString tag, TString name){
  vector<TString> namePieces;
  string search = tag.Data();
  string nameStub = name.Data();
  string token;
  size_t pos = 0;
  while ((pos = nameStub.find(search)) != string::npos){
    token = nameStub.substr(0,pos);
    namePieces.push_back(token.c_str());
    nameStub.erase(0,pos+search.length());
  }
  namePieces.push_back(nameStub.c_str());
  return namePieces;
}

bool PathExists(string path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

int main(int argc, char* argv[]) {

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc != 5) {
    cout << "Check usage..." << endl;
    return 2;
  }

  const TString inFileName = argv[1];
  const TString outFileName = argv[2];
  const TString univName = argv[3];
  const int univ = atoi(argv[4]);

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
	if (!(classNameInt.Contains("MnvH"))) continue;
	else if (classNameInt.Contains("MnvH2")){
	  MnvH2D* h2D = SwapSysUniverse((MnvH2D*)(inFile->Get(nameObj+"/"+nameObjInt)),univName,univ);
	  newOutDir->cd();
	  h2D->Write();
	  delete h2D;
	}
	else if (classNameInt.Contains("MnvH1")){
	  MnvH1D* h1D = SwapSysUniverse((MnvH1D*)(inFile->Get(nameObj+"/"+nameObjInt)),univName,univ);
	  newOutDir->cd();
	  h1D->Write();
	  delete h1D;
	}
	else {
	  cout << "HUH Inside?" << endl;
	}
      }
    }
    
    else if (!(className.Contains("MnvH") || className == "TParameter<double>") || nameObj.Contains("MYBins")) continue;
    else if (className == "TParameter<double>"){
      cout << "Working with: " << nameObj << endl;
      TParameter<double>* tPar = (TParameter<double>*)(inFile->Get(nameObj));
      outFile->cd();
      tPar->Write();
      delete tPar;
    }
    else if (className.Contains("MnvH2")){
      MnvH2D* h2D = SwapSysUniverse((MnvH2D*)(inFile->Get(nameObj)),univName,univ);
      outFile->cd();
      h2D->Write();
      delete h2D;
    }
    else if (className.Contains("MnvH1")){
      cout << "Working with: " << nameObj << endl;
      MnvH1D* h1D = SwapSysUniverse((MnvH1D*)(inFile->Get(nameObj)),univName,univ);
      outFile->cd();
      h1D->Write();
      delete h1D;
    }
    else{
      cout << "HUH?" << endl;
    }
  }

  outFile->cd();
  cout << "Closing outFile" << endl;
  outFile->Close();

  cout << "Closing input files" << endl;
  inFile->Close();
  cout << "HEY YOU DID IT" << endl;
  return 0;
}
