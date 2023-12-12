//File: ProjectOuterPlastic.cxx
//Info: This script takes an input MC file and a scale factors file to produce a copy with the right histos scaled
//
//Usage: ProjectOuterPlastic <inFile> <outFile> <loBin> <hiBin>
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

  const int loUSBin = 5;
  const int hiUSBin = 8;
  const int loDSBin = 2;
  const int hiDSBin = 5;

  //Pass an input file name to this script now
  if (argc != 2) {
    cout << "Check usage..." << endl;
    return 1;
  }

  TString inFileName = argv[1];
  if (inFileName.Contains("ProjectedOuterPlastics")){
    cout << "Already projected for this file. Skipping." << endl;
    return 2;
  }

  TString newFileName = "ProjectedOuterPlastics_"+inFileName;

  TFile* inFile = new TFile(inFileName, "READ");
  if (!inFile){
    cout << "Bad infile" << endl;
    return 3;
  }

  if(!inFile->Cp(newFileName)){
    cout << "Couldn't copy file" << endl;
    return 4;
  }

  inFile->Close();
  TFile* newFile = new TFile(newFileName, "Update");

  TDirectoryFile* dir = (TDirectoryFile*)newFile->Get("TwoD");

  TList* keyList = dir->GetListOfKeys();
  if (!keyList){
    cout << "issue with input file" << endl;
    return 5;
  }

  TIter next(keyList);
  TKey* key;
  while ( key = (TKey*)next() ){
    TString className = (TString)key->GetClassName();
    TString nameObj = (TString)key->GetName();

    if (!className.Contains("MnvH2D")) continue;
    if (nameObj.Contains("OuterUSPlastic")){
      MnvH2D* h2D = (MnvH2D*)newFile->Get("TwoD/"+nameObj);
      vector<TString> namePieces = BreakName("OuterUSPlastic",nameObj);
      MnvH1D* h1D = h2D->ProjectionX("pTmu_OuterUSPlastic"+namePieces.at(1),loUSBin,hiUSBin);
      newFile->cd();
      h1D->Write();
      delete h1D;
    } 
    else if (nameObj.Contains("OuterDSPlastic")){
      MnvH2D* h2D = (MnvH2D*)newFile->Get("TwoD/"+nameObj);
      vector<TString> namePieces = BreakName("OuterDSPlastic",nameObj);
      MnvH1D* h1D = h2D->ProjectionX("pTmu_OuterDSPlastic"+namePieces.at(1),loDSBin,hiDSBin);
      newFile->cd();
      h1D->Write();
      delete h1D;
    }
    else continue;
  }

  newFile->Close();

  return 0;
}
