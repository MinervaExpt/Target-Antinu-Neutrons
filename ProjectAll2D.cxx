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

  //Pass an input file name to this script now
  if (argc < 2 || argc > 3) {
    cout << "Check usage..." << endl;
    return 1;
  }

  TString inFileName = argv[1];
  if (inFileName.Contains("TwoDProjected")){
    cout << "Already projected for this file. Skipping." << endl;
    return 2;
  }

  TString newFileName = "TwoDProjected_"+inFileName;

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
  TFile* newFile = new TFile(newFileName, "UPDATE");

  TDirectoryFile* dir = (TDirectoryFile*)newFile->Get("TwoD");

  bool doBinByBin = (argc==2) ? false : (atoi(argv[2]) > 0);
  
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
    MnvH2D* h2D = (MnvH2D*)newFile->Get("TwoD/"+nameObj);
    vector<TString> namePieces = BreakName("_v_",nameObj);

    int nBins = h2D->GetNbinsX();
    
    vector<TString> varPieces = BreakName("_",namePieces.at(0));
    TString varName = varPieces.at(1);

    vector<TString> tagPieces = BreakName("pT_",namePieces.at(1));    
    if (tagPieces.size()!=2) continue;
    TString nametag = tagPieces.at(1);
    
    MnvH1D* h1D = h2D->ProjectionY("Projected_"+varName+"_"+nametag,1,nBins);
    newFile->cd();
    h1D->Write();
    delete h1D;
    if (doBinByBin){
      cout << "Doing bin by bin" << endl;
      for (int iBin=1; iBin <= nBins; ++iBin){
	TString binTag = (TString)(to_string(iBin));
	MnvH1D* h1D = h2D->ProjectionY("Projected_bin_"+binTag+"_"+varName+"_"+nametag,iBin,iBin);
	newFile->cd();
	h1D->Write();
	delete h1D;
      }
    }
  }

  newFile->Close();

  return 0;
}
