//File: ScaleHistos.cxx
//Info: This script takes an input MC file and a scale factors file to produce a copy with the right histos scaled
//
//Usage: ScaleHistos <inFile> <scaleFile> <fitName> <outFile>
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
  TString scaleFileName = argv[2];
  TString fitName = argv[3];
  TString outFileName = argv[4];

  TFile* inFile = new TFile(inFileName,"READ");
  TFile* scaleFile = new TFile(scaleFileName,"READ");
  TFile* outFile = new TFile(outFileName,"RECREATE");

  TList* scaleList = scaleFile->GetListOfKeys();
  if (!scaleList){
    cout << "issue with scale File" << endl;
    return 3;
  }

  TIter nextScale(scaleList);
  TKey* scale;
  map<TString,vector<TString>> varTagsMap;
  map<TString,TString> varNameMap;
  while ( scale = (TKey*)nextScale() ){
    TString scaleName = (TString)scale->GetName();
    if (!scaleName.Contains("fit_"+fitName)) continue;
    cout << scaleName << endl;
    vector<TString> fitPieces = BreakName("fit_"+fitName+"_t_",scaleName);
    vector<TString> varNamePieces = BreakName("_",fitPieces.at(0));
    TString varName = varNamePieces.at(0);
    cout << varName << endl;
    TString tagChunk = fitPieces.at(1);
    vector<TString> tags = BreakName("_t_",tagChunk);
    //TODO: Make able to handle more complicated case of not just using the fit from 500-1000 of just the background categories.
    for(auto tag: tags) varTagsMap[varName].push_back(tag);
    varNameMap[varName+tags.at(0)] = scaleName;
    cout << "" << endl;
  }

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
    if (!className.Contains("MnvH") || nameObj.Contains("MYBins")) continue;
    else if (className.Contains("MnvH2")){
      MnvH2D* h2D = (MnvH2D*)(inFile->Get(nameObj))->Clone(nameObj);
      outFile->cd();
      h2D->Write();
      delete h2D;
    }
    else if (className.Contains("MnvH1")){
      bool scaled = false;
      MnvH1D* h1D = (MnvH1D*)(inFile->Get(nameObj))->Clone(nameObj);
      for (auto varName: varTagsMap){
	if (!nameObj.Contains(varName.first) || scaled) continue;
	for (auto tag: varName.second){
	  if(!nameObj.Contains(tag)) continue;
	  cout << "Scaling: " << nameObj << endl;
	  TString nameOfScale = varNameMap[varName.first+varName.second.at(0)];
	  MnvH1D* hScale = (MnvH1D*)scaleFile->Get(nameOfScale);
	  h1D->Multiply(h1D,hScale);
	  delete hScale;
	  scaled = true;
	  break;
	}
      }
      outFile->cd();
      h1D->Write();
      delete h1D;
    }
    else{
      cout << "HUH?" << endl;
    }
  }

  outFile->cd();
  cout << "Getting POT" << endl;
  TParameter<double>* POT = (TParameter<double>*)(inFile->Get("POTUsed"))->Clone("POTUSED");
  cout << "Writing POT" << endl;
  POT->Write();

  cout << "Closing outFile" << endl;
  outFile->Close();

  cout << "Closing input files" << endl;
  inFile->Close();
  scaleFile->Close();

  cout << "HEY YOU DID IT" << endl;
  return 0;
}
