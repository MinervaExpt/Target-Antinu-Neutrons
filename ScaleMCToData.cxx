//File: ScaleMCToData.cxx
//Info: This script takes an input MC file and a scale factors file to produce a copy with the right histos scaled
//
//Usage: ScaleMCToData <mcFile> <dataFile> <outFile>
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
  if (argc != 4) {
    cout << "Check usage..." << endl;
    return 2;
  }

  TString mcFileName = argv[1];
  TString dataFileName = argv[2];
  TString outFileName = argv[3];

  TFile* mcFile = new TFile(mcFileName,"READ");
  TFile* dataFile = new TFile(dataFileName,"READ");
  TFile* outFile = new TFile(outFileName,"RECREATE");

  double mcPOT = ((TParameter<double>*)mcFile->Get("POTUsed"))->GetVal();
  double dataPOT = ((TParameter<double>*)dataFile->Get("POTUsed"))->GetVal();

  double scale = dataPOT/mcPOT;

  TList* keyList = mcFile->GetListOfKeys();
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
      TDirectoryFile* dirInt = (TDirectoryFile*)mcFile->Get(nameObj);
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
	  MnvH2D* h2D = (MnvH2D*)(mcFile->Get(nameObj+"/"+nameObjInt))->Clone(nameObjInt);
	  if(nameObjInt.Contains("reweightedflux")){
	    h2D->Scale(dataPOT);
	  }
	  else{
	    h2D->Scale(scale);
	  }
	  newOutDir->cd();
	  h2D->Write();
	  delete h2D;
	}
	else if (classNameInt.Contains("MnvH1")){
	  MnvH1D* h1D = (MnvH1D*)(mcFile->Get(nameObj+"/"+nameObjInt))->Clone(nameObjInt);
	  if(nameObjInt.Contains("reweightedflux")){
	    h1D->Scale(dataPOT);
	  }
	  else{
	    h1D->Scale(scale);
	  }
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
      if (!nameObj.Contains("POT")){
	TParameter<double>* tPar = (TParameter<double>*)(mcFile->Get(nameObj))->Clone(nameObj);
	outFile->cd();
	tPar->Write();
	delete tPar;
      }
      else{
	TParameter<double>* tPar = (TParameter<double>*)(dataFile->Get(nameObj))->Clone(nameObj);
	outFile->cd();
	tPar->Write();
	delete tPar;
      }
    }
    else if (className.Contains("MnvH2")){
      MnvH2D* h2D = (MnvH2D*)(mcFile->Get(nameObj))->Clone(nameObj);
      if(nameObj.Contains("reweightedflux")){
	h2D->Scale(dataPOT);
      }
      else{
	h2D->Scale(scale);
      }
      outFile->cd();
      h2D->Write();
      delete h2D;
    }
    else if (className.Contains("MnvH1")){
      MnvH1D* h1D = (MnvH1D*)(mcFile->Get(nameObj))->Clone(nameObj);
      if(nameObj.Contains("reweightedflux")){
	h1D->Scale(dataPOT);
      }
      else{
	h1D->Scale(scale);
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
  /*
  cout << "Getting POT" << endl;
  TParameter<double>* POT = (TParameter<double>*)(inFile->Get("POTUsed"))->Clone("POTUSED");
  cout << "Writing POT" << endl;
  POT->Write();
  */

  cout << "Closing outFile" << endl;
  outFile->Close();

  cout << "Closing input files" << endl;
  mcFile->Close();
  dataFile->Close();

  cout << "HEY YOU DID IT" << endl;
  return 0;
}
