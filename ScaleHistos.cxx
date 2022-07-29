//File: CombineTgts.cxx
//Info: This is a script to add up all the histograms across various files so that the materials in the targets are combined.
//
//Usage: CombineTgts <data/MC> <Material> <outDir> <Tgt1 FileName>
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

  bool isMC = ((TString)(argv[1]) == "MC") ? true : false;
  TString material = argv[2];
  TString matName = (material != "Water")? "_"+material : "";
  string outDir=string(argv[3]);
  string baseFileName=argv[4];

  if (material != "Fe" && material != "C" && material != "Pb" && material != "Water"){
    cout << "Can't handle that material yet." << endl;
    return -999;
  }

  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory doesn't exist. Exiting" << endl;
    return 3;
  }

  string rootExt = ".root";
  string slash = "/";
  string token;
  string fileNameStub = baseFileName;
  size_t pos=0;

  if (isMC){
    if ((pos = baseFileName.find("MC")) == string::npos){
      cout << "File is Not MC, and should be." << endl;
      return 4;
    }
  }
  else if ((pos = baseFileName.find("Data")) == string::npos){
    cout << "File is Not Data, and should be." << endl;
    return 4;
  }

  pos = 0;
  //cout << sigNameStub << endl;
  while ((pos = fileNameStub.find(slash)) != string::npos){
    //cout << sigNameStub << endl;
    token = fileNameStub.substr(0, pos);
    //cout << token << endl;
    fileNameStub.erase(0, pos+slash.length());
  }
  //cout << sigNameStub << endl;
  if ((pos=fileNameStub.find(rootExt)) == string::npos){
    cout << "Input need be .root file." << endl;
    return 5;
  }

  cout << "Input file name parsed to: " << fileNameStub << endl;

  TString tgtFind = "Tgt1";

  vector<TString> outFileNamePieces = BreakName(tgtFind,fileNameStub);
  vector<TString> fileNamePieces = BreakName(tgtFind,baseFileName);
  
  cout << "Input Name reclaimed as: " << MashNames("Tgt1",fileNamePieces) << endl;

  TString outName = MashNames(material,outFileNamePieces);
  TFile* outFile = new TFile(outDir+outName,"RECREATE");

  std::map<TString,TFile*> files;
  for (auto tgt : tgtsByMat[material]){
    TString tgtName = MashNames(tgt,fileNamePieces);
    files[tgt] = new TFile(tgtName,"READ");
  }

  TString firstName = tgtsByMat[material].at(0);
  TString dirBase = "ByTgt_";

  TDirectoryFile* dir = (TDirectoryFile*)files[firstName]->Get(dirBase+firstName+matName);
  if (!dir){
    cout << "Could not find material in the first Tgt. Exiting Now" << endl;
    return 18;
  }
 
  TList* keyList = dir->GetListOfKeys();
  if(!keyList){
    cout << "List of keys failed to get inside directory." << endl;
    return 19;
  }

  TIter nextKey(keyList);
  TKey* key;
  while ( key = (TKey*)nextKey() ){
    TString className = (TString)key->GetClassName();
    TString nameObj = (TString)key->GetName();
    cout << "AT Object: " << nameObj << endl;

    if (className == "TDirectoryFile"){
      TDirectory* newOutDir = outFile->mkdir(nameObj);
      TDirectoryFile* dirInt = (TDirectoryFile*)files[firstName]->Get(dirBase+firstName+matName+"/"+nameObj);
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
	cout << "AT Object: " << nameObjInt << endl;
	if (!classNameInt.Contains("PlotUtils::MnvH1")) continue;
	cout << "Getting: " << nameObjInt << endl;
	cout << "" << endl;
	vector<TString> nameObjIntPieces = BreakName(firstName+matName,nameObjInt);
	TString matNameObjInt = MashNames(material,nameObjIntPieces);
	
	if (nameObj.Contains("US_ByType")){
	  MnvH1D* h1D = nullptr;
	  bool first = true;
	  for (auto tgt : UStgtsByMat[material]){
	    cout << "US of Tgt: " << tgt << endl;
	    if (first){
	      TString newNameObjInt = MashNames(tgt+matName,nameObjIntPieces);
	      h1D = (MnvH1D*)(files[tgt]->Get(dirBase+tgt+matName+"/"+nameObj+"/"+newNameObjInt))->Clone(matNameObjInt);
	      first = false;
	    }
	    else{
	      TString newNameObjInt = MashNames(tgt+matName,nameObjIntPieces);
	      h1D->Add((MnvH1D*)(files[tgt]->Get(dirBase+tgt+matName+"/"+nameObj+"/"+newNameObjInt)));
	    }
	  }
	  newOutDir->cd();
	  h1D->Write();
	  delete h1D;
	}
	
	if (nameObj.Contains("DS_ByType")){
	  MnvH1D* h1D = nullptr;
	  bool first = true;
	  for (auto tgt:DStgtsByMat[material]){
	    cout << "DS of Tgt: " << tgt << endl;
	    if (first){
	      TString newNameObjInt = MashNames(tgt+matName,nameObjIntPieces);
	      h1D = (MnvH1D*)(files[tgt]->Get(dirBase+tgt+matName+"/"+nameObj+"/"+newNameObjInt))->Clone(matNameObjInt);
	      first = false;
	    }
	    else{
	      TString newNameObjInt = MashNames(tgt+matName,nameObjIntPieces);
	      h1D->Add((MnvH1D*)(files[tgt]->Get(dirBase+tgt+matName+"/"+nameObj+"/"+newNameObjInt)));
	    }
	  }
	  newOutDir->cd();
	  h1D->Write();
	  delete h1D;
	}
      }
      //continue;
    }

    else if (!className.Contains("PlotUtils::MnvH")) continue;
    cout << "Getting: " << nameObj << endl;
    cout << "" << endl;
    vector<TString> nameObjPieces = BreakName(firstName+matName,nameObj);
    TString matNameObj = MashNames(material,nameObjPieces);
    if (className.Contains("2")){
      MnvH2D* h2D = (MnvH2D*)(files[firstName]->Get(dirBase+firstName+matName+"/"+nameObj))->Clone(matNameObj);
      //cout << "First Histo has Integral: " << h2D->Integral() << endl;
      //cout << "" << endl;
      for(auto file : files){
	if (file.first == firstName) continue;
	//cout << "Getting: " << file.first << endl;
	TString newNameObj = MashNames(file.first+matName,nameObjPieces);
	h2D->Add((MnvH2D*)(file.second->Get(dirBase+file.first+matName+"/"+newNameObj)));
      }
      //h2D->SetDirectory(outFile);
      outFile->cd();
      h2D->Write();
      delete h2D;
    }
    else{
      MnvH1D* h1D = (MnvH1D*)(files[firstName]->Get(dirBase+firstName+matName+"/"+nameObj))->Clone(matNameObj);
      for(auto file : files){
	if (file.first == firstName) continue;
	TString newNameObj = MashNames(file.first+matName,nameObjPieces);
	h1D->Add((MnvH1D*)(file.second->Get(dirBase+file.first+matName+"/"+newNameObj)));
      }
      //h1D->SetDirectory(outFile);
      outFile->cd();
      h1D->Write();
      delete h1D;
    }
  }

  cout << "Writing File" << endl;
  //outFile->Write();

  cout << "File Written" << endl;
  /*
  TList* keyList = mcFile->GetListOfKeys();
  if (!keyList){
    cout << "List of keys failed to get." << endl;
    return 5;
  }

  TIter next(keyList);
  TKey* key;
  }
  */
  //cout << "Deleting the MnvPlotter." << endl;
  //delete plotter;

  //delete mcPOT;
  //delete
  cout << "Changing Dir. " << endl;
  outFile->cd();
  cout << "Getting POT" << endl;
  TParameter<double>* POT = (TParameter<double>*)(files[firstName]->Get("POTUsed"))->Clone("POTUsed");
  cout << "Writing POT" << endl;
  POT->Write();
  cout << "Closing outFile" << endl;
  outFile->Close();

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  for (auto file : files){
    file.second->Close();
  } 

  //mcFile->Close();
  //dataFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}