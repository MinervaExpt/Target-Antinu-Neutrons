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

MnvH1D* AdjustData(MnvH1D* hData, MnvH1D* hSig_PreScale, MnvH1D* hSig_PostScale, std::string name){
  MnvH1D* hOut = new MnvH1D(*hData);
  bool dataHasBand = hOut->HasErrorBand(name);
  bool preHasBand = hSig_PreScale->HasErrorBand(name);
  bool postHasBand =  hSig_PostScale->HasErrorBand(name);
  bool hasBand = dataHasBand && postHasBand;
  if (!hasBand){
    cout << "This should have this error band at this point..." << endl;
    cout << dataHasBand << ", " << preHasBand << ", " << postHasBand << endl;
    return hOut;
  }
  const auto univsData = hOut->GetVertErrorBand(name)->GetHists();
  const auto univsPostScale = hSig_PostScale->GetVertErrorBand(name)->GetHists();
  const auto univsPreScale = (preHasBand) ? hSig_PreScale->GetVertErrorBand(name)->GetHists() : hSig_PostScale->GetVertErrorBand(name)->GetHists();
  for (unsigned int i=0; i<univsData.size(); ++i){
    TH1D* hTmp = (preHasBand) ? new TH1D(*univsPreScale.at(i)) : new TH1D(hSig_PreScale->GetCVHistoWithStatError());
    univsData.at(i)->Add(hTmp,-1.0);
    univsData.at(i)->Add(univsPostScale.at(i),1.0);
    delete hTmp;
  }

  return hOut;
}

MnvH2D* AdjustData(MnvH2D* hData, MnvH2D* hSig_PreScale, MnvH2D* hSig_PostScale, std::string name){
  MnvH2D* hOut = new MnvH2D(*hData);
  bool dataHasBand = hOut->HasErrorBand(name);
  bool postHasBand =  hSig_PostScale->HasErrorBand(name);
  bool preHasBand = hSig_PreScale->HasErrorBand(name);
  bool hasBand = dataHasBand && postHasBand;
  if (!hasBand){
    cout << "These should all have this error band at this point..." << endl;
    return hOut;
  }
  const auto univsData = hOut->GetVertErrorBand(name)->GetHists();
  const auto univsPostScale = hSig_PostScale->GetVertErrorBand(name)->GetHists();
  const auto univsPreScale = (preHasBand) ? hSig_PreScale->GetVertErrorBand(name)->GetHists() : hSig_PostScale->GetVertErrorBand(name)->GetHists();
  for (unsigned int i=0; i<univsData.size(); ++i){
    TH2D* hTmp = (preHasBand) ? new TH2D(*univsPreScale.at(i)) : new TH2D (hSig_PreScale->GetCVHistoWithStatError());
    univsData.at(i)->Add(hTmp,-1.0);
    univsData.at(i)->Add(univsPostScale.at(i),1.0);
    delete hTmp;
  }

  return hOut;
}

MnvH1D* ScaleErrorBand(MnvH1D* h1D, std::string name, double scale, bool doScale=true){
  MnvH1D* hOut = new MnvH1D(*h1D);
  bool hasBand = hOut->HasErrorBand(name);

  if (!hasBand){
    if (!doScale) scale = 0.0;
    std::vector<TH1D*> univs;
    TH1D* hUp = new TH1D(hOut->GetCVHistoWithStatError());
    hUp->SetDirectory(nullptr);
    TH1D* hDown = new TH1D(*hUp);
    hUp->Scale(1+scale);
    univs.push_back(hUp);
    hDown->Scale(1-scale);
    univs.push_back(hDown);
    hOut->AddVertErrorBand(name, univs);
    delete hUp;
    delete hDown;
    univs.clear();
  }
  else{
    if (!doScale) scale = 1.0;
    TH1D* hCV = new TH1D(hOut->GetCVHistoWithStatError());
    const auto univs = hOut->GetVertErrorBand(name)->GetHists();
    for (auto hist:univs){
      hist->Add(hCV,-1.0);
      hist->Scale(scale);
      hist->Add(hCV,1.0);
    }
    delete hCV;
  }
  return hOut;
}

MnvH2D* ScaleErrorBand(MnvH2D* h2D, std::string name, double scale, bool doScale=true){
  MnvH2D* hOut = new MnvH2D(*h2D);
  bool hasBand = hOut->HasErrorBand(name);

  if (!hasBand){
    if (!doScale) scale=0.0;
    std::vector<TH2D*> univs;
    TH2D* hUp = new TH2D(hOut->GetCVHistoWithStatError());
    hUp->SetDirectory(nullptr);
    TH2D* hDown = new TH2D (*hUp);
    hUp->Scale(1+scale);
    univs.push_back(hUp);
    hDown->Scale(1-scale);
    univs.push_back(hDown);
    hOut->AddVertErrorBand(name, univs);
    delete hUp;
    delete hDown;
    univs.clear();
  }
  else{
    if (!doScale) scale=1.0;
    TH2D* hCV = new TH2D(hOut->GetCVHistoWithStatError());
    const auto univs = hOut->GetVertErrorBand(name)->GetHists();
    for (auto hist:univs){
      hist->Add(hCV,-1.0);
      hist->Scale(scale);
      hist->Add(hCV,1.0);
    }
    delete hCV;
  }
  return hOut;
}

vector<TString> selSigNames = {"_sig_","_selected_signal_reco","_efficiency_numerator","_migration"};

bool IsSelSig(TString nameObj){
  bool isSelSig = false;
  for(auto name:selSigNames){
    if(nameObj.Contains(name)){
      isSelSig = true;
      break;
    }
  }
  return isSelSig;
}

void RecursePerform(TList* keyList, TFile* inFile, TDirectory* outFile, TString nameDir, std::string bandName, double scale, bool selSigOnly){
  TIter next(keyList);
  TKey* key;
  while ( key = (TKey*)next() ){
    TString className = (TString)key->GetClassName();
    TString nameObj = (TString)key->GetName();
    TString nameDirFixed = (nameDir != "") ? nameDir+"/" : "";
    if (nameObj.Contains("ByTgt") && (nameObj.Contains("_Buffer") || nameObj.Contains("_Other"))) continue;//Removing these to save time. Shouldn't be needed by next steps.
    if (className == "TDirectoryFile"){
      cout << "Directory File!" << endl;
      TString nameDirInt = nameDirFixed+nameObj;
      TDirectory* newOutDir = outFile->mkdir(nameDirInt);
      TDirectoryFile* dirInt = (TDirectoryFile*)inFile->Get(nameDirInt);
      TList*  keyListInt = dirInt->GetListOfKeys();
      cout << "Entering New Recursion" << endl;
      RecursePerform(keyListInt, inFile, outFile, nameDirInt, bandName, scale, selSigOnly);
    }
    else if (!(className.Contains("MnvH") || className == "TParameter<double>") || nameObj.Contains("MYBins")) continue;
    else if (className == "TParameter<double>"){
      TParameter<double>* tPar = (TParameter<double>*)(inFile->Get(nameDirFixed+nameObj)->Clone(nameObj));
      outFile->cd(nameDirFixed);
      tPar->Write();
      delete tPar;
    }
    else if (className.Contains("MnvH2")){
      if (selSigOnly && nameObj.Contains("_data")) continue;
      bool scaleThisObj = (selSigOnly) ? IsSelSig(nameObj) : true;
      MnvH2D* h2D = (MnvH2D*)(inFile->Get(nameDirFixed+nameObj)->Clone(nameObj));
      MnvH2D* h2D2 = ScaleErrorBand(h2D,bandName,scale,scaleThisObj);
      outFile->cd(nameDirFixed);
      h2D2->Write();
      if (selSigOnly && nameObj.Contains("selected_signal_reco")){
	TString nameBase = BreakName("_selected_signal_reco", nameObj)[0];
	cout << nameBase << endl;
	MnvH2D* h2D_Data = (MnvH2D*)(inFile->Get(nameDirFixed+nameBase+"_data")->Clone(nameBase+"_data"));
	MnvH2D* h2D2_Data = ScaleErrorBand(h2D_Data,bandName,scale,false);
	MnvH2D* h2D3_Data = AdjustData(h2D2_Data,h2D,h2D2,bandName);
	cout << "Writing: " << nameDirFixed+nameBase+"_data" << endl;
	outFile->cd(nameDirFixed);
	h2D3_Data->Write();
	delete h2D_Data;
	delete h2D2_Data;
	delete h2D3_Data;
      }
      delete h2D;
      delete h2D2;
    }
    else if (className.Contains("MnvH1")){
      if (selSigOnly && nameObj.Contains("_data")) continue;
      bool scaleThisObj = (selSigOnly) ? IsSelSig(nameObj) : true;
      MnvH1D* h1D = (MnvH1D*)(inFile->Get(nameDirFixed+nameObj)->Clone(nameObj));
      MnvH1D* h1D2 = ScaleErrorBand(h1D,bandName,scale,scaleThisObj);
      outFile->cd(nameDirFixed);
      cout << "Writing: " << nameDirFixed+nameObj << endl;
      h1D2->Write();
      if (selSigOnly && nameObj.Contains("selected_signal_reco")){
	TString nameBase = BreakName("_selected_signal_reco", nameObj)[0];
	cout << nameBase << endl;
	MnvH1D* h1D_Data = (MnvH1D*)(inFile->Get(nameDirFixed+nameBase+"_data")->Clone(nameBase+"_data"));
	MnvH1D* h1D2_Data = ScaleErrorBand(h1D_Data,bandName,scale,false);
	MnvH1D* h1D3_Data = AdjustData(h1D2_Data,h1D,h1D2,bandName);
	cout << "Writing: " << nameDirFixed+nameBase+"_data" << endl;
	outFile->cd(nameDirFixed);
	h1D3_Data->Write();
	delete h1D_Data;
	delete h1D2_Data;
	delete h1D3_Data;
      }
      delete h1D;
      delete h1D2;
    }
    else{
      cout << "HUH?" << endl;
    }
  }
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
  if (argc != 6) {
    cout << "Check usage..." << endl;
    return 2;
  }

  TString inFileName = argv[1];
  std::string bandName = argv[2];
  double scale = atof(argv[3]);
  TString outFileName = argv[4];
  bool selSigOnly = (bool)(atoi(argv[5]));

  TFile* inFile = new TFile(inFileName,"READ");
  TFile* outFile = new TFile(outFileName,"RECREATE");

  TList* keyList = inFile->GetListOfKeys();
  if (!keyList){
    cout << "issue with input file" << endl;
    return 4;
  }

  TString nameDir = "";
  RecursePerform(keyList, inFile, outFile, nameDir, bandName, scale, selSigOnly);

  outFile->cd();

  cout << "Closing outFile" << endl;
  outFile->Close();

  cout << "Closing input file" << endl;
  inFile->Close();

  cout << "HEY YOU DID IT" << endl;
  return 0;
}
