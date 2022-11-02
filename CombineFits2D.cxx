//File: CombineFits2D.cxx
//Info: This is a script to combine fit results across 2D pT/pz bins into a MnvH2D* scaling histo to apply scale factors/errors by mulitplication.
//
//Usage: CombineFits2D <inFileName> <refFileName> <refHistoFilePath>
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

vector<int> binsIn1DRange(TAxis* ax, double xMin, double xMax){
  vector<int> outVec;
  for (int iBin=0; iBin <= ax->GetNbins()+1; ++iBin){
    if (ax->GetBinLowEdge(iBin) >= xMin && ax->GetBinUpEdge(iBin) <= xMax) outVec.push_back(iBin);
  }
  return outVec;
}

vector<int> binsIn2DRanges(MnvH2D* hRef, double xMin, double xMax, double yMin, double yMax){
  vector<int> outVec;

  TAxis* xAxis = hRef->GetXaxis();
  vector<int> xBins = binsIn1DRange(xAxis,xMin,xMax);

  TAxis* yAxis = hRef->GetYaxis();
  vector<int> yBins = binsIn1DRange(yAxis,yMin,yMax);
  for (auto yBin:yBins){
    for (auto xBin:xBins){
      //cout << "xBin: " << xBin;
      //cout << ", yBin: " << yBin;
      outVec.push_back(hRef->GetBin(xBin,yBin));
      //cout << ", bin: " << outVec.back() << endl;
    }
  }
  cout << "" << endl;
  return outVec;
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

  TString inFileName= argv[1];
  TString refFileName = argv[2];
  TString refName = argv[3];
  TString refString = "_bin_0_";

  if (!inFileName.Contains(".root")){
    cout << "Input file is not a root file." << endl;
    return 3;
  }

  if (!refFileName.Contains(".root")){
    cout << "Reference file is not a root file." << endl;
    return 4;
  }

  TFile* inFile = new TFile(inFileName,"READ");
  TFile* refFile = new TFile(refFileName,"READ");

  MnvH2D* refHisto = (MnvH2D*)refFile->Get(refName);

  map<int,vector<int>> mapOfFitBinsToFillBins;

  vector<int> binsLost;
  for (int yBin=0; yBin<=refHisto->GetNbinsY()+1; ++yBin){
    binsLost.push_back(refHisto->GetBin(0,yBin));
  }

  for (int yBin=0; yBin<=refHisto->GetNbinsY()+1; ++yBin){
    binsLost.push_back(refHisto->GetBin(refHisto->GetNbinsX()+1,yBin));
  }

  for (int xBin=1; xBin<=refHisto->GetNbinsX(); ++xBin){
    binsLost.push_back(refHisto->GetBin(xBin,0));
  }

  for (int xBin=1; xBin<=refHisto->GetNbinsX(); ++xBin){
    binsLost.push_back(refHisto->GetBin(xBin,refHisto->GetNbinsY()+1));
  }

  mapOfFitBinsToFillBins[0] = binsIn2DRanges(refHisto, 1.5,5.0, 0.0,0.25);//bin 0
  mapOfFitBinsToFillBins[1] = binsIn2DRanges(refHisto, 1.5,5.0, 0.25,0.4);//bin 1
  mapOfFitBinsToFillBins[2] = binsIn2DRanges(refHisto, 1.5,5.0, 0.4,0.7);//bin 2
  mapOfFitBinsToFillBins[3] = binsIn2DRanges(refHisto, 1.5,5.0, 0.7,0.85);//bin 3
  mapOfFitBinsToFillBins[4] = binsIn2DRanges(refHisto, 1.5,5.0, 0.85,1.0);//bin 4
  mapOfFitBinsToFillBins[5] = binsIn2DRanges(refHisto, 1.5,5.0, 1.0,2.5);//bin 5
  mapOfFitBinsToFillBins[6] = binsIn2DRanges(refHisto, 5.0,8.0, 0.0,0.25);//bin 6
  mapOfFitBinsToFillBins[7] = binsIn2DRanges(refHisto, 5.0,8.0, 0.25,0.4);//bin 7
  mapOfFitBinsToFillBins[8] = binsIn2DRanges(refHisto, 5.0,8.0, 0.4,0.7);//bin 8
  mapOfFitBinsToFillBins[9] = binsIn2DRanges(refHisto, 5.0,8.0, 0.7,0.85);//bin 9
  mapOfFitBinsToFillBins[10] = binsIn2DRanges(refHisto, 5.0,8.0, 0.85,1.0);//bin 10
  mapOfFitBinsToFillBins[11] = binsIn2DRanges(refHisto, 5.0,8.0, 1.0,2.5);//bin 11
  mapOfFitBinsToFillBins[12] = binsIn2DRanges(refHisto, 8.0,15.0, 0.0,0.55);//bin 12
  mapOfFitBinsToFillBins[13] = binsIn2DRanges(refHisto, 8.0,15.0, 0.55,2.5);//bin 13

  int nBins = mapOfFitBinsToFillBins.size();

  vector<TString> inFilePieces = BreakName(".root",inFileName);

  TString outFileName = inFilePieces.at(0)+"2D.root";
  TFile* outFile = new TFile(outFileName,"RECREATE");

  TList* keyList = inFile->GetListOfKeys();
  if(!keyList){
    cout << "List of keys failed to get inside inFile." << endl;
    return 5;
  }

  TIter nextKey(keyList);
  TKey* key;
  while ( key = (TKey*)nextKey() ){
    TString nameObj = (TString)key->GetName();
    if (nameObj.Contains(refString)){
      cout << "Creating new 2D histo starting with: " << nameObj << endl;
      vector<TString> nameObjPieces = BreakName("_fit_",nameObj);
      TString fitName = nameObjPieces.back();
      vector<TString> fitNamePieces = BreakName(refString,fitName);
      fitName = MashNames("_",fitNamePieces);
      cout << "Fit Name Parsed To: " << fitName << endl;
      
      map<int,MnvH1D*> fitHistos;
      map<int,double> fitCVResult;
      map<int,double> fitCVError;

      cout << "Getting: " << nameObj << endl;
      fitHistos[0] = (MnvH1D*)inFile->Get(nameObj);
      fitCVResult[0] = fitHistos[0]->GetBinContent(1);
      fitCVError[0] = fitHistos[0]->GetCVHistoWithStatError().GetBinError(1);

      nameObjPieces = BreakName(refString,nameObj);
      for (int iFit=1; iFit < nBins;++iFit){
	TString binNameObj = MashNames("_bin_"+to_string(iFit)+"_",nameObjPieces);
	cout << "Getting: " << binNameObj << endl;
	fitHistos[iFit] = (MnvH1D*)inFile->Get(binNameObj);
	fitCVResult[iFit] = fitHistos[iFit]->GetBinContent(1);
	fitCVError[iFit] = fitHistos[iFit]->GetCVHistoWithStatError().GetBinError(1);
      }

      TString refOutName = BreakName("data",refHisto->GetName()).at(0)+fitName;
      cout << refOutName << endl;
      MnvH2D* hOut = new MnvH2D(refOutName,"",refHisto->GetNbinsX(),refHisto->GetXaxis()->GetXbins()->GetArray(),refHisto->GetNbinsY(),refHisto->GetYaxis()->GetXbins()->GetArray());

      for (int iFit=0; iFit<fitHistos.size(); ++iFit){
	for (auto bin: mapOfFitBinsToFillBins[iFit]){
	  hOut->SetBinContent(bin,fitCVResult[iFit]);
	  hOut->SetBinError(bin,fitCVError[iFit]);
	}
      }

      for (auto bin: binsLost){
	hOut->SetBinContent(bin,1.0);
	hOut->SetBinError(bin,0.0);
      }

      hOut->AddMissingErrorBandsAndFillWithCV(*refHisto);

      const auto errorBandNames = hOut->GetErrorBandNames();
      for (const auto&  bandName: errorBandNames){
	const auto univs = hOut->GetVertErrorBand(bandName)->GetHists();
	for (size_t whichUniv=0; whichUniv < univs.size(); ++ whichUniv){
	  for (int iFit=0; iFit<fitHistos.size(); ++iFit){
	    for (auto bin: mapOfFitBinsToFillBins[iFit]){
	      hOut->GetVertErrorBand(bandName)->GetHist(whichUniv)->SetBinContent(bin,fitHistos[iFit]->GetVertErrorBand(bandName)->GetHist(whichUniv)->GetBinContent(1));
	    }
	  }
	}
      }

      outFile->cd();
      hOut->Write();
      delete hOut;
      cout << "" << endl;
    }
  }


  inFile->Close();
  refFile->Close();
  outFile->Close();
  return 0;
}

  /*
  string rootExt = ".root";
  string slash = "/";
  string token;
  string fileNameStub = inFileName;
  size_t pos=0;

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

  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory doesn't exist. Exiting" << endl;
    return 3;
  }

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
	if (!classNameInt.Contains("PlotUtils::MnvH")) continue;
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
	
	else if (nameObj.Contains("DS_ByType")){
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

	else if (classNameInt.Contains("2")){
	  MnvH2D* h2D = (MnvH2D*)(files[firstName]->Get(dirBase+firstName+matName+"/"+nameObj+"/"+nameObjInt))->Clone(matNameObjInt);
	  for(auto file : files){
	    if (file.first == firstName) continue;
	    TString newNameObjInt = MashNames(file.first+matName,nameObjIntPieces);
	    h2D->Add((MnvH2D*)(file.second->Get(dirBase+file.first+matName+"/"+nameObj+"/"+newNameObjInt)));
	  }

	  newOutDir->cd();
	  h2D->Write();
	  delete h2D;
	}

	else{
	  MnvH1D* h1D = (MnvH1D*)(files[firstName]->Get(dirBase+firstName+matName+"/"+nameObj+"/"+nameObjInt))->Clone(matNameObjInt);
	  for(auto file : files){
	    if (file.first == firstName) continue;
	    TString newNameObjInt = MashNames(file.first+matName,nameObjIntPieces);
	    h1D->Add((MnvH1D*)(file.second->Get(dirBase+file.first+matName+"/"+nameObj+"/"+newNameObjInt)));
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

  TList* keyList = mcFile->GetListOfKeys();
  if (!keyList){
    cout << "List of keys failed to get." << endl;
    return 5;
  }

  TIter next(keyList);
  TKey* key;
  }

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
  */
