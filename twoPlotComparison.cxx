//File: twoPlotComparison.cxx
//Info: This script lives to compare the same plots from two different files with stat error only.  
//
//Usage: twoPlotComparison.cxx <file1> <file2> <output_directory> <plot_label>
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

//TODO: FIX SEGFAULT ISSUE AT END OF EXECUTION... UNCLEAR WHY THAT'S HAPPENING AND IT DOESN'T SEEM TO AFFECT ANYTHING... MAYBE NEED TO CLOSE FILES? Cloning maybe tries to add keys to the file and it doesn't close well when that's the case?

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

TCanvas* DrawRatio(string name, TFile* file1, TFile* file2){

  MnvH1D* m1 = (MnvH1D*)file1->Get((TString)name);
  if (!m1){
    cout << "Not in file 1." << endl;
    return nullptr;
  }
  MnvH1D* m2 = (MnvH1D*)file2->Get((TString)name);
  if (!m2){
    cout << "Not in file 2." << endl;
    return nullptr;
  }

  TH1D* h1 = (TH1D*)m1->GetCVHistoWithStatError().Clone();
  TH1D* h2 = (TH1D*)m2->GetCVHistoWithStatError().Clone();

  cout << "Maximum: Pre-Divide" << h1->GetMaximum() << endl;

  h1->Divide(h2);

  cout << "Maximum: Post-Divide" << h1->GetMaximum() << endl;

  h1->SetTitle("Ratio for "+(TString)name+" of file 1 to file 2.");
  h1->GetYaxis()->SetRangeUser(0.5,1.5);

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);

  c1->cd();
  h1->Draw();
  c1->Update();
  
  return c1;
}

bool PathExists(string path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

int main(int argc, char* argv[]) {

  gStyle->SetOptStat(0);

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc != 5) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string file1Name=string(argv[1]);
  string file2Name=string(argv[2]);
  string outDir=string(argv[3]);
  TString label=argv[4];

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
  string fileNameStub = file1Name;
  size_t pos=0;

  //cout << sigNameStub << endl;
  while ((pos = fileNameStub.find(slash)) != string::npos){
    //cout << sigNameStub << endl;
    token = fileNameStub.substr(0, pos);
    //cout << token << endl;
    fileNameStub.erase(0, pos+slash.length());
  }
  //cout << sigNameStub << endl;
  if ((pos=fileNameStub.find(rootExt)) == string::npos){
    cout << "Input files need be .root files." << endl;
    return 4;
  }

  cout << "Input file name 1 parsed to: " << fileNameStub << endl;

  rootExt = ".root";
  slash = "/";
  token = "";
  fileNameStub = file2Name;
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
    cout << "Input files need be .root files." << endl;
    return 5;
  }

  cout << "Input file name 2 parsed to: " << fileNameStub << endl;

  //cout << "Setting up MnvPlotter" << endl;
  //MnvPlotter* plotter = new MnvPlotter(kCCQEAntiNuStyle);

  TFile* file1 = new TFile(file1Name.c_str(),"READ");
  TFile* file2 = new TFile(file2Name.c_str(),"READ");

  TList* keyList = file1->GetListOfKeys();
  if (!keyList){
    cout << "List of keys failed to get." << endl;
    return 5;
  }

  TIter next(keyList);
  TKey* key;
  while ( key = (TKey*)next() ){
    TString nameClass = (TString)(key->GetClassName());
    if (nameClass == "TDirectoryFile")
    {
      TDirectoryFile* dirInt = (TDirectoryFile*)file1->Get(key->GetName());
      TList* keyListInt = dirInt->GetListOfKeys();
      if (!keyListInt){
	cout << "List of keys failed to get inside directory." << endl;
	return 5;
      }
      
      TIter nextInt(keyListInt);
      TKey* keyInt;
      while ( keyInt = (TKey*)nextInt() ){
	TString nameClassInt = (TString)(keyInt->GetClassName());
	if (nameClassInt.Contains("2") || nameClassInt.Contains("TParameter")) continue; 
	string nameInt = (string)keyInt->GetName();
	string name = (string)key->GetName() + "/"+nameInt;
	cout << name << " of class type: " << nameClassInt <<endl;
	/*
	pos=0;
	if ((pos=nameInt.find("TwoD")) != string::npos) continue;
	else if((pos = name.find("_sig_IntType_QE")) != string::npos){
	  TCanvas* c1 = DrawIntType(name,mcFile,dataFile,label,scale);
	  TPad* top = (TPad*)c1->GetPrimitive("Overlay");
	  name.erase(name.length()-15,name.length());
	  c1->Print((TString)outDir+(TString)nameInt+"_IntType_stacked.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_IntType_stacked.png");
	  top->SetLogy();
	  c1->Update();
	  c1->Print((TString)outDir+(TString)nameInt+"_IntType_stacked_log.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_IntType_stacked_log.png");
	  cout << "" << endl;
	  delete c1;
	}
	else if ((pos = name.find("_sig_TargetType_C")) != string::npos){
	  TCanvas* c1 = DrawTargetType(name,mcFile,dataFile,label,scale);
	  TPad* top = (TPad*)c1->GetPrimitive("Overlay");
	  nameInt.erase(nameInt.length()-17,nameInt.length());
	  c1->Print((TString)outDir+(TString)nameInt+"_TargetType_stacked.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_TargetType_stacked.png");
	  top->SetLogy();
	  c1->Update();
	  c1->Print((TString)outDir+(TString)nameInt+"_TargetType_stacked_log.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_TargetType_stacked_log.png");
	  cout << "" << endl;
	  delete c1;
	}
	else if ((pos = name.find("_sig_LeadBlobType_neut")) != string::npos){
	  TCanvas* c1 = DrawLeadBlobType(name,mcFile,dataFile,label,scale);
	  TPad* top = (TPad*)c1->GetPrimitive("Overlay");
	  nameInt.erase(nameInt.length()-22,nameInt.length());
	  c1->Print((TString)outDir+(TString)nameInt+"_LeadBlobType_stacked.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_LeadBlobType_stacked.png");
	  top->SetLogy();
	  c1->Update();
	  c1->Print((TString)outDir+(TString)nameInt+"_LeadBlobType_stacked_log.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_LeadBlobType_stacked_log.png");
	  cout << "" << endl;
	  delete c1;
	}
	else if ((pos = name.find("_selected_signal_reco")) != string::npos){
	  TCanvas* c1 = DrawBKGCateg(name,mcFile,dataFile,label,scale);
	  TPad* top = (TPad*)c1->GetPrimitive("Overlay");
	  nameInt.erase(nameInt.length()-21,nameInt.length());
	  c1->Print((TString)outDir+(TString)nameInt+"_BKG_stacked.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_BKG_stacked.png");
	  top->SetLogy();
	  c1->Update();
	  c1->Print((TString)outDir+(TString)nameInt+"_BKG_stacked_log.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_BKG_stacked_log.png");
	  cout << "" << endl;
	  delete c1;
	}
	*/	
      }
    }
    else if (nameClass.Contains("2") || nameClass.Contains("TParameter")) continue;
    else {
      string name=(string)key->GetName();
      //cout << name << " of class type: " << nameClass << endl;
      //if ((pos=name.find("SB")) != string::npos) continue;
      //cout << "Plotting error summary for: " << name << endl;
      TCanvas* c1 = DrawRatio(name,file1,file2);
      if (!c1){
	cout << "See above for why ratio couldn't be plotted." << endl;
	continue;
      }
      c1->Print((TString)outDir+(TString)name+"_ratio_for_"+label+".pdf");
      c1->Print((TString)outDir+(TString)name+"_ratio_for_"+label+".png");
      cout << "" << endl;
      delete c1;
    }
  }

  //cout << "Deleting the MnvPlotter." << endl;
  //delete plotter;

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  file1->Close();
  file2->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
