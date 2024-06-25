//File: allErrSummaries.C
//Info: Quick macro to make initial pmu sideband/signal region plots.
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

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"

using namespace std;
using namespace PlotUtils;

void allErrSummaries(TString fileName, TString outDir) {

  gROOT->SetBatch();

  MnvPlotter* plotter = new MnvPlotter(kCCQEAntiNuStyle);
  TFile* inFile = new TFile(fileName,"READ");

  TList* keyList = inFile->GetListOfKeys();
  if (!keyList){
    cout << "List of keys failed to get." << endl;
    return;
  }

  TIter next(keyList);
  TKey* key;
  while ( key = (TKey*)next() ){
    //cout << key->GetName() << endl;
    if ((TString)(key->GetClassName()) == "TDirectoryFile"){
      TDirectoryFile* dirInt = (TDirectoryFile*)inFile->Get(key->GetName());
      TList* keyListInt = dirInt->GetListOfKeys();
      if (!keyListInt){
	cout << "List of keys failed to get inside directory." << endl;
	return;
      }
      
      TIter nextInt(keyListInt);
      TKey* keyInt;
      size_t pos=0;
      string name=(string)key->GetName();
      if ((pos = name.find("TwoD")) != string::npos ) continue;

      while ( keyInt = (TKey*)nextInt() ){
	string nameInt = (string)keyInt->GetName();
	name = (string)key->GetName() + "/"+nameInt;
	if((pos = nameInt.find("_data")) != string::npos){
	  pos=0;
	  if ((pos = nameInt.find("SB")) != string::npos || (pos = nameInt.find("vtxZ")) != string::npos ) continue;
	  cout << "Plotting error summary for: " << nameInt << endl;      
	  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
	  MnvH1D* h_mc_data = (MnvH1D*)inFile->Get((TString)name);      
	  name.erase(name.length()-5,name.length());      
	  nameInt.erase(nameInt.length()-5,nameInt.length());
	  plotter->DrawErrorSummary(h_mc_data);
	  //plotter->DrawErrorSummary(h_mc_data, "TR", true, true, 1e-5, false, "Neutron Detection");
	  c1->Print((TString)outDir+(TString)nameInt+"_err_summary.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_err_summary.png");
	  cout << "" << endl;
	  delete c1;
	}
	else if ((pos = nameInt.find("_efficiency_numerator")) != string::npos){
	  pos=0;
	  if ((pos = nameInt.find("SB")) != string::npos || (pos = nameInt.find("vtxZ")) != string::npos ) continue;
	  cout << "Plotting error summary for: " << nameInt << endl;
	  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
	  MnvH1D* h_mc_eff_num = (MnvH1D*)inFile->Get((TString)name);
	  plotter->DrawErrorSummary(h_mc_eff_num);
	  //plotter->DrawErrorSummary(h_mc_eff_num, "TR", true, true, 1e-5, false, "Neutron Detection");
	  c1->Print((TString)outDir+(TString)nameInt+"_err_summary.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_err_summary.png");
	  name.erase(name.length()-10,name.length());
	  string nameDen = name+"_denominator";
	  nameInt.erase(nameInt.length()-10,nameInt.length());
	  string nameIntDen = nameInt+"_denominator";
	  MnvH1D* h_mc_eff_den = (MnvH1D*)inFile->Get((TString)nameDen);
	  h_mc_eff_den->AddMissingErrorBandsAndFillWithCV(*h_mc_eff_num);
	  plotter->DrawErrorSummary(h_mc_eff_den);
	  //plotter->DrawErrorSummary(h_mc_eff_den, "TR", true, true, 1e-5, false, "Neutron Detection");
	  c1->Print((TString)outDir+(TString)nameIntDen+"_err_summary.pdf");
	  c1->Print((TString)outDir+(TString)nameIntDen+"_err_summary.png");
	  h_mc_eff_num->Divide(h_mc_eff_num,h_mc_eff_den);
	  plotter->DrawErrorSummary(h_mc_eff_num);
	  //plotter->DrawErrorSummary(h_mc_eff_num, "TR", true, true, 1e-5, false, "Neutron Detection");
	  c1->Print((TString)outDir+(TString)nameInt+"_err_summary.pdf");
	  c1->Print((TString)outDir+(TString)nameInt+"_err_summary.png");      
	  cout << "" << endl;
	}	
      }
    }

    size_t pos=0;
    string name=(string)key->GetName();
    if((pos = name.find("_data")) != string::npos){
      pos=0;
      if ((pos = name.find("SB")) != string::npos || (pos = name.find("vtxZ")) != string::npos ) continue;
      cout << "Plotting error summary for: " << name << endl;      
      TCanvas* c1 = new TCanvas("c1","c1",1200,800);
      MnvH1D* h_mc_data = (MnvH1D*)inFile->Get((TString)name);      
      name.erase(name.length()-5,name.length());      
      plotter->DrawErrorSummary(h_mc_data);
      //plotter->DrawErrorSummary(h_mc_data, "TR", true, true, 1e-5, false, "Neutron Detection");
      c1->Print((TString)outDir+(TString)name+"_err_summary.pdf");
      c1->Print((TString)outDir+(TString)name+"_err_summary.png");
      cout << "" << endl;
      delete c1;
    }
    else if ((pos = name.find("_efficiency_numerator")) != string::npos){
      pos=0;
      if ((pos = name.find("SB")) != string::npos || (pos = name.find("vtxZ")) != string::npos ) continue;
      cout << "Plotting error summary for: " << name << endl;
      TCanvas* c1 = new TCanvas("c1","c1",1200,800);
      MnvH1D* h_mc_eff_num = (MnvH1D*)inFile->Get((TString)name);
      plotter->DrawErrorSummary(h_mc_eff_num);
      //plotter->DrawErrorSummary(h_mc_eff_num, "TR", true, true, 1e-5, false, "Neutron Detection");
      c1->Print((TString)outDir+(TString)name+"_err_summary.pdf");
      c1->Print((TString)outDir+(TString)name+"_err_summary.png");
      name.erase(name.length()-10,name.length());
      string nameDen = name+"_denominator";
      MnvH1D* h_mc_eff_den = (MnvH1D*)inFile->Get((TString)nameDen);
      h_mc_eff_den->AddMissingErrorBandsAndFillWithCV(*h_mc_eff_num);
      plotter->DrawErrorSummary(h_mc_eff_den);
      //plotter->DrawErrorSummary(h_mc_eff_den, "TR", true, true, 1e-5, false, "Neutron Detection");
      c1->Print((TString)outDir+(TString)nameDen+"_err_summary.pdf");
      c1->Print((TString)outDir+(TString)nameDen+"_err_summary.png");
      h_mc_eff_num->Divide(h_mc_eff_num,h_mc_eff_den);
      plotter->DrawErrorSummary(h_mc_eff_num);
      //plotter->DrawErrorSummary(h_mc_eff_num, "TR", true, true, 1e-5, false, "Neutron Detection");
      c1->Print((TString)outDir+(TString)name+"_err_summary.pdf");
      c1->Print((TString)outDir+(TString)name+"_err_summary.png");      
      cout << "" << endl;
    }
  }

  cout << "deleting plotter" << endl;
  delete plotter;

  cout << "HEY YOU DID IT!!!" << endl;
  return;
}
