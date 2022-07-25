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
      c1->Print((TString)outDir+(TString)name+"_err_summary.pdf");
      c1->Print((TString)outDir+(TString)name+"_err_summary.png");
      cout << "" << endl;
      delete c1;
    }
  }

  cout << "deleting plotter" << endl;
  delete plotter;

  cout << "HEY YOU DID IT!!!" << endl;
  return;
}
