//File: AddUpPOT.cxx
//Info: This pairs with a call to madd to replace the POT calculation that madd doesn't handle.
//
//Usage: AddUpPOT <outFile> <added files>
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

int main(int argc, char* argv[]) {

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  if (argc < 3) {
    cout << "Check usage... there needs to be an extra file..." << endl;
    return 2;
  }

  TString outFileName = argv[1];
  vector<TString> fileNames;
  for (int i=2; i<argc; ++i) fileNames.push_back(argv[i]);

  double POTsum = 0.0;

  for (auto name: fileNames){
    TFile* file = new TFile(name,"READ");
    TParameter<double>* POT = (TParameter<double>*)file->Get("POTUsed");
    cout << "POT for file " << name << ": " << POT->GetVal() << endl;
    POTsum += POT->GetVal();
    file->Close();
  }

  TFile* outFile = new TFile(outFileName,"UPDATE");
  outFile->cd();
  auto outPOT = new TParameter<double>("POTUsed", POTsum);
  outPOT->Write();
  cout << "Closing outFile" << endl;
  outFile->Close();

  cout << "HEY YOU DID IT" << endl;
  return 0;
}
