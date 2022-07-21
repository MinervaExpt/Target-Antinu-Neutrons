//File: QuickPlotMacro.C
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

using namespace std;
using namespace PlotUtils;

map<int, vector<TString>> mats = {{1,{"_Fe","_Pb","_Buffer"}},
				  {2,{"_Fe","_Pb","_Buffer"}},
				  {3,{"_C","_Fe","_Pb","_Buffer"}},
				  {4,{"_Pb"}},
				  {5,{"_Fe","_Pb","_Buffer"}},
				  {6,{""}}};

vector<TString> TgtTypes = {"C","Fe","Pb","Water","Plastic","USPlastic","DSPlastic","Other"};

void ScanPlasticBKGs(TString fileName, int Tgt) {

  TFile* inFile = new TFile(fileName,"READ");

  TString TgtName = "Tgt"+to_string(Tgt);
  if (Tgt==6) TgtName = "WaterTgt";

  MnvH1D* h;

  for (auto matName: mats[Tgt]){
    for (auto typeName: TgtTypes){
      TString USplotName = "ByTgt_"+TgtName+matName+"/US_ByType_"+typeName+"/vtxZ_ByTgt_"+TgtName+matName+"_PreRecoilCut_"+typeName+"_data";
      TString DSplotName = "ByTgt_"+TgtName+matName+"/DS_ByType_"+typeName+"/vtxZ_ByTgt_"+TgtName+matName+"_PreRecoilCut_"+typeName+"_data";
      h = (MnvH1D*)inFile->Get(USplotName);
      cout << "US " << TgtName+matName << ", " << typeName << ": " << h->Integral(0,-1) << endl;
      h = (MnvH1D*)inFile->Get(DSplotName);
      cout << "DS " << TgtName+matName << ", " << typeName << ": " << h->Integral(0,-1) << endl;
    }
  }

  return;
}
