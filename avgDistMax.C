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
#include "TParameter.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"

using namespace std;
using namespace PlotUtils;

double DistAvg(vector<double> vals) {
  vector<double> weights;
  double weightSum=0.0;
  for (auto val:vals){
    double distSum=0.0;
    for (auto val2:vals){
      distSum += fabs(val-val2);
    }
    double weight = (distSum) ? 100.0/distSum : -1.0;
    weights.push_back(weight);
    weightSum += weight;
  }
  double avg = 0.0;
  for (int i=0; i<vals.size(); ++i){
    avg += vals.at(i)*weights.at(i)/weightSum;
  }
  return avg;
}

TH1D* avgDistMax(vector<TFile*> files, TString name, string band){
  vector<MnvH1D*> mnvs;
  vector<TH1D*> syst;
  vector<double> POTs;
  double POTsum=0.0;

  for (auto file:files){
    mnvs.push_back((MnvH1D*)file->Get(name));
    syst.push_back(mnvs.back()->GetVertErrorBand(band)->GetHist(0));
    TParameter<double>* POT = (TParameter<double>*)file->Get("POTUsed");
    POTsum += POT->GetVal();
    POTs.push_back(POT->GetVal());
  }

  TH1D* hOut = new TH1D(mnvs.at(0)->GetCVHistoWithStatError());
  for (int i=0; i <= hOut->GetNbinsX(); ++i){
    vector<double> vals;
    for (int j=0; j<syst.size(); ++j) vals.push_back(syst.at(j)->GetBinContent(i)*POTsum/POTs.at(j));
    hOut->SetBinContent(i, DistAvg(vals));
  }

  return hOut;
}
