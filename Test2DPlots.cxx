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

//David includes
#include "plotTools/CategoryPlots.h"

using namespace std;
using namespace PlotUtils;

int main(){//(TString mcName, TString dataName, TString outDir){

  TString mcName = "runEventLoopMC_SkippedSyst_MnvTuneV1_FVregion_Tracker_noNeut.root";
  TString dataName = "runEventLoopData_SkippedSyst_MnvTuneV1_FVregion_Tracker_noNeut.root";

  TString outDir = "TEST/";

  TFile* mcFile = new TFile(mcName,"READ");
  TFile* dataFile = new TFile(dataName,"READ");

  double mcPOT = ((TParameter<double>*)mcFile->Get("POTUsed"))->GetVal();
  double dataPOT = ((TParameter<double>*)dataFile->Get("POTUsed"))->GetVal();

  double scale = dataPOT/mcPOT;

  string name = "TwoD/TwoD_pmu2D_selected_signal_reco";
  TString nameToSave = outDir+"test_pmu2D_yaxis_projections";

  TESTDraw2DBKGCategLines(name, mcFile, dataFile, "", scale, nameToSave);

  return 0;
}
