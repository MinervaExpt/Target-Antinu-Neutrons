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

void Plot(PlotUtils::MnvH1D& hist)
{
  TString prefix = hist.GetName();
  
  TCanvas c1("c1");

  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kCCQEAntiNuStyle);
  plotter.axis_maximum = 0.4;

  plotter.DrawErrorSummary(&hist);
  c1.Print(prefix  + "_uncertaintySummary.png");
  c1.Print(prefix  + "_uncertaintySummary.pdf");
  c1.Print(prefix  + "_uncertaintySummary.C");

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Neutron Interactions");
  c1.Print(prefix  + "_neutronUncertainties.png");
  c1.Print(prefix  + "_neutronUncertainties.pdf");
  c1.Print(prefix  + "_neutronUncertainties.C");

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Cross Section Models");
  c1.Print(prefix  + "_xSecUncertainties.png");
  c1.Print(prefix  + "_xSecUncertainties.pdf");
  c1.Print(prefix  + "_xSecUncertainties.C");

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Calorimetric Response");
  c1.Print(prefix  + "_RecoilUncertainties.png");
  c1.Print(prefix  + "_RecoilUncertainties.pdf");
  c1.Print(prefix  + "_RecoilUncertainties.C");

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "FSI Models");
  c1.Print(prefix  + "_FSIUncertainties.png");
  c1.Print(prefix  + "_FSIUncertainties.pdf");
  c1.Print(prefix  + "_FSIUncertainties.C");

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Muon Reconstruction");
  c1.Print(prefix  + "_MuonUncertainties.png");
  c1.Print(prefix  + "_MuonUncertainties.pdf");
  c1.Print(prefix  + "_MuonUncertainties.C");

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "GEANT4 Charged");
  c1.Print(prefix  + "_GEANTUncertainties.png");
  c1.Print(prefix  + "_GEANTUncertainties.pdf");
  c1.Print(prefix  + "_GEANTUncertainties.C");

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Flux & Norm.");
  c1.Print(prefix  + "_NormalizationUncertainties.png");
  c1.Print(prefix  + "_NormalizationUncertainties.pdf");
  c1.Print(prefix  + "_NormalizationUncertainties.C");

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "MnvTune V2");
  c1.Print(prefix  + "_MnvTuneUncertainties.png");
  c1.Print(prefix  + "_MnvTuneUncertainties.pdf");
  c1.Print(prefix  + "_MnvTuneUncertainties.C");
}

void xSecErrSummaries(TString fileName) {

  gROOT->SetBatch();

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
    if ((TString)(key->GetClassName()) != "PlotUtils::MnvH1D") continue;
    MnvH1D* h = (MnvH1D*)inFile->Get((TString)(key->GetName()));
    Plot(*h);
  }

  cout << "HEY YOU DID IT!!!" << endl;
  return;
}
