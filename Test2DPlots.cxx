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

int main(int argc, char* argv[]){//(TString mcName, TString dataName, TString outDir){

  if (argc != 5){
    cout << "Not enough args." << endl;
    return 17;
  }

  TString mcName = argv[1];
  TString dataName = argv[2];
  TString material = argv[3];
  bool doFitScale = (bool)(atoi(argv[4]));

  gStyle->SetOptStat(0);

  //TString mcName = "runEventLoopMC_SkippedSyst_MnvTuneV1_FVregion_Tracker_noNeutCuts.root";
  //TString mcName = "runEventLoopMC_MnvTuneV1_FVregion_SingleTarget_Fe_noNeutCuts_noVtx.root";
    //TString dataName = "runEventLoopData_SkippedSyst_MnvTuneV1_FVregion_Tracker_noNeutCuts.root";
  //TString dataName = "runEventLoopData_MnvTuneV1_FVregion_SingleTarget_Fe_noNeutCuts_noVtx.root";

  TString outDir = "Plots2D_"+material+"/";

  TFile* mcFile = new TFile(mcName,"READ");
  TFile* dataFile = new TFile(dataName,"READ");

  double mcPOT = ((TParameter<double>*)mcFile->Get("POTUsed"))->GetVal();
  double dataPOT = ((TParameter<double>*)dataFile->Get("POTUsed"))->GetVal();

  double scale = dataPOT/mcPOT;

  TString name = "TwoD/TwoD_pmu2D_"+material+"_selected_signal_reco";
  TString nameQE = "TwoD/TwoD_pmu2D_"+material+"_sig_IntType_QE";

  if (material == "Tracker"){
    name = "TwoD/TwoD_pmu2D_selected_signal_reco";
    nameQE = "TwoD/TwoD_pmu2D_sig_IntType_QE";
  }
  
  TString nameToSave = outDir+"pmu2D_projections_"+material;

  Draw2DBKGCategStack(name.Data(), mcFile, dataFile, material, scale, nameToSave);
  //Draw2DIntTypeStack(nameQE, mcFile, dataFile, material, scale, nameToSave);

  if (doFitScale){
    TString scaleName = "test_Fe/fitResults_6_252D.root";
    TString scaleName2 = "../onlyVtx/testVtx_Fe/USvtxFitResults2D.root";
    TString scaleName3 = "../onlyVtx/testVtx_Fe/DSvtxFitResults2D.root";
    TString scalerName = "TwoD_pmu2D_Fe_recoilE_fit1B_low_6_hi_25_bkg_background_BKG";
    TString scalerName2 = "TwoD_pmu2D_Fe_NPlanesUS_fitUS_low_5_hi_8_Plastic";
    TString scalerName3 = "TwoD_pmu2D_Fe_NPlanesDS_fitDS_low_2_hi_5_Plastic";

    TFile* scaleFile = new TFile(scaleName,"READ");
    TFile* scaleFile2 = new TFile(scaleName2,"READ");
    TFile* scaleFile3 = new TFile(scaleName3,"READ");
    MnvH2D* scaler = (MnvH2D*)scaleFile->Get(scalerName);
    MnvH2D* scaler2 = (MnvH2D*)scaleFile2->Get(scalerName2);
    MnvH2D* scaler3 = (MnvH2D*)scaleFile3->Get(scalerName3);

    Draw2DBKGCategStackWithBKGScale(name.Data(), mcFile, dataFile, material, scale, nameToSave+"_BKGScaled", scaler);
    Draw2DPlotInGridStatOnly(scaler, "recoilFit_Fe_Scale_yaxis", "Muon Transverse Momentum [GeV/c]","Scale Factor", "p_{//}");
    Draw2DPlotInGridStatOnly(scaler2, "USFit_Fe_Scale_yaxis", "Muon Transverse Momentum [GeV/c]","Scale Factor", "p_{//}");
    Draw2DPlotInGridStatOnly(scaler3, "DSFit_Fe_Scale_yaxis", "Muon Transverse Momentum [GeV/c]","Scale Factor", "p_{//}");
  }

  return 0;
}
