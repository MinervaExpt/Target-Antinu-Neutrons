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

void MakeRenormsScript(TString fileName_NO, TString fileName_YES, TString playlist) {
  TFile* fileNO = TFile::Open(fileName_NO);
  TFile* fileYES = TFile::Open(fileName_YES);

  MnvH2D* sigQE_NO = (MnvH2D*)fileNO->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_sigQE");
  MnvH2D* sigQE_YES = (MnvH2D*)fileYES->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_sigQE");
  TH2D* h_sigQE_NO = new TH2D(sigQE_NO->GetCVHistoWithStatError());
  TH2D* h_sigQE_YES = new TH2D(sigQE_YES->GetCVHistoWithStatError());
  h_sigQE_YES->Rebin2D(10,2);
  h_sigQE_NO->Rebin2D(10,2);

  MnvH2D* sig2p2h_NO = (MnvH2D*)fileNO->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_sig2p2h");
  MnvH2D* sig2p2h_YES = (MnvH2D*)fileYES->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_sig2p2h");
  TH2D* h_sig2p2h_NO = new TH2D(sig2p2h_NO->GetCVHistoWithStatError());
  TH2D* h_sig2p2h_YES = new TH2D(sig2p2h_YES->GetCVHistoWithStatError());
  h_sig2p2h_YES->Rebin2D(10,2);
  h_sig2p2h_NO->Rebin2D(10,2);

  MnvH2D* sigOther_NO = (MnvH2D*)fileNO->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_sigOther");
  MnvH2D* sigOther_YES = (MnvH2D*)fileYES->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_sigOther");
  TH2D* h_sigOther_NO = new TH2D(sigOther_NO->GetCVHistoWithStatError());
  TH2D* h_sigOther_YES = new TH2D(sigOther_YES->GetCVHistoWithStatError());
  h_sigOther_YES->Rebin2D(10,2);
  h_sigOther_NO->Rebin2D(10,2);

  MnvH2D* TrackableProt_NO = (MnvH2D*)fileNO->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_TrackableProt");
  MnvH2D* TrackableProt_YES = (MnvH2D*)fileYES->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_TrackableProt");
  TH2D* h_TrackableProt_NO = new TH2D(TrackableProt_NO->GetCVHistoWithStatError());
  TH2D* h_TrackableProt_YES = new TH2D(TrackableProt_YES->GetCVHistoWithStatError());
  h_TrackableProt_YES->Rebin2D(10,2);
  h_TrackableProt_NO->Rebin2D(10,2);

  MnvH2D* SubThresh_NO = (MnvH2D*)fileNO->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_SubThresh");
  MnvH2D* SubThresh_YES = (MnvH2D*)fileYES->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_SubThresh");
  TH2D* h_SubThresh_NO = new TH2D(SubThresh_NO->GetCVHistoWithStatError());
  TH2D* h_SubThresh_YES = new TH2D(SubThresh_YES->GetCVHistoWithStatError());
  h_SubThresh_YES->Rebin2D(10,2);
  h_SubThresh_NO->Rebin2D(10,2);

  MnvH2D* OneChargePi_NO = (MnvH2D*)fileNO->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_1chargePi");
  MnvH2D* OneChargePi_YES = (MnvH2D*)fileYES->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_1chargePi");
  TH2D* h_1chargePi_NO = new TH2D(OneChargePi_NO->GetCVHistoWithStatError());
  TH2D* h_1chargePi_YES = new TH2D(OneChargePi_YES->GetCVHistoWithStatError());
  h_1chargePi_YES->Rebin2D(10,2);
  h_1chargePi_NO->Rebin2D(10,2);
  
  MnvH2D* OneNeutPi_NO = (MnvH2D*)fileNO->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_1neutPi");
  MnvH2D* OneNeutPi_YES = (MnvH2D*)fileYES->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_1neutPi");
  TH2D* h_1neutPi_NO = new TH2D(OneNeutPi_NO->GetCVHistoWithStatError());
  TH2D* h_1neutPi_YES = new TH2D(OneNeutPi_YES->GetCVHistoWithStatError());
  h_1neutPi_YES->Rebin2D(10,2);
  h_1neutPi_NO->Rebin2D(10,2);

  MnvH2D* NPi_NO = (MnvH2D*)fileNO->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_NPi");
  MnvH2D* NPi_YES = (MnvH2D*)fileYES->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_NPi");
  TH2D* h_NPi_NO = new TH2D(NPi_NO->GetCVHistoWithStatError());
  TH2D* h_NPi_YES = new TH2D(NPi_YES->GetCVHistoWithStatError());
  h_NPi_YES->Rebin2D(10,2);
  h_NPi_NO->Rebin2D(10,2);

  MnvH2D* Other_NO = (MnvH2D*)fileNO->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_Other");
  MnvH2D* Other_YES = (MnvH2D*)fileYES->Get("NeutronInelHists/TrueEvRate_pT_v_LeadTn_Other");
  TH2D* h_Other_NO = new TH2D(Other_NO->GetCVHistoWithStatError());
  TH2D* h_Other_YES = new TH2D(Other_YES->GetCVHistoWithStatError());
  h_Other_YES->Rebin2D(10,2);
  h_Other_NO->Rebin2D(10,2);

  TFile* file = new TFile("NeutronRenorm_"+playlist+".root","RECREATE");
  file->cd();
  TH2D* NRW_sigQE = (TH2D*)h_sigQE_YES->Clone("NeutronRenormWeight_sigQE");
  NRW_sigQE->Divide(h_sigQE_NO);
  NRW_sigQE->Write();

  TH2D* NRW_sig2p2h = (TH2D*)h_sig2p2h_YES->Clone("NeutronRenormWeight_sig2p2h");
  NRW_sig2p2h->Divide(h_sig2p2h_NO);
  NRW_sig2p2h->Write();

  TH2D* NRW_sigOther = (TH2D*)h_sigOther_YES->Clone("NeutronRenormWeight_sigOther");
  NRW_sigOther->Divide(h_sigOther_NO);
  NRW_sigOther->Write();

  TH2D* NRW_TrackableProt = (TH2D*)h_TrackableProt_YES->Clone("NeutronRenormWeight_TrackableProt");
  NRW_TrackableProt->Divide(h_TrackableProt_NO);
  NRW_TrackableProt->Write();
  
  TH2D* NRW_SubThresh = (TH2D*)h_SubThresh_YES->Clone("NeutronRenormWeight_SubThresh");
  NRW_SubThresh->Divide(h_SubThresh_NO);
  NRW_SubThresh->Write();
  
  TH2D* NRW_1chargePi = (TH2D*)h_1chargePi_YES->Clone("NeutronRenormWeight_1chargePi");
  NRW_1chargePi->Divide(h_1chargePi_NO);
  NRW_1chargePi->Write();

  TH2D* NRW_1neutPi = (TH2D*)h_1neutPi_YES->Clone("NeutronRenormWeight_1neutPi");
  NRW_1neutPi->Divide(h_1neutPi_NO);
  NRW_1neutPi->Write();

  TH2D* NRW_NPi = (TH2D*)h_NPi_YES->Clone("NeutronRenormWeight_NPi");
  NRW_NPi->Divide(h_NPi_NO);
  NRW_NPi->Write();

  TH2D* NRW_Other = (TH2D*)h_Other_YES->Clone("NeutronRenormWeight_Other");
  NRW_Other->Divide(h_Other_NO);
  NRW_Other->Write();

  file->Close();
}
