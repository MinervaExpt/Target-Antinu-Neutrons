//File: signalBKGStack.cxx
//Info: This is a script to run a loop over all MC breakdown plots in a single histos file and save nice plots from them. Primarily used for validation against older plots.
//
//Usage: signalBKGStack <mc_file> <data_file> <output_directory> <plot_label>
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

TCanvas* DrawBKGCateg(string name, TFile* mcFile, TFile* dataFile, TString sample, double scale){

  bool primPar = false;

  TString sampleName = sample;

  MnvH1D* h_Sig = (MnvH1D*)mcFile->Get((TString)name);
  h_Sig->Scale(scale);
  MnvH1D* mcSum = (MnvH1D*)h_Sig->Clone();
  h_Sig->SetLineColor(TColor::GetColor("#999933"));
  h_Sig->SetFillColor(TColor::GetColor("#999933"));

  cout << "Handling: " << name << endl;
  string title = (string)h_Sig->GetTitle();
  TString Xtitle = h_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_Sig->GetYaxis()->GetTitle();

  string name_bkg = name;
  name_bkg.erase(name_bkg.length()-21,name_bkg.length());

  MnvH1D* h_1PiC_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_1chargePi");
  h_1PiC_Bkg->Scale(scale);
  mcSum->Add(h_1PiC_Bkg);
  h_1PiC_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_1PiC_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));


  MnvH1D* h_1Pi0_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_1neutPi");
  h_1Pi0_Bkg->Scale(scale);
  mcSum->Add(h_1Pi0_Bkg);
  h_1Pi0_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_1Pi0_Bkg->SetFillColor(TColor::GetColor("#117733"));

  MnvH1D* h_NPi_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_NPi");
  h_NPi_Bkg->Scale(scale);
  mcSum->Add(h_NPi_Bkg);
  h_NPi_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_NPi_Bkg->SetFillColor(TColor::GetColor("#CC6677"));

  MnvH1D* h_Other_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_background_Other");
  h_Other_Bkg->Scale(scale);
  mcSum->Add(h_Other_Bkg);
  h_Other_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_data = (MnvH1D*)dataFile->Get((TString)name_bkg+"_data");
  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone();
  h_data->AddMissingErrorBandsAndFillWithCV(*h_Sig);

  THStack* h = new THStack();
  h->Add((TH1D*)h_Other_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_NPi_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_1Pi0_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_1PiC_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Sig->GetCVHistoWithError().Clone());

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  cout << "areaScale: " << areaScale << endl;

  h->Draw("hist");
  c1->Update();

  size_t pos = 0;
  if ((pos=name.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"Other");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.04);
    primPar = true;
  }

  h->SetTitle(sampleName);//+" "+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitle("Events");
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.75);
  if (!primPar)h->SetMaximum((dataHist->GetMaximum())*1.05);

  pos=0;
  if ((pos=name.find("_ENHitSB")) != string::npos){
    sampleName = "EM Blob E/NHit SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name.find("_NBlobsSB")) != string::npos){
    sampleName = "N EM Blob SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name.find("_MichelSB")) != string::npos){
    sampleName = "Michel SideBand " + sample;
    h->SetTitle(sampleName);
  }

  /*
  if (Xtitle.Contains("pmu")){
    h->GetXaxis()->SetTitle("p_{#mu} [GeV/c]");
    h->SetTitle("Muon Momentum "+sampleName);
  }

  if (Xtitle.Contains("vtxZ")){
    h->GetXaxis()->SetTitle("vtx. Z [mm]");
    h->SetTitle("Vertex Z "+sampleName);
  }
  */

  //h->Draw("hist");
  c1->Update();
  /* This is useful for debugging whether systematics actuall changed.
  TH1D* h_Tot = (TH1D*)h->GetStack()->Last()->Clone();
  h_Tot->SetLineColor(kRed);
  h_Tot->SetFillColorAlpha(kPink + 1, 0.4);
  h_Tot->Draw("E2 SAME");
  */
  dataHist->Draw("same");
  c1->Update();

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry(h_Sig,"Signal");
  leg->AddEntry(h_1PiC_Bkg,"single #pi^{#pm}");
  leg->AddEntry(h_1Pi0_Bkg,"single #pi^{0}");
  leg->AddEntry(h_NPi_Bkg,"N#pi");
  leg->AddEntry(h_Other_Bkg,"Other");

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  MnvH1D* ratio = (MnvH1D*)h_data->Clone();
  ratio->Divide(ratio,mcSum);

  TH1D* mcRatio = new TH1D(mcSum->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin, 1);
  }

  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(3);
  ratio->SetTitle("");
  //ratio->SetTitleSize(0);
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetTitleSize(0.05*areaScale);
  ratio->GetYaxis()->SetTitleOffset(0.75/areaScale);
  ratio->GetYaxis()->SetLabelSize(ratio->GetYaxis()->GetLabelSize()*areaScale);

  ratio->GetXaxis()->SetLabelSize(ratio->GetXaxis()->GetLabelSize()*areaScale);
  ratio->GetXaxis()->SetTitleSize(0.04*areaScale);
  ratio->SetMinimum(0.5);
  ratio->SetMaximum(1.5);
  ratio->Draw();

  mcRatio->SetLineColor(kRed);
  mcRatio->SetLineWidth(3);
  mcRatio->SetFillColorAlpha(kPink + 1, 0.4);
  mcRatio->Draw("E2 SAME");

  TH1D* straightLine = (TH1D*)mcRatio->Clone();
  straightLine->SetFillStyle(0);
  straightLine->Draw("HIST SAME");

  c1->Update();
  return c1;
}

TCanvas* DrawIntType(string name_QE, TFile* mcFile, TFile* dataFile, TString sample, double scale){

  bool primPar = false;

  TString sampleName = sample;

  MnvH1D* h_QE_Sig = (MnvH1D*)mcFile->Get((TString)name_QE);
  h_QE_Sig->Scale(scale);
  MnvH1D* mcSum = (MnvH1D*)h_QE_Sig->Clone();
  h_QE_Sig->SetLineColor(TColor::GetColor("#88CCEE"));
  h_QE_Sig->SetFillColor(TColor::GetColor("#88CCEE"));


  string name_sig = (string)h_QE_Sig->GetName();
  name_sig.erase(name_sig.length()-3,name_sig.length());
  string name_bkg = name_sig;
  name_bkg.erase(name_bkg.length()-12,name_bkg.length());

  cout << "Handling: " << name_sig << endl;
  string title = (string)h_QE_Sig->GetTitle();
  TString Xtitle = h_QE_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_QE_Sig->GetYaxis()->GetTitle();
  /*
  cout << title << endl;
  title.erase(0, 8);
  cout << title << endl;
  cout << "" << endl;
  */

  MnvH1D* h_RES_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_RES");
  h_RES_Sig->Scale(scale);
  mcSum->Add(h_RES_Sig);
  h_RES_Sig->SetLineColor(TColor::GetColor("#117733"));
  h_RES_Sig->SetFillColor(TColor::GetColor("#117733"));

  MnvH1D* h_DIS_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_DIS");
  h_DIS_Sig->Scale(scale);
  mcSum->Add(h_DIS_Sig);
  h_DIS_Sig->SetLineColor(TColor::GetColor("#CC6677"));
  h_DIS_Sig->SetFillColor(TColor::GetColor("#CC6677"));

  MnvH1D* h_2p2h_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_2p2h");
  h_2p2h_Sig->Scale(scale);
  mcSum->Add(h_2p2h_Sig);
  h_2p2h_Sig->SetLineColor(TColor::GetColor("#44AA99"));
  h_2p2h_Sig->SetFillColor(TColor::GetColor("#44AA99"));

  MnvH1D* h_Other_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_Other");
  h_Other_Sig->Scale(scale);
  mcSum->Add(h_Other_Sig);
  h_Other_Sig->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_QE_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_QE"); 
  h_QE_Bkg->Scale(scale);
  mcSum->Add(h_QE_Bkg);
  h_QE_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_QE_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_QE_Bkg->SetFillStyle(3444);

  MnvH1D* h_RES_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_RES");
  h_RES_Bkg->Scale(scale);
  mcSum->Add(h_RES_Bkg);
  h_RES_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_RES_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_RES_Bkg->SetFillStyle(3444);

  MnvH1D* h_DIS_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_DIS");
  h_DIS_Bkg->Scale(scale);
  mcSum->Add(h_DIS_Bkg);
  h_DIS_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_DIS_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_DIS_Bkg->SetFillStyle(3444);

  MnvH1D* h_2p2h_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_2p2h");
  h_2p2h_Bkg->Scale(scale);
  mcSum->Add(h_2p2h_Bkg);
  h_2p2h_Bkg->SetLineColor(TColor::GetColor("#44AA99"));
  h_2p2h_Bkg->SetFillColor(TColor::GetColor("#44AA99"));
  h_2p2h_Bkg->SetFillStyle(3444);

  MnvH1D* h_Other_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_IntType_Other");
  h_Other_Bkg->Scale(scale);
  mcSum->Add(h_Other_Bkg);
  h_Other_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_Other_Bkg->SetFillStyle(3444);

  MnvH1D* h_data = (MnvH1D*)dataFile->Get((TString)name_bkg+"_data");
  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone();
  h_data->AddMissingErrorBandsAndFillWithCV(*h_QE_Sig);

  THStack* h = new THStack();
  h->Add((TH1D*)h_Other_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_2p2h_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_DIS_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_RES_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_QE_Bkg->GetCVHistoWithError().Clone());

  h->Add((TH1D*)h_Other_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_2p2h_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_DIS_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_RES_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_QE_Sig->GetCVHistoWithError().Clone());

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  h->Draw("hist");
  c1->Update();

  size_t pos = 0;
  if ((pos=name_sig.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"Other");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.04);
    primPar = true;
  }

  h->SetTitle(sampleName);//+" "+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitle("Events");
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.75);
  if (!primPar) h->SetMaximum((dataHist->GetMaximum())*1.05);
  
  pos=0;
  if ((pos=name_sig.find("_ENHitSB")) != string::npos){
    sampleName = "EM Blob E/NHit SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_NBlobsSB")) != string::npos){
    sampleName = "N EM Blob SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_MichelSB")) != string::npos){
    sampleName = "Michel SideBand " + sample;
    h->SetTitle(sampleName);
  }

  /*
  if (Xtitle.Contains("pmu")){
    h->GetXaxis()->SetTitle("p_{#mu} [GeV/c]");
    h->SetTitle("Muon Momentum "+sampleName);
  }

  if (Xtitle.Contains("vtxZ")){
    h->GetXaxis()->SetTitle("vtx. Z [mm]");
    h->SetTitle("Vertex Z "+sampleName);
  }
  */

  //h->Draw("hist");
  c1->Update();
  /* This is useful for debugging whether systematics actuall changed.
  TH1D* h_Tot = (TH1D*)h->GetStack()->Last()->Clone();
  h_Tot->SetLineColor(kRed);
  h_Tot->SetFillColorAlpha(kPink + 1, 0.4);
  h_Tot->Draw("E2 SAME");
  */
  dataHist->Draw("same");
  c1->Update();

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  leg->SetNColumns(2);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry((TObject*)0,"","");

  leg->AddEntry(h_QE_Sig,"Sig. + QE");
  leg->AddEntry(h_QE_Bkg,"Bkg. + QE");

  leg->AddEntry(h_RES_Sig,"Sig. + RES");
  leg->AddEntry(h_RES_Bkg,"Bkg. + RES");

  leg->AddEntry(h_DIS_Sig,"Sig. + DIS");
  leg->AddEntry(h_DIS_Bkg,"Bkg. + DIS");

  leg->AddEntry(h_2p2h_Sig,"Sig. + 2p2h");
  leg->AddEntry(h_2p2h_Bkg,"Bkg. + 2p2h");

  leg->AddEntry(h_Other_Sig,"Sig. + Other");
  leg->AddEntry(h_Other_Bkg,"Bkg. + Other");

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  MnvH1D* ratio = (MnvH1D*)h_data->Clone();
  ratio->Divide(ratio,mcSum);

  TH1D* mcRatio = new TH1D(mcSum->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin, 1);
  }

  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(3);
  ratio->SetTitle("");
  //ratio->SetTitleSize(0);
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetTitleSize(0.05*areaScale);
  ratio->GetYaxis()->SetTitleOffset(0.75/areaScale);
  ratio->GetYaxis()->SetLabelSize(ratio->GetYaxis()->GetLabelSize()*areaScale);

  ratio->GetXaxis()->SetLabelSize(ratio->GetXaxis()->GetLabelSize()*areaScale);
  ratio->GetXaxis()->SetTitleSize(0.04*areaScale);
  ratio->SetMinimum(0.5);
  ratio->SetMaximum(1.5);
  ratio->Draw();

  mcRatio->SetLineColor(kRed);
  mcRatio->SetLineWidth(3);
  mcRatio->SetFillColorAlpha(kPink + 1, 0.4);
  mcRatio->Draw("E2 SAME");

  TH1D* straightLine = (TH1D*)mcRatio->Clone();
  straightLine->SetFillStyle(0);
  straightLine->Draw("HIST SAME");

  c1->Update();
  return c1;
}

TCanvas* DrawTargetType(string name_C, TFile* mcFile, TFile* dataFile, TString sample, double scale){

  bool primPar = false;

  TString sampleName = sample;

  MnvH1D* h_C_Sig = (MnvH1D*)mcFile->Get((TString)name_C);
  h_C_Sig->Scale(scale);
  MnvH1D* mcSum = (MnvH1D*)h_C_Sig->Clone();
  h_C_Sig->SetLineColor(TColor::GetColor("#88CCEE"));
  h_C_Sig->SetFillColor(TColor::GetColor("#88CCEE"));

  string name_sig = (string)h_C_Sig->GetName();
  name_sig.erase(name_sig.length()-2,name_sig.length());
  string name_bkg = name_sig;
  name_bkg.erase(name_bkg.length()-15,name_bkg.length());

  cout << "Handling: " << name_sig << endl;
  string title = (string)h_C_Sig->GetTitle();
  TString Xtitle = h_C_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_C_Sig->GetYaxis()->GetTitle();
  /*
  cout << title << endl;
  title.erase(0, 7);
  cout << title << endl;
  cout << "" << endl;
  */

  MnvH1D* h_Fe_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_Fe");
  h_Fe_Sig->Scale(scale);
  mcSum->Add(h_Fe_Sig);
  h_Fe_Sig->SetLineColor(TColor::GetColor("#882255"));
  h_Fe_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_Pb_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_Pb");
  h_Pb_Sig->Scale(scale);
  mcSum->Add(h_Pb_Sig);
  h_Pb_Sig->SetLineColor(TColor::GetColor("#117733"));
  h_Pb_Sig->SetFillColor(TColor::GetColor("#117733"));

  MnvH1D* h_O_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_O");
  h_O_Sig->Scale(scale);
  mcSum->Add(h_O_Sig);
  h_O_Sig->SetLineColor(TColor::GetColor("#332288"));
  h_O_Sig->SetFillColor(TColor::GetColor("#332288"));

  MnvH1D* h_H_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_H");
  h_H_Sig->Scale(scale);
  mcSum->Add(h_H_Sig);
  h_H_Sig->SetLineColor(TColor::GetColor("#DDCC77"));
  h_H_Sig->SetFillColor(TColor::GetColor("#DDCC77"));

  //h_Prot_Sig->SetFillColor(TColor::GetColor("#999933"));

  MnvH1D* h_Other_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_Other");
  h_Other_Sig->Scale(scale);
  mcSum->Add(h_Other_Sig);
  h_Other_Sig->SetLineColor(TColor::GetColor("#CC6677"));
  h_Other_Sig->SetFillColor(TColor::GetColor("#CC6677"));

  //h_None_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_C_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_C");
  h_C_Bkg->Scale(scale);
  mcSum->Add(h_C_Bkg);
  h_C_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_C_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_C_Bkg->SetFillStyle(3444);

  MnvH1D* h_Fe_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_Fe");
  h_Fe_Bkg->Scale(scale);
  mcSum->Add(h_Fe_Bkg);
  h_Fe_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_Fe_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_Fe_Bkg->SetFillStyle(3444);

  MnvH1D* h_Pb_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_Pb");
  h_Pb_Bkg->Scale(scale);
  mcSum->Add(h_Pb_Bkg);
  h_Pb_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_Pb_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_Pb_Bkg->SetFillStyle(3444);

  MnvH1D* h_O_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_O");
  h_O_Bkg->Scale(scale);
  mcSum->Add(h_O_Bkg);
  h_O_Bkg->SetLineColor(TColor::GetColor("#332288"));
  h_O_Bkg->SetFillColor(TColor::GetColor("#332288"));
  h_O_Bkg->SetFillStyle(3444);

  MnvH1D* h_H_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_H");
  h_H_Bkg->Scale(scale);
  mcSum->Add(h_H_Bkg);
  h_H_Bkg->SetLineColor(TColor::GetColor("#DDCC77"));
  h_H_Bkg->SetFillColor(TColor::GetColor("#DDCC77"));
  h_H_Bkg->SetFillStyle(3444);

  MnvH1D* h_Other_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_TargetType_Other");
  h_Other_Bkg->Scale(scale);
  mcSum->Add(h_Other_Bkg);
  h_Other_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_Other_Bkg->SetFillStyle(3444);

  MnvH1D* h_data = (MnvH1D*)dataFile->Get((TString)name_bkg+"_data");
  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone();
  h_data->AddMissingErrorBandsAndFillWithCV(*h_C_Sig);

  THStack* h = new THStack();
  h->Add((TH1D*)h_Other_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_H_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_O_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Pb_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Fe_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_C_Bkg->GetCVHistoWithError().Clone());

  h->Add((TH1D*)h_Other_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_H_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_O_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Pb_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Fe_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_C_Sig->GetCVHistoWithError().Clone());

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  h->Draw("hist");
  c1->Update();

  size_t pos=0;
  if ((pos=name_sig.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"None");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.04);
    primPar = true;
  }

  h->SetTitle(sampleName);//+" "+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitle("Events");
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.75);
  if (!primPar)h->SetMaximum((dataHist->GetMaximum())*1.05);
  
  pos=0;
  if ((pos=name_sig.find("_ENHitSB")) != string::npos){
    sampleName = "EM Blob E/NHit SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_NBlobsSB")) != string::npos){
    sampleName = "N EM Blob SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_MichelSB")) != string::npos){
    sampleName = "Michel SideBand " + sample;
    h->SetTitle(sampleName);
  }

  /*
  if (Xtitle.Contains("pmu")){
    h->GetXaxis()->SetTitle("p_{#mu} [GeV/c]");
    h->SetTitle("Muon Momentum "+sampleName);
  }

  if (Xtitle.Contains("vtxZ")){
    h->GetXaxis()->SetTitle("vtx. Z [mm]");
    h->SetTitle("Vertex Z "+sampleName);
  }
  */

  //h->Draw("hist");
  c1->Update();
  /* This is useful for debugging whether systematics actuall changed.
  TH1D* h_Tot = (TH1D*)h->GetStack()->Last()->Clone();
  h_Tot->SetLineColor(kRed);
  h_Tot->SetFillColorAlpha(kPink + 1, 0.4);
  h_Tot->Draw("E2 SAME");
  */
  dataHist->Draw("same");
  c1->Update();

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  leg->SetNColumns(2);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry((TObject*)0,"","");

  leg->AddEntry(h_C_Sig,"Sig. + C");
  leg->AddEntry(h_C_Bkg,"Bkg. + C");

  leg->AddEntry(h_Fe_Sig,"Sig. + Fe");
  leg->AddEntry(h_Fe_Bkg,"Bkg. + Fe");

  leg->AddEntry(h_Pb_Sig,"Sig. + Pb");  
  leg->AddEntry(h_Pb_Bkg,"Bkg. + Pb");
  
  leg->AddEntry(h_O_Sig,"Sig. + O");
  leg->AddEntry(h_O_Bkg,"Bkg. + O");

  leg->AddEntry(h_H_Sig,"Sig. + H");
  leg->AddEntry(h_H_Bkg,"Bkg. + H");

  leg->AddEntry(h_Other_Sig,"Sig. + Other");
  leg->AddEntry(h_Other_Bkg,"Bkg. + Other");

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  MnvH1D* ratio = (MnvH1D*)h_data->Clone();
  ratio->Divide(ratio,mcSum);

  TH1D* mcRatio = new TH1D(mcSum->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin, 1);
  }

  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(3);
  ratio->SetTitle("");
  //ratio->SetTitleSize(0);
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetTitleSize(0.05*areaScale);
  ratio->GetYaxis()->SetTitleOffset(0.75/areaScale);
  ratio->GetYaxis()->SetLabelSize(ratio->GetYaxis()->GetLabelSize()*areaScale);

  ratio->GetXaxis()->SetLabelSize(ratio->GetXaxis()->GetLabelSize()*areaScale);  
  ratio->GetXaxis()->SetTitleSize(0.04*areaScale);
  ratio->SetMinimum(0.5);
  ratio->SetMaximum(1.5);
  ratio->Draw();

  mcRatio->SetLineColor(kRed);
  mcRatio->SetLineWidth(3);
  mcRatio->SetFillColorAlpha(kPink + 1, 0.4);
  mcRatio->Draw("E2 SAME");

  TH1D* straightLine = (TH1D*)mcRatio->Clone();
  straightLine->SetFillStyle(0);
  straightLine->Draw("HIST SAME");

  c1->Update();
  return c1;
}

TCanvas* DrawLeadBlobType(string name_Neut, TFile* mcFile, TFile* dataFile, TString sample, double scale){

  bool primPar = false;

  TString sampleName = sample;

  MnvH1D* h_Neut_Sig = (MnvH1D*)mcFile->Get((TString)name_Neut);
  h_Neut_Sig->Scale(scale);
  MnvH1D* mcSum = (MnvH1D*)h_Neut_Sig->Clone();
  h_Neut_Sig->SetLineColor(TColor::GetColor("#88CCEE"));
  h_Neut_Sig->SetFillColor(TColor::GetColor("#88CCEE"));

  string name_sig = (string)h_Neut_Sig->GetName();
  name_sig.erase(name_sig.length()-5,name_sig.length());
  string name_bkg = name_sig;
  name_bkg.erase(name_bkg.length()-17,name_bkg.length());
  cout << "Handling: " << name_sig << endl;

  string title = (string)h_Neut_Sig->GetTitle();
  TString Xtitle = h_Neut_Sig->GetXaxis()->GetTitle();
  TString Ytitle = h_Neut_Sig->GetYaxis()->GetTitle();

  /*
  cout << title << endl;
  title.erase(0, 10);
  cout << title << endl;
  cout << "" << endl;
  */

  MnvH1D* h_Mu_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_mu");
  h_Mu_Sig->Scale(scale);
  mcSum->Add(h_Mu_Sig);
  h_Mu_Sig->SetLineColor(TColor::GetColor("#44AA99"));
  h_Mu_Sig->SetFillColor(TColor::GetColor("#44AA99"));

  MnvH1D* h_Pi0_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_pi0");
  h_Pi0_Sig->Scale(scale);
  mcSum->Add(h_Pi0_Sig);
  h_Pi0_Sig->SetLineColor(TColor::GetColor("#117733"));
  h_Pi0_Sig->SetFillColor(TColor::GetColor("#117733"));

  MnvH1D* h_PiM_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_pim");
  h_PiM_Sig->Scale(scale);
  mcSum->Add(h_PiM_Sig);
  h_PiM_Sig->SetLineColor(TColor::GetColor("#332288"));
  h_PiM_Sig->SetFillColor(TColor::GetColor("#332288"));

  MnvH1D* h_PiP_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_pip");
  h_PiP_Sig->Scale(scale);
  mcSum->Add(h_PiP_Sig);
  h_PiP_Sig->SetLineColor(TColor::GetColor("#DDCC77"));
  h_PiP_Sig->SetFillColor(TColor::GetColor("#DDCC77"));

  MnvH1D* h_Prot_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_prot");
  h_Prot_Sig->Scale(scale);
  mcSum->Add(h_Prot_Sig);
  h_Prot_Sig->SetLineColor(TColor::GetColor("#999933"));
  h_Prot_Sig->SetFillColor(TColor::GetColor("#999933"));

  MnvH1D* h_Other_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_Other");
  h_Other_Sig->Scale(scale);
  mcSum->Add(h_Other_Sig);
  h_Other_Sig->SetLineColor(TColor::GetColor("#CC6677"));
  h_Other_Sig->SetFillColor(TColor::GetColor("#CC6677"));

  MnvH1D* h_None_Sig = (MnvH1D*)mcFile->Get((TString)name_sig+"_None");
  h_None_Sig->Scale(scale);
  mcSum->Add(h_None_Sig);
  h_None_Sig->SetLineColor(TColor::GetColor("#882255"));
  h_None_Sig->SetFillColor(TColor::GetColor("#882255"));

  MnvH1D* h_Neut_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_neut");
  h_Neut_Bkg->Scale(scale);
  mcSum->Add(h_Neut_Bkg);
  h_Neut_Bkg->SetLineColor(TColor::GetColor("#88CCEE"));
  h_Neut_Bkg->SetFillColor(TColor::GetColor("#88CCEE"));
  h_Neut_Bkg->SetFillStyle(3444);

  MnvH1D* h_Mu_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_mu");
  h_Mu_Bkg->Scale(scale);
  mcSum->Add(h_Mu_Bkg);
  h_Mu_Bkg->SetLineColor(TColor::GetColor("#44AA99"));
  h_Mu_Bkg->SetFillColor(TColor::GetColor("#44AA99"));
  h_Mu_Bkg->SetFillStyle(3444);

  MnvH1D* h_Pi0_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_pi0");
  h_Pi0_Bkg->Scale(scale);
  mcSum->Add(h_Pi0_Bkg);
  h_Pi0_Bkg->SetLineColor(TColor::GetColor("#117733"));
  h_Pi0_Bkg->SetFillColor(TColor::GetColor("#117733"));
  h_Pi0_Bkg->SetFillStyle(3444);

  MnvH1D* h_PiM_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_pim");
  h_PiM_Bkg->Scale(scale);
  mcSum->Add(h_PiM_Bkg);
  h_PiM_Bkg->SetLineColor(TColor::GetColor("#332288"));
  h_PiM_Bkg->SetFillColor(TColor::GetColor("#332288"));
  h_PiM_Bkg->SetFillStyle(3444);

  MnvH1D* h_PiP_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_pip");
  h_PiP_Bkg->Scale(scale);
  mcSum->Add(h_PiP_Bkg);
  h_PiP_Bkg->SetLineColor(TColor::GetColor("#DDCC77"));
  h_PiP_Bkg->SetFillColor(TColor::GetColor("#DDCC77"));
  h_PiP_Bkg->SetFillStyle(3444);

  MnvH1D* h_Prot_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_prot");
  h_Prot_Bkg->Scale(scale);
  mcSum->Add(h_Prot_Bkg);
  h_Prot_Bkg->SetLineColor(TColor::GetColor("#999933"));
  h_Prot_Bkg->SetFillColor(TColor::GetColor("#999933"));
  h_Prot_Bkg->SetFillStyle(3444);

  MnvH1D* h_Other_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_Other");
  h_Other_Bkg->Scale(scale);
  mcSum->Add(h_Other_Bkg);
  h_Other_Bkg->SetLineColor(TColor::GetColor("#CC6677"));
  h_Other_Bkg->SetFillColor(TColor::GetColor("#CC6677"));
  h_Other_Bkg->SetFillStyle(3444);

  MnvH1D* h_None_Bkg = (MnvH1D*)mcFile->Get((TString)name_bkg+"_bkg_LeadBlobType_None");
  h_None_Bkg->Scale(scale);
  mcSum->Add(h_None_Bkg);
  h_None_Bkg->SetLineColor(TColor::GetColor("#882255"));
  h_None_Bkg->SetFillColor(TColor::GetColor("#882255"));
  h_None_Bkg->SetFillStyle(3444);

  MnvH1D* h_data = (MnvH1D*)dataFile->Get((TString)name_bkg+"_data");
  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone();
  h_data->AddMissingErrorBandsAndFillWithCV(*h_Neut_Sig);

  THStack* h = new THStack();
  h->Add((TH1D*)h_None_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Other_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Prot_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_PiP_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_PiM_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Pi0_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Mu_Bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Neut_Bkg->GetCVHistoWithError().Clone());

  h->Add((TH1D*)h_None_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Other_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Prot_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_PiP_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_PiM_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Pi0_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Mu_Sig->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_Neut_Sig->GetCVHistoWithError().Clone());

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  h->Draw("hist");
  c1->Update();

  //ToDo: Get the naming of the axes fixed to be what it needs to be/make it easier to automate. This will need to be accompanied by a change to the Variable class to get the labels correct there. Current changes temporary in the interest of making plots for a talk on Oct. 14, 2021 in the exclusives meeting.

  size_t pos=0;
  if ((pos=name_sig.find("_primary_parent")) != string::npos){
    h->GetXaxis()->SetBinLabel(1,"None");
    h->GetXaxis()->SetBinLabel(2,"");
    h->GetXaxis()->SetBinLabel(3,"n");
    h->GetXaxis()->SetBinLabel(4,"p");
    h->GetXaxis()->SetBinLabel(5,"#pi^{0}");
    h->GetXaxis()->SetBinLabel(6,"#pi^{+}");
    h->GetXaxis()->SetBinLabel(7,"#pi^{-}");
    h->GetXaxis()->SetBinLabel(8,"#gamma");
    h->GetXaxis()->SetBinLabel(9,"e");
    h->GetXaxis()->SetBinLabel(10,"#mu");
    h->GetXaxis()->SetLabelSize(0.06);
    h->GetXaxis()->SetTitle("Blob Primary Parent");
    h->GetXaxis()->SetTitleSize(0.04);
    primPar = true;
  }

  h->SetTitle(sampleName);//+" "+title.c_str());
  h->GetXaxis()->SetTitle(Xtitle);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitle("Events");
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(0.75);
  if (!primPar)h->SetMaximum((dataHist->GetMaximum())*1.05);
  
  pos=0;
  if ((pos=name_sig.find("_ENHitSB")) != string::npos){
    sampleName = "EM Blob E/NHit SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_NBlobsSB")) != string::npos){
    sampleName = "N EM Blob SideBand " + sample;
    h->SetTitle(sampleName);
  }

  pos=0;
  if ((pos=name_sig.find("_MichelSB")) != string::npos){
    sampleName = "Michel SideBand " + sample;
    h->SetTitle(sampleName);
  }

  /*
  if (Xtitle.Contains("pmu")){
    h->GetXaxis()->SetTitle("p_{#mu} [GeV/c]");
    h->SetTitle("Muon Momentum "+sampleName);
  }

  if (Xtitle.Contains("vtxZ")){
    h->GetXaxis()->SetTitle("vtx. Z [mm]");
    h->SetTitle("Vertex Z "+sampleName);
  }
  */

  //  h->Draw("hist");
  c1->Update();
  /* This is useful for debugging whether systematics actuall changed.
  TH1D* h_Tot = (TH1D*)h->GetStack()->Last()->Clone();
  h_Tot->SetLineColor(kRed);
  h_Tot->SetFillColorAlpha(kPink + 1, 0.4);
  h_Tot->Draw("E2 SAME");
  */
  dataHist->Draw("same");
  c1->Update();

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  leg->SetNColumns(2);

  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry((TObject*)0,"","");

  leg->AddEntry(h_Neut_Sig,"Sig. + n");
  leg->AddEntry(h_Neut_Bkg,"Bkg. + n");

  leg->AddEntry(h_Mu_Sig,"Sig. + #mu");
  leg->AddEntry(h_Mu_Bkg,"Bkg. + #mu");

  leg->AddEntry(h_Pi0_Sig,"Sig. + #pi^{0}");  
  leg->AddEntry(h_Pi0_Bkg,"Bkg. + #pi^{0}");
  
  leg->AddEntry(h_PiM_Sig,"Sig. + #pi^{-}");
  leg->AddEntry(h_PiM_Bkg,"Bkg. + #pi^{-}");

  leg->AddEntry(h_PiP_Sig,"Sig. + #pi^{+}");
  leg->AddEntry(h_PiP_Bkg,"Bkg. + #pi^{+}");

  leg->AddEntry(h_Prot_Sig,"Sig. + p");
  leg->AddEntry(h_Prot_Bkg,"Bkg. + p");

  leg->AddEntry(h_Other_Sig,"Sig. + Other");
  leg->AddEntry(h_Other_Bkg,"Bkg. + Other");

  leg->AddEntry(h_None_Sig,"Sig. + None");
  leg->AddEntry(h_None_Bkg,"Bkg. + None");

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  MnvH1D* ratio = (MnvH1D*)h_data->Clone();
  ratio->Divide(ratio,mcSum);

  TH1D* mcRatio = new TH1D(mcSum->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin, 1);
  }

  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(3);
  ratio->SetTitle("");
  //ratio->SetTitleSize(0);
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetTitleSize(0.05*areaScale);
  ratio->GetYaxis()->SetTitleOffset(0.75/areaScale);
  ratio->GetYaxis()->SetLabelSize(ratio->GetYaxis()->GetLabelSize()*areaScale);

  ratio->GetXaxis()->SetLabelSize(ratio->GetXaxis()->GetLabelSize()*areaScale);
  ratio->GetXaxis()->SetTitleSize(0.04*areaScale);
  ratio->SetMinimum(0.5);
  ratio->SetMaximum(1.5);
  ratio->Draw();

  mcRatio->SetLineColor(kRed);
  mcRatio->SetLineWidth(3);
  mcRatio->SetFillColorAlpha(kPink + 1, 0.4);
  mcRatio->Draw("E2 SAME");

  TH1D* straightLine = (TH1D*)mcRatio->Clone();
  straightLine->SetFillStyle(0);
  straightLine->Draw("HIST SAME");

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

  string MCfileName=string(argv[1]);
  string DATAfileName=string(argv[2]);
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
  string fileNameStub = MCfileName;
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
    cout << "MC Input need be .root file." << endl;
    return 4;
  }

  cout << "Input MC file name parsed to: " << fileNameStub << endl;

  rootExt = ".root";
  slash = "/";
  token = "";
  fileNameStub = DATAfileName;
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
    cout << "DATA Input need be .root file." << endl;
    return 5;
  }

  cout << "Input Data file name parsed to: " << fileNameStub << endl;

  //cout << "Setting up MnvPlotter" << endl;
  //MnvPlotter* plotter = new MnvPlotter(kCCQEAntiNuStyle);

  TFile* mcFile = new TFile(MCfileName.c_str(),"READ");
  TFile* dataFile = new TFile(DATAfileName.c_str(),"READ");

  TParameter<double>* mcPOT = (TParameter<double>*)mcFile->Get("POTUsed");
  TParameter<double>* dataPOT = (TParameter<double>*)dataFile->Get("POTUsed");

  double scale = dataPOT->GetVal()/mcPOT->GetVal();
  cout << "POT scale factor: " << scale << endl;

  TList* keyList = mcFile->GetListOfKeys();
  if (!keyList){
    cout << "List of keys failed to get." << endl;
    return 5;
  }

  TIter next(keyList);
  TKey* key;
  while ( key = (TKey*)next() ){
    //cout << key->GetName() << endl;
    pos=0;
    string name=(string)key->GetName();
    if((pos=name.find("TwoD")) != string::npos) continue;
    else if((pos = name.find("_sig_IntType_QE")) != string::npos){
      TCanvas* c1 = DrawIntType(name,mcFile,dataFile,label,scale);
      TPad* top = (TPad*)c1->GetPrimitive("Overlay");
      name.erase(name.length()-15,name.length());
      c1->Print((TString)outDir+(TString)name+"_IntType_stacked.pdf");
      c1->Print((TString)outDir+(TString)name+"_IntType_stacked.png");
      top->SetLogy();
      c1->Update();
      c1->Print((TString)outDir+(TString)name+"_IntType_stacked_log.pdf");
      c1->Print((TString)outDir+(TString)name+"_IntType_stacked_log.png");     
      cout << "" << endl;
      delete c1;
    }
    else if ((pos = name.find("_sig_TargetType_C")) != string::npos){
      TCanvas* c1 = DrawTargetType(name,mcFile,dataFile,label,scale);
      TPad* top = (TPad*)c1->GetPrimitive("Overlay");
      name.erase(name.length()-17,name.length());
      c1->Print((TString)outDir+(TString)name+"_TargetType_stacked.pdf");
      c1->Print((TString)outDir+(TString)name+"_TargetType_stacked.png");
      top->SetLogy();
      c1->Update();
      c1->Print((TString)outDir+(TString)name+"_TargetType_stacked_log.pdf");
      c1->Print((TString)outDir+(TString)name+"_TargetType_stacked_log.png");           cout << "" << endl;
      delete c1;
    }
    else if ((pos = name.find("_sig_LeadBlobType_neut")) != string::npos){
      TCanvas* c1 = DrawLeadBlobType(name,mcFile,dataFile,label,scale);
      TPad* top = (TPad*)c1->GetPrimitive("Overlay");
      name.erase(name.length()-22,name.length());
      c1->Print((TString)outDir+(TString)name+"_LeadBlobType_stacked.pdf");
      c1->Print((TString)outDir+(TString)name+"_LeadBlobType_stacked.png");
      top->SetLogy();
      c1->Update();
      c1->Print((TString)outDir+(TString)name+"_LeadBlobType_stacked_log.pdf");
      c1->Print((TString)outDir+(TString)name+"_LeadBlobType_stacked_log.png");           cout << "" << endl;
      delete c1;
    }
    else if ((pos = name.find("_selected_signal_reco")) != string::npos){
      pos = 0;
      //Maybe can remove the following continues...
      //if ((pos = name.find("EMnBlobs_")) != string::npos) continue;
      //else if ((pos = name.find("EMBlobE_")) != string::npos) continue;
      //else if ((pos = name.find("EMBlobNHit_")) != string::npos) continue;
      //else if ((pos = name.find("EMBlobENHitRatio_")) != string::npos) continue;
      TCanvas* c1 = DrawBKGCateg(name,mcFile,dataFile,label,scale);
      TPad* top = (TPad*)c1->GetPrimitive("Overlay");
      name.erase(name.length()-21,name.length());
      c1->Print((TString)outDir+(TString)name+"_BKG_stacked.pdf");
      c1->Print((TString)outDir+(TString)name+"_BKG_stacked.png");
      top->SetLogy();
      c1->Update();
      c1->Print((TString)outDir+(TString)name+"_BKG_stacked_log.pdf");
      c1->Print((TString)outDir+(TString)name+"_BKG_stacked_log.png");           cout << "" << endl;
      delete c1;
    }
    /*
    else if ((pos = name.find("_data")) != string::npos){
      pos = 0;
      if ((pos=name.find("SB")) != string::npos) continue;
      cout << "Plotting error summary for: " << name << endl;
      TCanvas* c1 = new TCanvas("c1","c1",1200,800);
      MnvH1D* h_mc_data = (MnvH1D*)mcFile->Get((TString)name);
      name.erase(name.length()-5,name.length());
      plotter->DrawErrorSummary(h_mc_data);
      c1->Print((TString)outDir+(TString)name+"_err_summary.pdf");
      c1->Print((TString)outDir+(TString)name+"_err_summary.png");
      cout << "" << endl;
      delete c1;
    }
    */
  }

  //cout << "Deleting the MnvPlotter." << endl;
  //delete plotter;

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  mcFile->Close();
  dataFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
