//File: BKGFittingOverhaul.cxx
//Info: This script is intended to fit recoil/pTmu plots using TMinuit primarily for the neutron selected sample.
//
//Usage: BKGFitting <mc_file> <data_file> <outdir> <outFileTag> <recoilE/pTmu/vtxZ> <doSyst (only 0 means no)> <material> <breakdownInnerPlastic> optional: <mainTagName> <lowFitBinNum> <hiFitBinNum> <do fits in bins of muon momentum (only 0 means no)> TODO: Save the information beyond just printing it out
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

//TODO: Same SegFault Business From My Plotting Code... I'm assuming I just need to delete things carefully that I'm not yet.

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
#include "TObjString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"
#include "TParameter.h"
#include "TFractionFitter.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

using namespace std;
using namespace PlotUtils;

bool PathExists(string path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

void syncCVHistos1D(MnvH1D* hist){
  TH1D theCVHisto(*hist);
  theCVHisto.SetDirectory(0);
  auto bandNames = hist->GetErrorBandNames();

  for (auto bandName : bandNames){
    auto band = hist->GetVertErrorBand(bandName);
    theCVHisto.Copy(*band);
  }
}

double calcVarianceSig(double Ds, double Db, double Ss, double Sb, double Bs, double Bb, double errDs, double errDb, double errSs, double errSb, double errBs, double errBb){
  if (errDs==0 && errDb==0 && errSs==0 && errSb==0 && errBs==0 && errBb==0) return 0.0;
  
  double num=Ds*Bb-Db*Bs;
  double den=Ss*Bb-Sb*Bs;
  double varDs=errDs*Bb/den;
  varDs *= varDs;
  double varDb=-errDb*Bs/den;
  varDb *= varDb;
  double varSs=(-errSs*Bb*num)/(den*den);
  varSs *= varSs;
  double varSb=(errSb*Bs*num)/(den*den);
  varSb *= varSb;
  double varBs=((-Db/den)+((Sb*num)/(den*den)))*errBs;
  varBs *= varBs;
  double varBb=((Ds/den)+((-Ss*num)/(den*den)))*errBb;
  varBb *= varBb;
  return varDs+varDb+varSs+varSb+varBs+varBb;
}

double calcVarianceBKG(double Ds, double Db, double Ss, double Sb, double Bs, double Bb, double errDs, double errDb, double errSs, double errSb, double errBs, double errBb){
  if (errDs==0 && errDb==0 && errSs==0 && errSb==0 && errBs==0 && errBb==0) return 0.0;

  double num=Db*Ss-Ds*Sb;
  double den=Ss*Bb-Sb*Bs;
  double varDs=-errDs*Sb/den;
  varDs *= varDs;
  double varDb=errDb*Ss/den;
  varDb *= varDb;
  double varBs=(errBs*Sb*num)/(den*den);
  varBs *= varBs;
  double varBb=(-errBb*Ss*num)/(den*den);
  varBb *= varBb;
  double varSs=((Db/den)+((-Bb*num)/(den*den)))*errSs;
  varSs *= varSs;
  double varSb=((-Ds/den)+((Bs*num)/(den*den)))*errSb;
  varSb *= varSb;
  return varDs+varDb+varSs+varSb+varBs+varBb;
}

map<TString, MnvH1D*> PerformTwoFunctionSolution(map<TString, vector<MnvH1D*>> mcFitHistsINPUT, vector<MnvH1D*> dataFitHistsINPUT, vector<vector<MnvH1D*>> nonFitHists, TString varName, map<TString, vector<TString>> tagsToSave){
  map<TString, MnvH1D*> scalers;

  map<TString, vector<MnvH1D*>> mcFitHists;
  for (auto categ:mcFitHistsINPUT){
    for (unsigned int iHist=0; iHist < categ.second.size(); ++iHist){
      mcFitHists[categ.first].push_back((MnvH1D*)categ.second.at(iHist)->Clone());
    }
  }
    
  if (mcFitHists.size() != 2){
    cout << "Passed too many or too few MC histograms to the Solver." << endl;
    return scalers;
  }

  if (mcFitHists["BKG"].size() !=2 || mcFitHists["Signal"].size() != 2 || mcFitHists.size() != 2){
    cout << "Must pass signal and BKG to be fit in both regions." << endl;
    return scalers;
  }  

  vector<MnvH1D*> dataFitHists;
  for (unsigned int iHist=0; iHist < dataFitHistsINPUT.size(); ++iHist){
    dataFitHists.push_back((MnvH1D*)dataFitHistsINPUT.at(iHist)->Clone());
  }
  
  if (dataFitHists.size() != 2){
    cout << "There should be two data histograms." << endl;
    return scalers;
  }

  map<int, int> nonFitMap;
  nonFitMap[2] = 1;

  int nonFit = 0;
  
  for (auto hists:nonFitHists){
    nonFit = nonFitMap[hists.size()];
    if (nonFit == 0){
      cout << "There should be 2 nonFit histograms per unfit category." << endl;
      return scalers;
    }
  }
  if (nonFit==0 && nonFitHists.size() > 0){
    cout << "There were nonFit histograms which weren't the appropriate number of 2 per category." << endl;
    return scalers;
  }
  
  map<TString, TString> namesToSave;
  
  for (auto hists:mcFitHists){
    TString dumpTag = "";//tag which lets you know which histos to scale later.
    for (auto tag:tagsToSave[hists.first]){
      dumpTag = dumpTag+"_t_"+tag;
    }
    namesToSave[hists.first]= varName+"_fit_Exact"+dumpTag;
  }

  //Guide <D/S/B/_<i> is data(-nonFit BKGs), signal, background, respectively in region i (s for signal, b for sideband)
  MnvH1D* denom1 = (MnvH1D*)mcFitHists["Signal"].at(0)->Clone(); //S_s
  denom1->Multiply(denom1,mcFitHists["BKG"].at(1));//S_s*B_b
  MnvH1D* denom2 = (MnvH1D*)mcFitHists["Signal"].at(1)->Clone();//S_b
  denom2->Multiply(denom2,mcFitHists["BKG"].at(0));//S_b*B_s
  MnvH1D* denom = (MnvH1D*)denom1->Clone();
  denom->Add(denom2,-1.0);//S_s*B_b-S_b*B_s
  
  //Do this at the last second before manipulating these data hists with MC hists... If any later stage throws an error about missing error bands, then it raises a question about why denom and others don't have the same error bands.
  dataFitHists.at(0)->AddMissingErrorBandsAndFillWithCV(*denom);
  dataFitHists.at(1)->AddMissingErrorBandsAndFillWithCV(*denom);
  
  //Subtracting nonFit BKGS so that the data used is expected to be Signal + fit backgrounds
  for(auto hists:nonFitHists){
    dataFitHists.at(0)->Add(hists.at(0),-1.0);//D_s
    dataFitHists.at(1)->Add(hists.at(1),-1.0); //D_b
  }
  
  scalers["Signal"] = (MnvH1D*)dataFitHists.at(0)->Clone(namesToSave["Signal"]);//D_s
  
  scalers["BKG"] = (MnvH1D*)dataFitHists.at(1)->Clone(namesToSave["BKG"]);//D_b
  
  scalers["Signal"]->Multiply(scalers["Signal"],mcFitHists["BKG"].at(1));//D_s*B_b
  MnvH1D* sigSub = (MnvH1D*)dataFitHists.at(1)->Clone();//D_b
  sigSub->Multiply(sigSub,mcFitHists["BKG"].at(0));//D_b*B_s
  scalers["Signal"]->Add(sigSub,-1.0);//D_s*B_b-D_b*B_s
  scalers["Signal"]->Divide(scalers["Signal"],denom);//(D_s*B_b-D_b*B_s)/(S_s*B_b-S_b*B_s)

  scalers["BKG"]->Multiply(scalers["BKG"],mcFitHists["Signal"].at(0));//D_b*S_s
  MnvH1D* bkgSub = (MnvH1D*)dataFitHists.at(0)->Clone();//D_s
  bkgSub->Multiply(bkgSub,mcFitHists["Signal"].at(1));//D_s*S_b
  scalers["BKG"]->Add(bkgSub,-1.0);//D_b*S_s-D_s*S_b
  scalers["BKG"]->Divide(scalers["BKG"],denom);//(D_b*S_s-D_s*S_b)/(S_s*B_b-S_b*B_s)

  TH1D* resultBKG=(TH1D*)(scalers["BKG"]->GetCVHistoWithStatError().Clone());
  TH1D* resultSig=(TH1D*)(scalers["Signal"]->GetCVHistoWithStatError().Clone());

  TH1D* Ds=(TH1D*)dataFitHists.at(0)->GetCVHistoWithStatError().Clone();
  TH1D* Db=(TH1D*)dataFitHists.at(1)->GetCVHistoWithStatError().Clone();
  TH1D* Ss=(TH1D*)mcFitHists["Signal"].at(0)->GetCVHistoWithStatError().Clone();
  TH1D* Sb=(TH1D*)mcFitHists["Signal"].at(1)->GetCVHistoWithStatError().Clone();
  TH1D* Bs=(TH1D*)mcFitHists["BKG"].at(0)->GetCVHistoWithStatError().Clone();
  TH1D* Bb=(TH1D*)mcFitHists["BKG"].at(1)->GetCVHistoWithStatError().Clone();

  int nBins=(Ds->GetNbinsX()+1);  
  for (int iBin=0;iBin<=nBins;++iBin){
    double errSig=sqrt(calcVarianceSig(Ds->GetBinContent(iBin),Db->GetBinContent(iBin),Ss->GetBinContent(iBin),Sb->GetBinContent(iBin),Bs->GetBinContent(iBin),Bb->GetBinContent(iBin),Ds->GetBinError(iBin),Db->GetBinError(iBin),Ss->GetBinError(iBin),Sb->GetBinError(iBin),Bs->GetBinError(iBin),Bb->GetBinError(iBin)));
    double errBKG=sqrt(calcVarianceBKG(Ds->GetBinContent(iBin),Db->GetBinContent(iBin),Ss->GetBinContent(iBin),Sb->GetBinContent(iBin),Bs->GetBinContent(iBin),Bb->GetBinContent(iBin),Ds->GetBinError(iBin),Db->GetBinError(iBin),Ss->GetBinError(iBin),Sb->GetBinError(iBin),Bs->GetBinError(iBin),Bb->GetBinError(iBin)));
    scalers["Signal"]->SetBinError(iBin, errSig);
    scalers["BKG"]->SetBinError(iBin, errBKG);
  }

  syncCVHistos1D(scalers["Signal"]);
  syncCVHistos1D(scalers["BKG"]);
  
  return scalers;//This does nothing to correct for anything like non-solvable equations or negative factors. Want to see if the negative factors are even an issue. The underflow will be and I'm not sure exactly what it will do to the histogram if it tries to divide by the 0 content of 0*0-0*0 lol...
}

int main(int argc, char* argv[]) {

  gStyle->SetOptStat(0);

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc < 9 || argc > 13) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string MCfileName = string(argv[1]);
  string DATAfileName = string(argv[2]);
  string outDir = string(argv[3]);
  string outFileTag = string(argv[4]);
  TString varName= argv[5];
  bool doSyst = (bool)(atoi(argv[6]));
  TString material = argv[7];
  TString materialTag = material;
  if (materialTag.Contains("LINE")){
    if (materialTag.Contains("Carbon")){
      material="C";
    }
    else if (materialTag.Contains("Iron")){
      material="Fe";
    }
    else if (materialTag.Contains("Lead")){
      material="Pb";
    }
    else if (materialTag.Contains("Water")){
      material="Water";
    }
    else material = "";
  }
  material = (material != "") ? "_"+material : material;
  
  if (materialTag == ""){
    cout << "Fitting Tracker" << endl;
  }
  else {
    cout << "Fitting " << materialTag << endl;
  }

  bool breakInner = (bool)(atoi(argv[8]));
  int fitMuonBins = 0;
  
  //vector<TString> namesToSave = {"pTmu","recoilE"};
  vector<TString> namesToSave = {varName};//, "NPlanes"};
  //vector<TString> namesToSave = {"pTmu","recoilE","NPlanes"};
  //vector<TString> namesToSave = {};

  int lowBin = 1;//Will be truncated later. Lowest allowed value for pT or recoil fits.
  int hiBin = 50;//Will be truncated later. Highest allowed value for pT or recoil fits.
  TString mainTag = "";
  if (argc > 10) lowBin = max(atoi(argv[10]),1);//Not allowed lower than 1
  if (argc > 11) hiBin = min(50, atoi(argv[11]));//Not allowed higher than 50
  if (argc > 12) fitMuonBins = atoi(argv[12]);

  string rootExt = ".root";
  string slash = "/";
  string token;
  string fileNameStub = MCfileName;
  size_t pos=0;

  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory: " << outDir << " doesn't exist. Exiting" << endl;
    return 3;
  }

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

  if (varName.Contains("pTmu")){
    if (lowBin >= 9){
      cout << "Not a valid fitting range for pTmu. Maximum bin is 14." << endl;
      return 100;
    }
    hiBin = min(9,hiBin);
    mainTag = "_RecoilSB";
  }
  else if (varName == "recoilE"){
    if (argc < 9) lowBin = 26;
    mainTag = "_PreRecoilCut";
  }
  else {
    cout << "Not a valid variable name to fit." << endl;
    return 111;
  }
  if (argc > 8) mainTag = argv[9];  

  //NEED TO CODE IN SOMETHING THAT HANDLES THE VERTEX PLOTS IN THE TARGETS.

  //cout << "Setting up MnvPlotter" << endl;
  //MnvPlotter* plotter = new MnvPlotter(kCCQEAntiNuStyle);

  TFile* outFile = new TFile((outDir+"fitResults_"+outFileTag+".root").c_str(),"RECREATE");
  TFile* mcFile = new TFile(MCfileName.c_str(),"READ");
  TFile* dataFile = new TFile(DATAfileName.c_str(),"READ");

  TParameter<double>* mcPOT = (TParameter<double>*)mcFile->Get("POTUsed");
  TParameter<double>* dataPOT = (TParameter<double>*)dataFile->Get("POTUsed");

  vector<TString> tags = {mainTag};

  /*
  if (Tgts){
    tags.push_back(mainTag+"_Tgt1");
    tags.push_back(mainTag+"_Tgt2");
    tags.push_back(mainTag+"_Tgt3");
    tags.push_back(mainTag+"_Tgt4");
    tags.push_back(mainTag+"_Tgt5");
  }
  */

  if (fitMuonBins){
    tags.push_back((TString)("_bin_lost"));
    for (int iBin=0; iBin < 14; ++iBin){
      TString tag = "_bin_"+to_string(iBin);
      //cout << tag << endl;
      tags.push_back(tag);
    }
  }

  double POTscale = dataPOT->GetVal()/mcPOT->GetVal();
  //cout << "POT scale factor: " << scale << endl;
  
  //TString varName = "recoilE";

  for (int iTag=0; iTag < tags.size(); ++iTag){

    TString tag = tags.at(iTag);
    TString name = varName+tag+material;
    TString grabName, grabName_noTag;
    if (tag.Contains("_Tgt")){
      TObjArray* tagArr = tag.Tokenize("Tgt");
      grabName="ByTgt_Tgt"+((TObjString*)(tagArr->At(tagArr->GetEntries()-1)))->String()+"/"+name;
      grabName_noTag = "ByTgt_Tgt"+((TObjString*)(tagArr->At(tagArr->GetEntries()-1)))->String()+"/"+varName;
    }
    else{
      grabName = name;
      grabName_noTag = varName+material;
    }

    cout << "Performing Fitting and Scaling for: " << name << endl;
    cout << "Grabbing: " << grabName << endl;
    cout << "TESTING NO TAG: " << grabName_noTag << endl;

    //TODO:
    // 
    double coarseVtxBinsUS[11];//Modified to be a number of planes instead of vertex Z values.                                             
    for (int iBin=0; iBin<11;++iBin) coarseVtxBinsUS[iBin] = (double)iBin-10.5;

    double coarseVtxBinsDS[11];//Modified to be a number of planes instead of vertex Z values.                                             
    for (int iBin=11; iBin<22;++iBin) coarseVtxBinsDS[iBin-11] = (double)iBin-10.5;

    double binsNPlanes[22];
    for (unsigned int i=0; i < 22;++i) binsNPlanes[i]=1.0*i-10.5;

    map<TString,MnvH1D*> varsToSave = {};
    for (auto nameSave:namesToSave){
      //TODO:
      // 
      if (nameSave == "NPlanes"){
	MnvH1D* hNPlanesUS = new MnvH1D("vtxZ_ByTgt_US"+tag,"",10,coarseVtxBinsUS);
	varsToSave[nameSave+"US"+tag] = hNPlanesUS->Clone();
	MnvH1D* hNPlanesDS = new MnvH1D(nameSave+"DS"+tag,"",10,coarseVtxBinsDS);
	varsToSave[nameSave+"DS"+tag] = hNPlanesDS->Clone();
      }
      //
      else if (tag.Contains("_Tgt")){
	TObjArray* tagArr = tag.Tokenize("Tgt");
	TString grabNameSave = "ByTgt_Tgt"+((TObjString*)(tagArr->At(tagArr->GetEntries()-1)))->String()+"/"+nameSave+tag;
	varsToSave[nameSave+tag] = (MnvH1D*)(mcFile->Get(grabNameSave+"_selected_signal_reco"))->Clone();
      }
      else varsToSave[nameSave+tag] = (MnvH1D*)(mcFile->Get(nameSave+tag+material+"_selected_signal_reco"))->Clone();
    }

    cout << "DATA?" << endl;
    cout << "Name w/ tag: " << grabName+"_data" << endl;

    MnvH1D* dataHist = (MnvH1D*)(dataFile->Get(grabName+"_data"))->Clone();

    cout << "Name w/o tag: " << grabName_noTag+"_data" << endl;
    MnvH1D* dataHist_noTag = (MnvH1D*)(dataFile->Get(grabName_noTag+"_data"))->Clone();

    cout << "Sig?" << endl;

    MnvH1D* sigHist = (MnvH1D*)(mcFile->Get(grabName+"_selected_signal_reco"))->Clone();
    MnvH1D* sigHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_selected_signal_reco"))->Clone();
    sigHist->Scale(POTscale);
    sigHist_noTag->Scale(POTscale);

    cout << "Bkg?" << endl;

    //
    MnvH1D* chargePiHist = (MnvH1D*)(mcFile->Get(grabName+"_background_1chargePi"))->Clone();
    MnvH1D* chargePiHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_1chargePi"))->Clone();
    chargePiHist->Scale(POTscale);
    chargePiHist_noTag->Scale(POTscale);
    MnvH1D* neutPiHist = (MnvH1D*)(mcFile->Get(grabName+"_background_1neutPi"))->Clone();
    MnvH1D* neutPiHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_1neutPi"))->Clone();
    neutPiHist->Scale(POTscale);
    neutPiHist_noTag->Scale(POTscale);
    MnvH1D* NPiHist = (MnvH1D*)(mcFile->Get(grabName+"_background_NPi"))->Clone();
    MnvH1D* NPiHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_NPi"))->Clone();
    NPiHist->Scale(POTscale);
    NPiHist_noTag->Scale(POTscale);
    MnvH1D* otherHist = (MnvH1D*)(mcFile->Get(grabName+"_background_Other"))->Clone();
    MnvH1D* otherHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_Other"))->Clone();
    otherHist->Scale(POTscale);
    otherHist_noTag->Scale(POTscale);

    cout << "Wrong Nucl?" << endl;

    MnvH1D* wrongNuclHist = (MnvH1D*)(mcFile->Get(grabName+"_background_Wrong_Nucleus"))->Clone();
    MnvH1D* wrongNuclHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_Wrong_Nucleus"))->Clone();
    wrongNuclHist->Scale(POTscale);
    wrongNuclHist_noTag->Scale(POTscale);

    cout << "Material : " << material << endl;

    MnvH1D* USHist = (MnvH1D*)(mcFile->Get(grabName+"_background_USPlastic"))->Clone();
    MnvH1D* USHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_USPlastic"))->Clone();
    if (breakInner){
      cout << "Entering check!" << endl;
      USHist = (MnvH1D*)(mcFile->Get(varName+"_InnerUSPlastic"+tag+material+"_selected_signal_reco"))->Clone();
      USHist->Add((MnvH1D*)(mcFile->Get(varName+"_InnerUSPlastic"+tag+material+"_background_1chargePi")));
      USHist->Add((MnvH1D*)(mcFile->Get(varName+"_InnerUSPlastic"+tag+material+"_background_1neutPi")));
      USHist->Add((MnvH1D*)(mcFile->Get(varName+"_InnerUSPlastic"+tag+material+"_background_NPi")));
      USHist->Add((MnvH1D*)(mcFile->Get(varName+"_InnerUSPlastic"+tag+material+"_background_Other")));
      USHist_noTag = (MnvH1D*)(mcFile->Get(varName+"_InnerUSPlastic"+material+"_selected_signal_reco"))->Clone();
      USHist_noTag->Add((MnvH1D*)(mcFile->Get(varName+"_InnerUSPlastic"+material+"_background_1chargePi")));
      USHist_noTag->Add((MnvH1D*)(mcFile->Get(varName+"_InnerUSPlastic"+material+"_background_1neutPi")));
      USHist_noTag->Add((MnvH1D*)(mcFile->Get(varName+"_InnerUSPlastic"+material+"_background_NPi")));
      USHist_noTag->Add((MnvH1D*)(mcFile->Get(varName+"_InnerUSPlastic"+material+"_background_Other")));
    }
    USHist->Scale(POTscale);
    USHist_noTag->Scale(POTscale);

    cout << "Material : " << material << endl;

    MnvH1D* DSHist = (MnvH1D*)(mcFile->Get(grabName+"_background_DSPlastic"))->Clone();
    MnvH1D* DSHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_background_DSPlastic"))->Clone();
    if (breakInner){
      cout << "Entering Check!" << endl;
      DSHist = (MnvH1D*)(mcFile->Get(varName+"_InnerDSPlastic"+tag+material+"_selected_signal_reco"))->Clone();
      DSHist->Add((MnvH1D*)(mcFile->Get(varName+"_InnerDSPlastic"+tag+material+"_background_1chargePi")));
      DSHist->Add((MnvH1D*)(mcFile->Get(varName+"_InnerDSPlastic"+tag+material+"_background_1neutPi")));
      DSHist->Add((MnvH1D*)(mcFile->Get(varName+"_InnerDSPlastic"+tag+material+"_background_NPi")));
      DSHist->Add((MnvH1D*)(mcFile->Get(varName+"_InnerDSPlastic"+tag+material+"_background_Other")));
      DSHist_noTag = (MnvH1D*)(mcFile->Get(varName+"_InnerDSPlastic"+material+"_selected_signal_reco"))->Clone();
      DSHist_noTag->Add((MnvH1D*)(mcFile->Get(varName+"_InnerDSPlastic"+material+"_background_1chargePi")));
      DSHist_noTag->Add((MnvH1D*)(mcFile->Get(varName+"_InnerDSPlastic"+material+"_background_1neutPi")));
      DSHist_noTag->Add((MnvH1D*)(mcFile->Get(varName+"_InnerDSPlastic"+material+"_background_NPi")));
      DSHist_noTag->Add((MnvH1D*)(mcFile->Get(varName+"_InnerDSPlastic"+material+"_background_Other")));
    }
    DSHist->Scale(POTscale);
    DSHist_noTag->Scale(POTscale);

    MnvH1D* bkgNNeutPiHist = neutPiHist->Clone();
    MnvH1D* bkgNNeutPiHist_noTag = neutPiHist_noTag->Clone();
    bkgNNeutPiHist->Add(NPiHist);
    bkgNNeutPiHist_noTag->Add(NPiHist_noTag);

    MnvH1D* bkg1PiHist = neutPiHist->Clone();
    MnvH1D* bkg1PiHist_noTag = neutPiHist_noTag->Clone();
    bkg1PiHist->Add(chargePiHist);
    bkg1PiHist_noTag->Add(chargePiHist_noTag);

    //
    MnvH1D* QEHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_QE"))->Clone();
    MnvH1D* QEHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_bkg_IntType_QE"))->Clone();
    QEHist->Scale(POTscale);
    QEHist_noTag->Scale(POTscale);
    MnvH1D* RESHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_RES"))->Clone();
    MnvH1D* RESHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_bkg_IntType_RES"))->Clone();
    RESHist->Scale(POTscale);
    RESHist_noTag->Scale(POTscale);
    MnvH1D* DISHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_DIS"))->Clone();
    MnvH1D* DISHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_bkg_IntType_DIS"))->Clone();
    DISHist->Scale(POTscale);
    DISHist_noTag->Scale(POTscale);
    MnvH1D* MECHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_2p2h"))->Clone();
    MnvH1D* MECHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_bkg_IntType_2p2h"))->Clone();
    MECHist->Scale(POTscale);
    MECHist_noTag->Scale(POTscale);
    MnvH1D* OtherIntTypeHist = (MnvH1D*)(mcFile->Get(grabName+"_bkg_IntType_Other"))->Clone();
    MnvH1D* OtherIntTypeHist_noTag = (MnvH1D*)(mcFile->Get(grabName_noTag+"_bkg_IntType_Other"))->Clone();
    OtherIntTypeHist->Scale(POTscale);
    OtherIntTypeHist_noTag->Scale(POTscale);

    MnvH1D* bkgNonRESHist = DISHist->Clone();
    bkgNonRESHist->Add(QEHist);
    bkgNonRESHist->Add(MECHist);
    bkgNonRESHist->Add(OtherIntTypeHist);
    MnvH1D* bkgNonRESHist_noTag = DISHist_noTag->Clone();
    bkgNonRESHist_noTag->Add(QEHist_noTag);
    bkgNonRESHist_noTag->Add(MECHist_noTag);
    bkgNonRESHist_noTag->Add(OtherIntTypeHist_noTag);

    //Modifying so that the tuned backgrounds are what is grabbed...
    //MnvH1D* bkgTotHist_Int = bkgNonRESHist->Clone();
    //bkgTotHist_Int->Add(RESHist);
    //bkgTotHist->Add(wrongNuclHist);
    //MnvH1D* bkgTotHist_Int_noTag = bkgNonRESHist_noTag->Clone();
    //bkgTotHist_Int_noTag->Add(RESHist_noTag);
    //bkgTotHist_noTag->Add(wrongNuclHist_noTag);

    MnvH1D* bkgTotHist = bkg1PiHist->Clone();
    bkgTotHist->Add(NPiHist);
    bkgTotHist->Add(otherHist);
    MnvH1D* bkgTotHist_noTag = bkg1PiHist_noTag->Clone();
    bkgTotHist_noTag->Add(NPiHist_noTag);
    bkgTotHist_noTag->Add(otherHist_noTag);

    ////const map<TString, vector<MnvH1D*>> mcFitHistsINPUT, const vector<MnvH1D*> dataFitHists, vector<vector<MnvH1D*>> nonFitHists, const map<TString, vector<TString>> tagsToSave
    map<TString, vector<MnvH1D*>> mcFitHists;
    mcFitHists["Signal"].push_back((MnvH1D*)sigHist_noTag->Clone());
    mcFitHists["Signal"].push_back((MnvH1D*)sigHist->Clone());
    mcFitHists["BKG"].push_back((MnvH1D*)bkgTotHist_noTag->Clone());
    mcFitHists["BKG"].push_back((MnvH1D*)bkgTotHist->Clone());

    vector<MnvH1D*> dataHists;
    dataHists.push_back((MnvH1D*)dataHist_noTag->Clone());
    dataHists.push_back((MnvH1D*)dataHist->Clone());
    
    vector<vector<MnvH1D*>> nonFitHists;
    vector<MnvH1D*> tmp;
    nonFitHists.push_back(tmp);
    nonFitHists.at(nonFitHists.size()-1).push_back((MnvH1D*)USHist_noTag->Clone());
    nonFitHists.at(nonFitHists.size()-1).push_back((MnvH1D*)USHist->Clone());

    nonFitHists.push_back(tmp);
    nonFitHists.at(nonFitHists.size()-1).push_back((MnvH1D*)DSHist_noTag->Clone());
    nonFitHists.at(nonFitHists.size()-1).push_back((MnvH1D*)DSHist->Clone());
    
    nonFitHists.push_back(tmp);
    nonFitHists.at(nonFitHists.size()-1).push_back((MnvH1D*)wrongNuclHist_noTag->Clone());
    nonFitHists.at(nonFitHists.size()-1).push_back((MnvH1D*)wrongNuclHist->Clone());

    map<TString, vector<TString>> nameKeys;
    nameKeys["Signal"]={"sig","signal"};
    nameKeys["BKG"]={"1chargePi","1neutPi","NPi","Other"};
    
    cout << "Fitting" << endl;
    map<TString,map<TString,MnvH1D*>> result;
    ////const map<TString, vector<MnvH1D*>> mcFitHistsINPUT, const vector<MnvH1D*> dataFitHists, const vector<MnvH1D*> nonFitHists, const map<TString, vector<TString>> tagsToSave
    map<TString, MnvH1D*> TESTresult = PerformTwoFunctionSolution(mcFitHists, dataHists, nonFitHists, name, nameKeys);

    /*
    map<TString, MnvH1D*> tmpMap;
    vector<map<TString, MnvH1D*>> scaledHistsAndNamesTEST(dataHists.size(),map<TString,MnvH1D*>(tmpMap));
    vector<map<TString, MnvH1D*>> unfitHistsAndNamesTEST(dataHists.size(),map<TString,MnvH1D*>(tmpMap));
    for (auto hists:fitTEST1A){
      for (auto hist:hists.second){
	for (unsigned int iRegion=0; iRegion<hist.second.size(); ++iRegion){
	  if (hists.first=="NonFit"){
	    unfitHistsAndNamesTEST.at(iRegion)[hist.first]=(MnvH1D*)hist.second.at(iRegion)->Clone();
	  }
	  else{
	    scaledHistsAndNamesTEST.at(iRegion)[hist.first]=(MnvH1D*)hist.second.at(iRegion)->Clone();
	    scaledHistsAndNamesTEST.at(iRegion)[hist.first]->Multiply(scaledHistsAndNamesTEST.at(iRegion)[hist.first], TESTresult[hist.first]);
	  }
	}
      }
    }

    for (unsigned int iRegion=0; iRegion<dataHists.size(); ++iRegion){
      if (!PathExists((string)(outDir+name+"_fitTEST_"+to_string(breakInner)+"_Region_"+(TString)(to_string(iRegion))+"_postFit.pdf"))){
	////DrawFromMnvH1Ds(dataHists.at(iRegion),scaledHistsAndNamesTEST.at(iRegion),unfitHistsAndNamesTEST.at(iRegion),true,outDir+name+"_fitTEST_"+to_string(breakInner)+"_Region_"+(TString)(to_string(iRegion))+"_postFit");
	////DrawFromMnvH1Ds(dataHists.at(iRegion),scaledHistsAndNamesTEST.at(iRegion),unfitHistsAndNamesTEST.at(iRegion),false,outDir+name+"_fitTEST_"+to_string(breakInner)+"_Region_"+(TString)(to_string(iRegion))+"_postFit");
      }
    }
    */

    for (auto hist:TESTresult) hist.second->SetDirectory(outFile);    
    outFile->Write();
    for (auto hist:TESTresult) delete hist.second;
    TESTresult.clear();

    //Do good cleanup here in a minute... I now see why there were issues before actually haha. I didn't delete everything I needed to when I started fitting signal too almost definitely... will try fixing here and see if able to move into the "Overhaul" code as well successfully.
    delete dataHist;
    delete sigHist;
    delete chargePiHist;
    delete neutPiHist;
    delete NPiHist;
    delete otherHist;
    delete USHist;
    delete DSHist;
    delete wrongNuclHist;
    delete bkgNNeutPiHist;
    delete bkg1PiHist;
    delete QEHist;
    delete RESHist;
    delete DISHist;
    delete MECHist;
    delete OtherIntTypeHist;
    delete bkgNonRESHist;
    delete bkgTotHist;
    //for (auto hist:fitHists1A) delete hist.second;
    //for (auto hist:unfitHists1A) delete hist.second;
    //for (auto hist:scaledHists1A) delete hist.second;
    for (auto hist:varsToSave) delete hist.second;
    varsToSave.clear();
  }

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  mcFile->Close();
  dataFile->Close();
  outFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
