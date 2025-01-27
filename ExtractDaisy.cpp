//File: ExtractCrossSection.cpp
//Brief: Given data and MC files from analyses/studies/CrossSection.h, extract a 1D differential cross section.
//       Subtracts backgrounds, performs unfolding, applies efficiency x acceptance correction, and 
//       divides by flux and number of nucleons.  Writes a .root file with the cross section histogram.
//
//Usage: ExtractCrossSection <unfolding iterations> <data.root> <mc.root> <stop at efficiency correction> <varName> <tgtZ> <no. Flux Universes> <multiply by data POT> <background naming> <5A Fraction> : optional <flux_file> <fluxVarName>
//
//Author: Andrew Olivier aolivier@ur.rochester.edu

//util includes
#include "util/GetIngredient.h"

//UnfoldUtils includes
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#include "MinervaUnfold/MnvUnfold.h"

//PlotUtils includes
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/TargetUtils.h"
#include "PlotUtils/FluxReweighter.h"
#include "PlotUtils/TargetMassSystematics.h"
#pragma GCC diagnostic pop

//ROOT includes
#include "TH1D.h"
#include "TFile.h"
#include "TKey.h"
#include "TParameter.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TColor.h"

//Cintex is only needed for older ROOT versions like the GPVMs.
////Let CMake decide whether it's needed.
#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

//c++ includes
#include <iostream>
#include <exception>
#include <algorithm>
#include <numeric>

//Convince the STL to talk to TIter so I can use std::find_if()
namespace std
{
  template <>
  struct iterator_traits<TIter>
  {
    using value_type = TObject;
    using pointer = TObject*;
    using reference = TObject&;
    using iterator_category = forward_iterator_tag;
  };
}

std::map<int,TString> matTag={{-1, "#it{Tracker (CH)}"},{6, "#it{Carbon (C)}"},{8, "#it{Water (H_{2}O)}"},{26, "#it{Iron (Fe)}"},{82,"#it{Lead (Pb)}"}};

PlotUtils::MnvH1D* GetTargetMassSystHist(PlotUtils::MnvH1D* hTemp, int tgtZ){
  PlotUtils::MnvH1D* hOut = nullptr;
  if (tgtZ == -1) hOut = PlotUtils::GetNTargetsScintillatorHist<PlotUtils::MnvH1D>(1.0,hTemp);
  else if (tgtZ == 6) hOut = PlotUtils::GetNTargetsCarbonHist<PlotUtils::MnvH1D>(1.0,hTemp);
  else if (tgtZ == 26) hOut = PlotUtils::GetNTargetsIronHist<PlotUtils::MnvH1D>(1.0,hTemp);
  else if (tgtZ == 82) hOut = PlotUtils::GetNTargetsLeadHist<PlotUtils::MnvH1D>(1.0,hTemp);
  else if (tgtZ == 8) hOut = PlotUtils::GetNTargetsWaterHist<PlotUtils::MnvH1D>(1.0,hTemp);
  return hOut;
}

//Check if the object is meant to be subtracted as the background inner plastic.
//This is forcing the variable to be pTmu... could change but don't need to yet.
bool isInnerPlastic(std::string name, std::string bkgName){
  bool isGood = false;
  bool isInnerPlastic = (name.find("pTmu_InnerUSPlastic") != std::string::npos || name.find("pTmu_InnerDSPlastic") != std::string::npos);
  bool isNotSB = (name.find("_PreRecoilCut_") == std::string::npos);
  bool isSignal = (name.find("_selected_signal_reco") != std::string::npos);
  bool isCorrectBKG = (name.find(bkgName+"_") != std::string::npos);
  isGood = (isInnerPlastic && isNotSB && (isSignal || isCorrectBKG));
  return isGood;
}

//Plot a step in cross section extraction.
void Plot(PlotUtils::MnvH1D& hist, const std::string& stepName, const std::string& prefix, const int tgtZ)
{
  TCanvas can(stepName.c_str());
  hist.GetCVHistoWithError().Clone()->Draw();
  can.Print((prefix + "_" + stepName + ".png").c_str());
  can.Print((prefix + "_" + stepName + ".pdf").c_str());
  can.Print((prefix + "_" + stepName + ".C").c_str());

  //Uncertainty summary
  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kCCQEAntiNuStyle);
  plotter.axis_maximum = 0.4;
  TLatex* tex = new TLatex(0.7375,0.6,matTag[tgtZ]);
  tex->SetNDC();
  tex->SetTextColor(TColor::GetColor("#ff0000"));
  tex->SetTextFont(43);
  tex->SetTextSize(40);
  tex->SetLineWidth(3);
  tex->Draw();

  plotter.DrawErrorSummary(&hist);
  can.Print((prefix + "_" + stepName + "_uncertaintySummary.png").c_str());
  can.Print((prefix + "_" + stepName + "_uncertaintySummary.pdf").c_str());
  can.Print((prefix + "_" + stepName + "_uncertaintySummary.C").c_str());

  /*
  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Other");
  can.Print((prefix + "_" + stepName + "_otherUncertainties.png").c_str());
  can.Print((prefix + "_" + stepName + "_otherUncertainties.pdf").c_str());
  can.Print((prefix + "_" + stepName + "_otherUncertainties.C").c_str());
  */

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Neutron Detection");
  can.Print((prefix + "_" + stepName + "_neutronUncertainties.png").c_str());
  can.Print((prefix + "_" + stepName + "_neutronUncertainties.pdf").c_str());
  can.Print((prefix + "_" + stepName + "_neutronUncertainties.C").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Cross Section Models");
  can.Print((prefix + "_" + stepName + "_xSecUncertainties.png").c_str());
  can.Print((prefix + "_" + stepName + "_xSecUncertainties.pdf").c_str());
  can.Print((prefix + "_" + stepName + "_xSecUncertainties.C").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Calorimetric Response");
  can.Print((prefix + "_" + stepName + "_RecoilUncertainties.png").c_str());
  can.Print((prefix + "_" + stepName + "_RecoilUncertainties.pdf").c_str());
  can.Print((prefix + "_" + stepName + "_RecoilUncertainties.C").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "FSI Models");
  can.Print((prefix + "_" + stepName + "_FSIUncertainties.png").c_str());
  can.Print((prefix + "_" + stepName + "_FSIUncertainties.pdf").c_str());
  can.Print((prefix + "_" + stepName + "_FSIUncertainties.C").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Muon Reconstruction");
  can.Print((prefix + "_" + stepName + "_MuonUncertainties.png").c_str());
  can.Print((prefix + "_" + stepName + "_MuonUncertainties.pdf").c_str());
  can.Print((prefix + "_" + stepName + "_MuonUncertainties.C").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "GEANT4");
  can.Print((prefix + "_" + stepName + "_GEANTUncertainties.png").c_str());
  can.Print((prefix + "_" + stepName + "_GEANTUncertainties.pdf").c_str());
  can.Print((prefix + "_" + stepName + "_GEANTUncertainties.C").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Normalization");
  can.Print((prefix + "_" + stepName + "_NormalizationUncertainties.png").c_str());
  can.Print((prefix + "_" + stepName + "_NormalizationUncertainties.pdf").c_str());
  can.Print((prefix + "_" + stepName + "_NormalizationUncertainties.C").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "MnvTune V1");
  can.Print((prefix + "_" + stepName + "_MnvTuneUncertainties.png").c_str());
  can.Print((prefix + "_" + stepName + "_MnvTuneUncertainties.pdf").c_str());
  can.Print((prefix + "_" + stepName + "_MnvTuneUncertainties.C").c_str());
}

//Unfolding function from Aaron Bercelle
//TODO: Trim it down a little?  Remove that static?
PlotUtils::MnvH1D* UnfoldHist( PlotUtils::MnvH1D* h_folded, PlotUtils::MnvH2D* h_migration, int num_iter )
{
  static MinervaUnfold::MnvUnfold unfold;
  PlotUtils::MnvH1D* h_unfolded = nullptr;

  if (num_iter==0) return h_folded;
  
  //bool bUnfolded = false;

  TMatrixD dummyCovMatrix;
  if(!unfold.UnfoldHisto( h_unfolded, dummyCovMatrix, h_migration, h_folded, RooUnfold::kBayes, num_iter, true, false ))
    return nullptr;

  /////////////////////////////////////////////////////////////////////////////////////////  
  //No idea if this is still needed
  //Probably.  This gets your stat unfolding covariance matrix
  TMatrixD unfoldingCovMatrixOrig; 
  int correctNbins;
  int matrixRows;  
  TH1D* hUnfoldedDummy  = new TH1D(h_unfolded->GetCVHistoWithStatError());
  TH1D* hRecoDummy      = new TH1D(h_migration->ProjectionX()->GetCVHistoWithStatError());
  TH1D* hTruthDummy     = new TH1D(h_migration->ProjectionY()->GetCVHistoWithStatError());
  TH1D* hBGSubDataDummy = new TH1D(h_folded->GetCVHistoWithStatError());
  TH2D* hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());
  unfold.UnfoldHisto(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy,RooUnfold::kBayes, num_iter);//Stupid RooUnfold.  This is dummy, we don't need iterations

  correctNbins=hUnfoldedDummy->fN;
  matrixRows=unfoldingCovMatrixOrig.GetNrows();
  if(correctNbins!=matrixRows){
    std::cout << "****************************************************************************" << std::endl;
    std::cout << "*  Fixing unfolding matrix size because of RooUnfold bug. From " << matrixRows << " to " << correctNbins << std::endl;
    std::cout << "****************************************************************************" << std::endl;
    // It looks like this, since the extra last two bins don't have any content
    unfoldingCovMatrixOrig.ResizeTo(correctNbins, correctNbins);
  }
  
  for(int i=0; i<unfoldingCovMatrixOrig.GetNrows(); ++i) unfoldingCovMatrixOrig(i,i)=0;
  delete hUnfoldedDummy;
  delete hMigrationDummy;
  delete hRecoDummy;
  delete hTruthDummy;
  delete hBGSubDataDummy;
  
  h_unfolded->PushCovMatrix("unfoldingCov",unfoldingCovMatrixOrig);

  /////////////////////////////////////////////////////////////////////////////////////////  
  return h_unfolded;
}

double GetTotalScatteringCenters(int targetZ, bool isMC)
{
  // TARGET INFO
  PlotUtils::TargetUtils targetInfo;
  double Nucleons = 0.0;

  // Target 1 is generally excluded due to rock muon contamination (in the inclusive analysis), keeping for now...
  if(targetZ == 6){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ); // Target 3
  }
  else if(targetZ == 26){                                                                                                                                                                                              
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 1, targetZ, isMC ) // Target 1
      + targetInfo.GetPassiveTargetNNucleons( 2, targetZ, isMC ) // Target 2                                                                                                                                
      + targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ) // Target 3                                                                                                                                
      + targetInfo.GetPassiveTargetNNucleons( 5, targetZ, isMC );// Target 5
  }
  else if(targetZ == 82){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 1, targetZ, isMC ) // Target 2
      + targetInfo.GetPassiveTargetNNucleons( 2, targetZ, isMC ) // Target 2
      + targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ) // Target 3
      + targetInfo.GetPassiveTargetNNucleons( 4, targetZ, isMC ) // Target 4
      + targetInfo.GetPassiveTargetNNucleons( 5, targetZ, isMC );// Target 5
  }
  else if(targetZ == 8){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 6, targetZ, isMC );//Water
      //+ targetInfo.GetPassiveTargetNNucleons( 6, 1, isMC );//Water Hydrogen is handled above. This was wrong from before. Explains why it seemed about a factor of 2 low... essentially divided by the number of water nucleons twice...
      }
  else if(targetZ > 90 ){
    Nucleons = targetInfo.GetTrackerNNucleons(5980, 8422, isMC, 850);
  }
  return Nucleons;
}                                                                                                                                                                                                                 

double GetProtonScatteringCenters(int targetZ, bool isMC)
{
  // TARGET INFO
  PlotUtils::TargetUtils targetInfo;
  double Nucleons = 0.0;

  // Target 1 is generally excluded due to rock muon contamination (in the inclusive analysis), keeping for now...
  if(targetZ == 6){
    Nucleons = targetInfo.GetPassiveTargetNProtons( 3, targetZ, isMC ); // Target 3
  }
  else if(targetZ == 26){
    Nucleons = targetInfo.GetPassiveTargetNProtons( 1, targetZ, isMC ) // Target 1
      + targetInfo.GetPassiveTargetNProtons( 2, targetZ, isMC ) // Target 2
      + targetInfo.GetPassiveTargetNProtons( 3, targetZ, isMC ) // Target 3
      + targetInfo.GetPassiveTargetNProtons( 5, targetZ, isMC );// Target 5
  }
  else if(targetZ == 82){
    Nucleons = targetInfo.GetPassiveTargetNProtons( 1, targetZ, isMC ) // Target 1
      + targetInfo.GetPassiveTargetNProtons( 2, targetZ, isMC ) // Target 2
      + targetInfo.GetPassiveTargetNProtons( 3, targetZ, isMC ) // Target 3
      + targetInfo.GetPassiveTargetNProtons( 4, targetZ, isMC ) // Target 4
      + targetInfo.GetPassiveTargetNProtons( 5, targetZ, isMC );// Target 5                                                                                                                                       
  }
  else if(targetZ == 8){
    Nucleons = targetInfo.GetPassiveTargetNProtons( 6, targetZ, isMC );//Water
    //+ targetInfo.GetPassiveTargetNNucleons( 6, 1, isMC );//Water Hydrogen is handled above. This was wrong from before. Explains why it seemed about a factor of 2 low... essentially divided by the number of water nucleons twice...                                                                                                                                                                                           
  }
  else if(targetZ > 90 ){
    Nucleons = targetInfo.GetTrackerNProtons(5980, 8422, isMC, 850);
  }
  return Nucleons;
}

//The final step of cross section extraction: normalize by flux, bin width, POT, and number of targets
PlotUtils::MnvH1D* normalize(PlotUtils::MnvH1D* efficiencyCorrected, PlotUtils::MnvH1D* fluxIntegral, const double nNucleons, const double POT)
{
  std::cout << "Dividing" << std::endl;
  efficiencyCorrected->Divide(efficiencyCorrected, fluxIntegral);
  
  std::cout << "Scaling" << std::endl;
  efficiencyCorrected->Scale(1./nNucleons/POT);
  std::cout << "Units" << std::endl;
  efficiencyCorrected->Scale(1.e4); //Flux histogram is in m^-2, but convention is to report cm^2
  std::cout << "Bin Width" << std::endl;
  efficiencyCorrected->Scale(1., "width");

  std::cout << "Returning" << std::endl;
  return efficiencyCorrected;
}

int main(const int argc, const char** argv)
{
  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable(); //Needed to look up dictionaries for PlotUtils classes like MnvH1D
  #endif

  TH1::AddDirectory(kFALSE); //Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).

  if((argc < 18))
  {
    std::cerr << "Expected at least 12 daisy files and many other arguments \n";
    std::cerr << "Got " << argc << "\n";
    return 1;
  }

  TString POTFileName = argv[1];
  TFile* POTFile = TFile::Open(POTFileName);

  double POT = ((TParameter<double>*)(POTFile->Get("POTUsed")))->GetVal();
  
  std::string outFileName = argv[2];
  TFile* outFile = new TFile((TString)(outFileName)+".root","RECREATE");
  
  int tgtZ = atoi(argv[3]);
  int numFluxUniv = atoi(argv[4]);
  int nuPDG = -14; //hard-coded for my analyses
  const std::string project_dir = "targets_12345_jointNueIMD";//Copied from Anezka for target fluxes

  auto& frw6A = PlotUtils::flux_reweighter("minervame6A", nuPDG, true, numFluxUniv);//playlist hard-coded to 6A for all anti-nu. Anezka says conclusion was flux consistent enough that this is fine.
  auto& frw5A = PlotUtils::flux_reweighter("minervame5A", nuPDG, true, numFluxUniv);//playlist hard-coded to 6A for all anti-nu. Anezka says conclusion was flux consistent enough that this is fine.

  double frac5A = atof(argv[5]);

  PlotUtils::MnvH1D* flux;//Hard-coded for 1D. Might break... but currently necessary...
  
  std::string material;
  if (tgtZ == -1) material = "tracker";
  else if (tgtZ == 6) material = "carbon";
  else if (tgtZ == 26) material = "iron";
  else if (tgtZ == 82) material = "lead";
  else {
    std::cout << "This ain't a material you can do this with." << std::endl;
    return 123456;
  }  

  std::map<int,PlotUtils::MnvH1D*> simulatedEventsDaisy;
  std::map<int,PlotUtils::MnvH1D*> dataEffCorrDaisy;

  std::vector<TString> petalFiles;

  for (int i=6; i<argc; ++i){
    petalFiles.push_back(argv[i]);
  }

  int i=0;
  for (auto fileName: petalFiles){
    std::cout << fileName << std::endl;
    TFile* file = TFile::Open(fileName);
    simulatedEventsDaisy[i]=((PlotUtils::MnvH1D*)file->Get("simulatedEventRate")->Clone());
    dataEffCorrDaisy[i]=((PlotUtils::MnvH1D*)file->Get("efficiencyCorrected")->Clone());
    file->Close();
    delete file;
    ++i;
  }
  
  PlotUtils::MnvH1D* unfolded=frw6A.GetReweightedDaisySum(nuPDG, material, dataEffCorrDaisy, project_dir);
  PlotUtils::MnvH1D* simEventRate=frw6A.GetReweightedDaisySum(nuPDG, material, simulatedEventsDaisy, project_dir);
  
  flux = frw6A.GetIntegratedTargetFlux(nuPDG, material, unfolded, 0, 100, project_dir);
  PlotUtils::MnvH1D* flux5A = frw5A.GetIntegratedTargetFlux(nuPDG, material, unfolded, 0, 100, project_dir);
  flux->Scale((1.0-frac5A));
  flux5A->Scale(frac5A);
  flux->Add(flux5A);
  delete flux5A;

  double nNuke = 1.0;
  double nNukeMC = 1.0;
  nNukeMC = GetTotalScatteringCenters(99, true);
  nNuke = GetTotalScatteringCenters(99, false);
  
  double nProt = 1.0;
  double nProtMC = 1.0;
  nProtMC = GetProtonScatteringCenters(99, true);
  nProt = GetProtonScatteringCenters(99, false);
  
  std::cout << "No. of nucleons: " << nNuke << std::endl;
  std::cout << "No. of protons: " << nProt << std::endl;
  std::cout << "No. of antineutrinos: " << flux->GetBinContent(1) << std::endl;
  std::cout << "No. of antineutrinos multiplied: " << flux->GetBinContent(1)*POT << std::endl;

  auto unfoldedProt = unfolded->Clone();
  auto simEventRateProt = simEventRate->Clone();
  
  auto crossSection = normalize(unfolded, flux, nNuke, POT);
  auto MassSyst = GetTargetMassSystHist(crossSection, -1);
  crossSection->AddMissingErrorBandsAndFillWithCV(*MassSyst);
  crossSection->Multiply(crossSection,MassSyst);
  Plot(*crossSection, "crossSection", outFileName, tgtZ);
  outFile->cd();
  crossSection->Clone()->Write("crossSection");
  
  auto crossSectionProt = normalize(unfoldedProt, flux, nProt, POT);
  crossSectionProt->AddMissingErrorBandsAndFillWithCV(*MassSyst);
  crossSectionProt->Multiply(crossSectionProt,MassSyst);
  Plot(*crossSectionProt, "crossSectionProt", outFileName, tgtZ);
  outFile->cd();
  crossSectionProt->Clone()->Write("crossSectionProt");
  
  normalize(simEventRate, flux, nNukeMC, POT);

  Plot(*simEventRate, "simulatedCrossSection", outFileName, tgtZ);
  simEventRate->Write("simulatedCrossSection");

  normalize(simEventRateProt, flux, nProtMC, POT);

  Plot(*simEventRateProt, "simulatedCrossSectionProt", outFileName, tgtZ);
  simEventRateProt->Write("simulatedCrossSectionProt");
  outFile->Close();

  return 0;
}
