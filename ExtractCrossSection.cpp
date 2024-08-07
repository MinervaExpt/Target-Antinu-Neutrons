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

  if(!(argc == 12 || argc == 14))
  {
    std::cerr << "Expected 11 or 13 arguments, but I got " << argc-1 << ".\n"
              << "USAGE: ExtractCrossSection <unfolding iterations> <data.root> <mc.root> <stop at eff. corr.> <varName> ...\n";
    return 1;
  }

  const int nIterations = std::stoi(argv[1]);
  TString dataFileName = argv[2];
  bool isMC = false;
  if (dataFileName.Contains("runEventLoopData")){
    std::cout << "Data file is a data file" << std::endl;
  }
  else if (dataFileName.Contains("runEventLoopMC")){
    std::cout << "Data file is an MC file. Setting isMC to true" << std::endl;
    isMC = true;
  }
  else {
    std::cout << "Data file is not a proper input. Exiting now" << std::endl;
    return 304123;
  }

  auto dataFile = TFile::Open(argv[2], "READ");
  if(!dataFile)
  {
    std::cerr << "Failed to open data file " << argv[2] << ".\n";
    return 2;
  }

  auto mcFile = TFile::Open(argv[3], "READ");
  if(!mcFile)
  {
    std::cerr << "Failed to open MC file " << argv[3] << ".\n";
    return 3;
  }

  bool stopAtEffCorr = (bool)atoi(argv[4]);
  std::string varName = std::string(argv[5]);
  int tgtZ = atoi(argv[6]);
  int numFluxUniv = atoi(argv[7]);
  int nuPDG = -14; //hard-coded for my analyses
  const std::string project_dir = "targets_12345_jointNueIMD";//Copied from Anezka for target fluxes
  bool multPOT = (bool)atoi(argv[8]);
  std::string background_naming = std::string(argv[9]); 
  if (background_naming != "bkg_IntType") background_naming = "background";

  TFile* fluxFile = nullptr;
  std::string fluxVarName = "";

  auto& frw6A = PlotUtils::flux_reweighter("minervame6A", nuPDG, true, numFluxUniv);//playlist hard-coded to 6A for all anti-nu. Anezka says conclusion was flux consistent enough that this is fine.
  auto& frw5A = PlotUtils::flux_reweighter("minervame5A", nuPDG, true, numFluxUniv);//playlist hard-coded to 6A for all anti-nu. Anezka says conclusion was flux consistent enough that this is fine.

  double frac5A = atof(argv[10]);
  double unfoldingFactor = atof(argv[11]);

  if (argc == 14){
    fluxFile = TFile::Open(argv[12],"READ");
    fluxVarName =std::string(argv[13]);
  }

  std::vector<std::string> crossSectionPrefixes;
  for(auto key: *dataFile->GetListOfKeys())
  {
    const std::string keyName = key->GetName();
    const size_t endOfPrefix = keyName.find("_data");
    if(endOfPrefix != std::string::npos) crossSectionPrefixes.push_back(keyName.substr(0, endOfPrefix));
  }

  const double mcPOT = util::GetIngredient<TParameter<double>>(*mcFile, "POTUsed")->GetVal(),
               dataPOT = util::GetIngredient<TParameter<double>>(*dataFile, "POTUsed")->GetVal();

  crossSectionPrefixes.clear();
  crossSectionPrefixes.push_back(varName);

  for(const auto& prefix: crossSectionPrefixes)
  {
    try
    {
      auto folded = util::GetIngredient<PlotUtils::MnvH1D>(*dataFile, "data", prefix);
      Plot(*folded, "data", prefix, tgtZ);
      auto migration = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "migration", prefix);
      auto effNum = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_numerator", prefix);
      auto effDenom = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "efficiency_denominator", prefix);
      auto simEventRate = effDenom->Clone(); //Make a copy for later

      //Look for backgrounds with <prefix>_<analysis>_Background_<name>
      std::vector<PlotUtils::MnvH1D*> backgrounds;
      bool skipCH = true;//Testflag to compare doing the old subtraction with the new.
      for(auto key: *mcFile->GetListOfKeys())
      {
        if(std::string(key->GetName()).find(prefix + "_" + background_naming + "_") != std::string::npos)
        {
	  if((std::string(key->GetName()).find(background_naming + "_USPlastic") == std::string::npos && std::string(key->GetName()).find(background_naming + "_DSPlastic") == std::string::npos) || !skipCH){
	    backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, key->GetName()));
	  }
        }
	//This assumes the tuned backgrounds from the plastic are the same as those from the material.
	else if (isInnerPlastic(std::string(key->GetName()), background_naming) && skipCH){
	  backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, key->GetName()));
	}
      }

      //There are no error bands in the data, but I need somewhere to put error bands on the results I derive from it.
      folded->AddMissingErrorBandsAndFillWithCV(*migration);

      //Basing my unfolding procedure for a differential cross section on Alex's MINERvA 101 talk at https://minerva-docdb.fnal.gov/cgi-bin/private/RetrieveFile?docid=27438&filename=whatsACrossSection.pdf&version=1

      //TODO: Remove these debugging plots when done
      auto toSubtract = std::accumulate(std::next(backgrounds.begin()), backgrounds.end(), (*backgrounds.begin())->Clone(),
                                        [](auto sum, const auto hist)
                                        {
                                          sum->Add(hist);
                                          return sum;
                                        });
      Plot(*toSubtract, "BackgroundSum", prefix, tgtZ);

      auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                           [mcPOT, dataPOT](auto sum, const auto hist)
                                           {
                                             std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT/mcPOT << " from " << sum->GetName() << "\n";
                                             sum->Add(hist, -dataPOT/mcPOT);
                                             return sum;
                                           });
      Plot(*bkgSubtracted, "backgroundSubtracted", prefix, tgtZ);

      auto outFile = TFile::Open((prefix + "_crossSection.root").c_str(), "CREATE");
      if(!outFile)
      {
        std::cerr << "Could not create a file called " << prefix + "_crossSection.root" << ".  Does it already exist?\n";
        return 5;
      }

      bkgSubtracted->Write("backgroundSubtracted");

      //d'Aogstini unfolding
      auto unfolded = UnfoldHist(bkgSubtracted, migration, nIterations);
      if(!unfolded) throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
      unfolded->ModifyStatisticalUnc(unfoldingFactor,"unfoldingCov");
      std::cout << "Survived asjuting the unfolded histogram.\n" << std::flush;
      Plot(*unfolded, "unfolded", prefix, tgtZ);
      unfolded->Clone()->Write("unfolded"); //TODO: Seg fault first appears when I uncomment this line
      std::cout << "Survived writing the unfolded histogram.\n" << std::flush; //This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().

      effNum->Divide(effNum, effDenom); //Only the 2 parameter version of MnvH1D::Divide()
                                        //handles systematics correctly.
      Plot(*effNum, "efficiency", prefix, tgtZ);
      effNum->Clone()->Write("efficiency");

      unfolded->Divide(unfolded, effNum);
      Plot(*unfolded, "efficiencyCorrected", prefix, tgtZ);

      unfolded->Clone()->Write("efficiencyCorrected");

      if (!stopAtEffCorr){
	PlotUtils::MnvH1D* flux;//Hard-coded for 1D. Might break... but currently necessary...

	std::string material;
	if (tgtZ == 6) material = "carbon";
	else if (tgtZ == 26) material = "iron";
	else if (tgtZ == 82) material = "lead";
	else if (tgtZ == 8){
	  material = ""; //Not sure this is correct... Going to ignore for now... cause water issues... also I'm not doing the nucleon division correctly for water, so really ignoring it for now...
	  if (argc != 14){
	    std::cout << "Returning random value for code since water requires a tracker flux to be provided." << std::endl;
	    return -12039;
	  }
	}
	else material = "badTarget";

	if (argc == 14) flux = util::GetIngredient<PlotUtils::MnvH1D>(*fluxFile,"reweightedflux_integrated",fluxVarName);
	else if (tgtZ == -1) flux = util::GetIngredient<PlotUtils::MnvH1D>(*mcFile, "reweightedflux_integrated", prefix);
	else{ 
	  flux = frw6A.GetIntegratedTargetFlux(nuPDG, material, unfolded, 0, 100, project_dir);
	  PlotUtils::MnvH1D* flux5A = frw5A.GetIntegratedTargetFlux(nuPDG, material, unfolded, 0, 100, project_dir);
	  flux->Scale((1.0-frac5A));
	  flux5A->Scale(frac5A);
	  flux->Add(flux5A);
	  delete flux5A;
	}

	/*
	const auto fiducialFound = std::find_if(mcFile->GetListOfKeys()->begin(), mcFile->GetListOfKeys()->end(),
						[&prefix](const auto key)
						{
						  const std::string keyName = key->GetName();
						  const size_t fiducialEnd = keyName.find("_fiducial_nucleons");
						  return (fiducialEnd != std::string::npos) && (prefix.find(keyName.substr(0, fiducialEnd)) != std::string::npos);
						});
	if(fiducialFound == mcFile->GetListOfKeys()->end()) throw std::runtime_error("Failed to find a number of nucleons that matches prefix " + prefix);

	auto nNucleons = util::GetIngredient<TParameter<double>>(*mcFile, (*fiducialFound)->GetName()); //Dan: Use the same truth fiducial volume for all extractions.  The acceptance correction corrects data back to this fiducial even if the reco fiducial cut is different.
	double nNuke = nNucleons->GetVal();
	*/
	double nNuke=1.0;
	double nNukeMC = 1.0;
	if (tgtZ != -1){
	  nNukeMC = GetTotalScatteringCenters(tgtZ,true);
	  nNuke = GetTotalScatteringCenters(tgtZ, isMC);
	}
	else{
	  nNukeMC = GetTotalScatteringCenters(99,true);
	  nNuke = GetTotalScatteringCenters(99, isMC);
	}

	std::cout << "No. of nucleons: " << nNuke << std::endl;
	std::cout << "No. of antineutrinos: " << flux->GetBinContent(1) << std::endl;
	std::cout << "No. of antineutrinos multiplied: " << flux->GetBinContent(1)*dataPOT << std::endl;
	auto crossSection = normalize(unfolded, flux, nNuke, dataPOT);
	if (multPOT) crossSection->Scale(dataPOT);
	if (!isMC){
	  auto MassSyst = GetTargetMassSystHist(crossSection, tgtZ);
	  crossSection->AddMissingErrorBandsAndFillWithCV(*MassSyst);
	  crossSection->Multiply(crossSection,MassSyst);
	}
	Plot(*crossSection, "crossSection", prefix, tgtZ);
	outFile->cd();
	crossSection->Clone()->Write("crossSection");
      
	//Write a "simulated cross section" to compare to the data I just extracted.
	//If this analysis passed its closure test, this should be the same cross section as
	//what GENIEXSecExtract would produce.
	normalize(simEventRate, flux, nNukeMC, mcPOT);
	if (multPOT) simEventRate->Scale(dataPOT);  

	Plot(*simEventRate, "simulatedCrossSection", prefix, tgtZ);
	simEventRate->Write("simulatedCrossSection");
      }
    }
    catch(const std::runtime_error& e)
    {
      std::cerr << "Failed to extra a cross section for prefix " << prefix << ": " << e.what() << "\n";
      return 4;
    }
  }

  return 0;
}
