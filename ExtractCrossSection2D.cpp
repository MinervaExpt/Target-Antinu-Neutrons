//File: ExtractCrossSection2D.cpp
//Brief: Given data and MC files from analyses/studies/CrossSection.h, extract a 1D differential cross section.
//       Subtracts backgrounds, performs unfolding, applies efficiency x acceptance correction, and 
//       divides by flux and number of nucleons.  Writes a .root file with the cross section histogram.
//
//Usage: ExtractCrossSection2D <unfolding iterations> <data.root> <mc.root> <stop at efficiency correction> <varName> <tgtZ> : optional <flux_file> <fluxVarName>
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
#pragma GCC diagnostic pop

//ROOT includes
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TKey.h"
#include "TParameter.h"
#include "TCanvas.h"

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

//Plot a step in cross section extraction.
void Plot(PlotUtils::MnvH1D& hist, const std::string& stepName, const std::string& prefix)
{
  TCanvas can(stepName.c_str());
  hist.GetCVHistoWithError().Clone()->Draw();
  can.Print((prefix + "_" + stepName + ".png").c_str());

  //Uncertainty summary
  PlotUtils::MnvPlotter plotter;
  plotter.ApplyStyle(PlotUtils::kCCQEAntiNuStyle);
  plotter.axis_maximum = 0.4;

  plotter.DrawErrorSummary(&hist);
  can.Print((prefix + "_" + stepName + "_uncertaintySummary.png").c_str());

  plotter.DrawErrorSummary(&hist, "TR", true, true, 1e-5, false, "Other");
  can.Print((prefix + "_" + stepName + "_otherUncertainties.png").c_str());
}

//Unfolding function from Aaron Bercelle
PlotUtils::MnvH2D* UnfoldHist( PlotUtils::MnvH2D* h_folded, PlotUtils::MnvH2D* h_selSig, PlotUtils::MnvH2D* h_effNum, PlotUtils::MnvH2D* h_migration, int num_iter )
{
  static MinervaUnfold::MnvUnfold unfold;
  PlotUtils::MnvH2D* h_unfolded = nullptr;

  //bool bUnfolded = false;

  //TMatrixD dummyCovMatrix;
  if(!unfold.UnfoldHisto2D( h_unfolded, h_migration, h_selSig, h_effNum, h_folded, num_iter, true, false ))
    return nullptr;

  /////////////////////////////////////////////////////////////////////////////////////////  
  //No idea if this is correct in 2D... Think it's handled by the above MnvUnfold::UnfoldHisto2D()... code looks the same...
  //Probably.  This gets your stat unfolding covariance matrix
  /*
  TMatrixD unfoldingCovMatrixOrig; 
  int correctNbins;
  int matrixRows;  
  TH2D* hUnfoldedDummy  = new TH2D(h_unfolded->GetCVHistoWithStatError());
  TH2D* hRecoDummy      = new TH2D(h_selSig->GetCVHistoWithStatError());
  TH2D* hTruthDummy     = new TH2D(h_effNum->GetCVHistoWithStatError());
  TH2D* hBGSubDataDummy = new TH2D(h_folded->GetCVHistoWithStatError());
  TH2D* hMigrationDummy = new TH2D(h_migration->GetCVHistoWithStatError());
  unfold.UnfoldHisto2D(hUnfoldedDummy, unfoldingCovMatrixOrig, hMigrationDummy, hRecoDummy, hTruthDummy, hBGSubDataDummy, num_iter);//Stupid RooUnfold.  This is dummy, we don't need iterations

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
  */
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
  if(targetZ == 26){                                                                                                                                                                                              
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 1, targetZ, isMC ) // Target 1
      + targetInfo.GetPassiveTargetNNucleons( 2, targetZ, isMC ) // Target 2                                                                                                                                
      + targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ) // Target 3                                                                                                                                
      + targetInfo.GetPassiveTargetNNucleons( 5, targetZ, isMC );// Target 5
  }
  if(targetZ == 82){
    Nucleons = targetInfo.GetPassiveTargetNNucleons( 1, targetZ, isMC ) // Target 2
      + targetInfo.GetPassiveTargetNNucleons( 2, targetZ, isMC ) // Target 2
      + targetInfo.GetPassiveTargetNNucleons( 3, targetZ, isMC ) // Target 3
      + targetInfo.GetPassiveTargetNNucleons( 4, targetZ, isMC ) // Target 4
      + targetInfo.GetPassiveTargetNNucleons( 5, targetZ, isMC );// Target 5
  }
  if(targetZ > 90 ){
    Nucleons = targetInfo.GetTrackerNNucleons(5980, 8422, isMC, 850);
  }
  return Nucleons;
}                                                                                                                                                                                                                 

PlotUtils::MnvH2D* normalize(PlotUtils::MnvH2D* efficiencyCorrected, PlotUtils::MnvH2D* fluxIntegral, const double nNucleons, const double POT)
{
  efficiencyCorrected->Divide(efficiencyCorrected, fluxIntegral);

  efficiencyCorrected->Scale(1./nNucleons/POT);
  efficiencyCorrected->Scale(1.e4); //Flux histogram is in m^-2, but convention is to report cm^2
  efficiencyCorrected->Scale(1., "width");

  return efficiencyCorrected;
}

int main(const int argc, const char** argv)
{
  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable(); //Needed to look up dictionaries for PlotUtils classes like MnvH1D
  #endif

  TH1::AddDirectory(kFALSE); //Needed so that MnvH1D gets to clean up its own MnvLatErrorBands (which are TH1Ds).

  if(!(argc == 7 || argc == 9))
  {
    std::cerr << "Expected 6 or 8 arguments, but I got " << argc-1 << ".\n"
              << "USAGE: ExtractCrossSection <unfolding iterations> <data.root> <mc.root> <stop at eff. corr.> <varName> \n";
    return 1;
  }

  const int nIterations = std::stoi(argv[1]);
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
  std::string dirName = "TwoD";//Add argument for this later...
  std::string varName = std::string(argv[5]);
  std::cout << "At set: " << dirName+"/"+varName << std::endl;
  int tgtZ = atoi(argv[6]);

  TFile* fluxFile = nullptr;
  std::string fluxVarName = "";
  if (argc == 9){
    fluxFile = TFile::Open(argv[7],"READ");
    fluxVarName =std::string(argv[8]);
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
  crossSectionPrefixes.push_back(dirName+"/"+varName);

  for(const auto& prefix: crossSectionPrefixes)
  {
    try
    {
      auto folded = util::GetIngredient<PlotUtils::MnvH2D>(*dataFile, "data", prefix);
      //Plot(*folded, "data", prefix);
      auto selSig = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "selected_signal_reco", prefix);
      auto migration = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "migration", prefix);
      auto effNum = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "efficiency_numerator", prefix);
      auto effDenom = util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "efficiency_denominator", prefix);
      auto simEventRate = effDenom->Clone(); //Make a copy for later

      std::cout << "Check folded: " << folded->Integral() << std::endl;
      std::cout << "Check selSig: " << selSig->Integral() << std::endl;
      std::cout << "Check migration: " << migration->Integral() << std::endl;
      std::cout << "Check effNum: " << effNum->Integral() << std::endl;
      std::cout << "Check effDenom: " << effDenom->Integral() << std::endl;
      std::cout << "Check simEventRate: " << simEventRate->Integral() << std::endl;

      //Look for backgrounds with <prefix>_<analysis>_Background_<name>
      std::vector<PlotUtils::MnvH2D*> backgrounds;
      TDirectoryFile* dir = (TDirectoryFile*)mcFile->Get(dirName.c_str());
      for(auto key: *dir->GetListOfKeys())
      {
	std::cout << "Object: " << key->GetName() << std::endl;
        if(std::string(key->GetName()).find(varName + "_background_") != std::string::npos)
        {
	  std::cout << "Getting background: " << dirName+"/"+key->GetName() << std::endl;
          backgrounds.push_back(util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, dirName+"/"+key->GetName()));
        }
      }

      std::cout << "BKG size: " << backgrounds.size() << std::endl;

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

      //Plot(*toSubtract, "BackgroundSum", prefix);

      auto bkgSubtracted = std::accumulate(backgrounds.begin(), backgrounds.end(), folded->Clone(),
                                           [mcPOT, dataPOT](auto sum, const auto hist)
                                           {
                                             std::cout << "Subtracting " << hist->GetName() << " scaled by " << -dataPOT/mcPOT << " from " << sum->GetName() << "\n";
                                             sum->Add(hist, -dataPOT/mcPOT);
                                             return sum;
                                           });
      //Plot(*bkgSubtracted, "backgroundSubtracted", prefix);

      auto outFile = TFile::Open((varName + "_crossSection.root").c_str(), "CREATE");
      if(!outFile)
      {
        std::cerr << "Could not create a file called " << varName + "_crossSection.root" << ".  Does it already exist?\n";
        return 5;
      }

      bkgSubtracted->Write("backgroundSubtracted");

      //d'Aogstini unfolding
      auto unfolded = UnfoldHist(bkgSubtracted, selSig, effNum, migration, nIterations);
      if(!unfolded) throw std::runtime_error(std::string("Failed to unfold ") + folded->GetName() + " using " + migration->GetName());
      //Plot(*unfolded, "unfolded", prefix);
      unfolded->Clone()->Write("unfolded"); //TODO: Seg fault first appears when I uncomment this line
      std::cout << "Survived writing the unfolded histogram.\n" << std::flush; //This is evidence that the problem is on the final file Write() and not unfolded->Clone()->Write().

      std::cout << "Eff" << std::endl;
      effNum->Divide(effNum, effDenom); //Only the 2 parameter version of MnvH1D::Divide()
                                        //handles systematics correctly.
      //Plot(*effNum, "efficiency", prefix);

      std::cout << "Unfold/Eff" << std::endl;
      unfolded->Divide(unfolded, effNum);
      //Plot(*unfolded, "efficiencyCorrected", prefix);

      std::cout << "Unfold Clone" << std::endl;
      unfolded->Clone()->Write("efficiencyCorrected");

      std::cout << "Should be roughly the end..." << std::endl;

      if (!stopAtEffCorr){
	auto flux = (argc==9) ? util::GetIngredient<PlotUtils::MnvH2D>(*fluxFile,"reweightedflux_integrated",fluxVarName) : util::GetIngredient<PlotUtils::MnvH2D>(*mcFile, "reweightedflux_integrated", prefix);
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
	if (tgtZ != -1){
	  nNuke = GetTotalScatteringCenters(tgtZ,true);
	}
	else{
	  nNuke = GetTotalScatteringCenters(99,true);
	}

	auto crossSection = normalize(unfolded, flux, nNuke, dataPOT);
	//Plot(*crossSection, "crossSection", prefix);
	crossSection->Clone()->Write("crossSection");
      
	//Write a "simulated cross section" to compare to the data I just extracted.
	//If this analysis passed its closure test, this should be the same cross section as
	//what GENIEXSecExtract would produce.
	normalize(simEventRate, flux, nNuke, mcPOT);
      
	//Plot(*simEventRate, "simulatedCrossSection", prefix);
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
