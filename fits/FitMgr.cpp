//File: FitMgr.cpp
//Info: Class for handling different fit functions with TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop 
//      as of 07/21/2023.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/FitMgr.h"

namespace fit{
  FitMgr::FitMgr(const std::vector<fit::Fit*> fits, const std::vector<TH1D*> dataHists):IBaseFunctionMultiDimTempl<double>(), fFits(fits), 
											fDataHists(dataHists), fFirstBin(1), fLastBin(-1), fNRegions(0), fDoFit(true)
  {
    if (fFits.size() == 0){
      std::cout << "Function has no histos to fit... Setting the fit condition to false." << std::endl;
      fDoFit = false;
    }

    int nBins = 0;
    int nRegions = 0;
    if (fDoFit){
      nBins = fFits.at(0)->GetNBins();
      nRegions = fFits.at(0)->GetNRegions();
    }

    if (nRegions != fDataHists.size()){
      std::cout << "No. of data hists provided does not match the number of regions expected to be fit." << std::endl;
      fDoFit = false;
    }

    for(unsigned int iFit=0; iFit < fFits.size(); ++iFit){
      if (!fDoFit) break;
      if (fFits.at(iFit)->GetNBins() != nBins || !fFits.at(iFit)->CheckFit()){
	std::cout << "No. of fit bins not consistent or fit says it isn't ready... Setting the fit condition to false." << std::endl;
	fDoFit = false;
      }
    }

    for (unsigned int iRegion=0; iRegion < fDataHists.size(); ++iRegion){
      if (!fDoFit) break;
      if (fDataHists.at(iRegion)->GetNbinsX() != nBins){
	std::cout << "No. of data bins not consistent... Setting the fit condition to false." << std::endl;
	fDoFit = false;
      }
    }

    //Getting the fit range from the provided fit objects, and checking that they are the same.
    if (fDoFit){
      fFirstBin = fFits.at(0)->GetFirstFitBin();
      fLastBin = fFits.at(0)->GetLastFitBin();
      for(unsigned int iFit=0; iFit < fFits.size(); ++iFit){
	if (!fDoFit) break;
	if (fFits.at(iFit)->GetFirstFitBin() != fFirstBin || fFits.at(iFit)->GetLastFitBin() != fLastBin){
	  std::cout << "No. of fit bins not consistent among fits... Setting the fit condition to false." << std::endl;
	  fDoFit = false;
	}
      }
    }

    //Avoid fitting overflow, this should never do anything since this is checked and changed for every Fit.
    if (fDoFit) {
      if (fLastBin < 0 || fLastBin > nBins) fLastBin=nBins;
      fNRegions = nRegions;
    }
  }

  unsigned int FitMgr::NDim() const{
    if (!fDoFit) return 0;
    unsigned int nDimSum = 0;
    for (unsigned int iFit=0; iFit < fFits.size(); ++iFit){
      nDimSum += fFits.at(iFit)->NDim();
    }
    return nDimSum;
  }
  
  double FitMgr::DoEval(const double* parameters) const{
    //For now, just copying in the chi^2 function more or less directly from Andrew
    //std::cout << "Bin lo in DoEval: " << fFirstBin << std::endl;
    //std::cout << "Bin hi in DoEval: " << fLastBin << std::endl;
    double chi2 = 0.0;
    if (!fDoFit){
      std::cout << "Not fitting..." << std::endl;
      return chi2;
    }
    
    for (int whichBin = fFirstBin; whichBin <= fLastBin; ++whichBin){
      //std::cout << "whichBin Global: " << whichBin << std::endl; 
      int whichParam = 0;
      std::vector<double> fitSums(fNRegions, 0.0);

      for(unsigned int whichFit=0; whichFit < fFits.size(); ++whichFit){
	std::vector<double> values = fFits.at(whichFit)->GetVals(parameters, whichParam, whichBin);
	if (values.size() != fNRegions){
	  //std::cout << "Issue getting vals!" << std::endl;
	  //std::cout << "whichFit: " << whichFit << std::endl;
	  //std::cout << "only got: " << values.size() << std::endl;
	  return 0.0;
	}
	for(int iRegion=0; iRegion < fNRegions; ++iRegion) fitSums.at(iRegion) += values.at(iRegion);
	whichParam += fFits.at(whichFit)->NDim();
      }

      for (int iRegion=0; iRegion < fNRegions; ++iRegion){
	double dataContent = fDataHists.at(iRegion)->GetBinContent(whichBin);
	double dataErr = fDataHists.at(iRegion)->GetBinError(whichBin);
	double diff = fitSums.at(iRegion)-dataContent;
	//std::cout << "Fit Sum: " << fitSums.at(iRegion) << ", Data: " << dataContent << ", Difference: " << diff << ", Error: " << dataErr << std::endl;
	//std::cout << "How the chi2 should change: " << (diff*diff)/(dataErr*dataErr) << std::endl;
	//std::cout << "chi2 before update: " << std::endl;
	if (dataErr > 1e-10) chi2 += (diff*diff)/(dataErr*dataErr);
	//std::cout << "Updated chi2: " << chi2 << std::endl;
      }
    }

    //std::cout << "About to return chi2 of: " << chi2 << std::endl;

    return chi2;
  }

  //Check if everything is good to fit.
  bool FitMgr::CheckFit() { return fDoFit; }

  //Required for ROOT fittable function base class :( (indeed this is sad, Andrew)
  ROOT::Math::IBaseFunctionMultiDimTempl<double>* FitMgr::Clone() const{
    return new FitMgr(fFits, fDataHists);
  }

}
