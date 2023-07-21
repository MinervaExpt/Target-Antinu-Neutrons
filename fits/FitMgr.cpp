//File: FitMgr.cpp
//Info: Class for handling different fit functions with TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop 
//      as of 07/21/2023.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/FitMgr.h"

namespace fit{
  FitMgr::FitMgr(const std::vector<fit::Fit*> fits, const std::vector<TH1D*> unfitHists, const TH1D* dataHist):IBaseFunctionMultiDimTempl<double>(), fFits(fits), 
								    fUnfitHists(unfitHists), fDataHist(dataHist), fFirstBin(1), fLastBin(-1), fDoFit(true)
  {
    if (fFits.size() == 0){
      std::cout << "Function has no histos to fit... Setting the fit condition to false." << std::endl;
      fDoFit = false;
    }

    int nBins = 0;
    if (fDoFit) nBins = fFits.at(0)->GetNBins();

    for(unsigned int iFit=0; iFit < fFits.size(); ++iFit){
      if (!fDoFit) break;
      if (fFits.at(iFit)->GetNBins() != nBins){
	std::cout << "No. of bins not consistent... Setting the fit condition to false." << std::endl;
	fDoFit = false;
      }
    }

    for(unsigned int iFit=0; iFit < fUnfitHists.size(); ++iFit){
      if (!fDoFit) break;
      if (fUnfitHists.at(iFit)->GetNbinsX() != nBins){
	std::cout << "No. of bins not consistent... Setting the fit condition to false." << std::endl;
	fDoFit = false;
      }
    }

    if (fDataHist->GetNbinsX() != nBins){
      std::cout << "No. of bins not consistent... Setting the fit condition to false." << std::endl;
      fDoFit = false;
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

    //Avoid fitting overflow
    if (fDoFit && (fLastBin < 0 || fLastBin > nBins)) fLastBin=nBins;
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
    //For now just copying in the chi^2 function more or less directly from Andrew
    double chi2 = 0.0;
    if (!fDoFit) return chi2;
    
    for (int whichBin = fFirstBin; whichBin <= fLastBin; ++whichBin){
      int whichParam = 0;
      double fitSum = 0.0;

      for(unsigned int whichFit=0; whichFit < fFits.size(); ++whichFit){
	fitSum += fFits.at(whichFit)->GetVal(parameters, whichParam, whichBin);
	whichParam += fFits.at(whichFit)->NDim();
      }

      for (unsigned int whichUnfit=0; whichUnfit < fUnfitHists.size(); ++whichUnfit){
	fitSum += fUnfitHists.at(whichUnfit)->GetBinContent(whichBin);
      }

      double dataContent = fDataHist->GetBinContent(whichBin);
      double dataErr = fDataHist->GetBinError(whichBin);
      double diff = fitSum-dataContent;
      //std::cout << "Fit Sum: " << fitSum << ", Data: " << dataContent << ", Difference: " << diff << ", Error: " << dataErr << std::endl;
      //std::cout << "How the chi2 should change: " << (diff*diff)/(dataErr*dataErr) << std::endl;
      if (dataErr > 1e-10) chi2 += (diff*diff)/(dataErr*dataErr);
      //std::cout << "Updated chi2: " << chi2 << std::endl;
    }

    //std::cout << "About to return chi2 of: " << chi2 << std::endl;

    return chi2;
  }

  //Required for ROOT fittable function base class :( (indeed this is sad, Andrew)
  ROOT::Math::IBaseFunctionMultiDimTempl<double>* FitMgr::Clone() const{
    return new FitMgr(fFits, fUnfitHists, fDataHist);
  }

}
