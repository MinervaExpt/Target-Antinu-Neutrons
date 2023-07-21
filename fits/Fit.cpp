//File: Fit.cpp
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/Fit.h"

namespace fit{
  Fit::Fit(const TH1D* fitHist, int firstBin, int lastBin): fFitHist(fitHist), fFirstBin(firstBin), fLastBin(lastBin), fDoFit(true)
  {
    //Need to replace the checks for bins and stuff somewhere in the new structure. 
    //TODO: When this holds different fit regions, need to check that the different regions have the same bins and such.

    /*
    if (fFitHists.size() == 0){
      std::cout << "Function has no histos to fit... Setting the fit condition to false." << std::endl;
      fDoFit = false;
    }

    int nBins = 0;
    if (fDoFit) nBins = fFitHists.at(0)->GetNbinsX();

    for(unsigned int iFit; iFit < fFitHists.size(); ++iFit){
      if (!fDoFit) break;
      if (fFitHists.at(iFit)->GetNbinsX() != nBins){
	std::cout << "No. of bins not consistent... Setting the fit condition to false." << std::endl;
	fDoFit = false;
      }
    }

    for(unsigned int iFit; iFit < fUnfitHists.size(); ++iFit){
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

    //Avoid fitting overflow
    if (fDoFit && fLastBin < 0 || fLastBin > nBins) fLastBin=nBins;
    */
  }

  unsigned int Fit::NDim() const{
    return 0;//This base class should not be used without a specified fit function to be overwritten in derived classes, which also therefore define the number of relevant parameters.
  }
  
  double Fit::GetVal(const double* parameters, int whichParam, int whichBin) const{
    return -999;//This base class should not be used without a specified fit dunction to be overwritten in derived classes.
  }

  int Fit::GetNBins() { return fFitHist->GetNbinsX(); }

  int Fit::GetFirstFitBin() { return fFirstBin; }
  int Fit::GetLastFitBin() { return fLastBin; }

  bool Fit::CheckFit() { return fDoFit; }
  void Fit::SetFit(bool doFit){ fDoFit = doFit; }

}
