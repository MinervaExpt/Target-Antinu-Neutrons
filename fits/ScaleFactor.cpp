//File: ScaleFactor.cpp
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/ScaleFactor.h"

namespace fit{
  ScaleFactor::ScaleFactor(const TH1D* fitHist, int firstBin, int lastBin): fFitHist(fitHist), fFirstBin(firstBin), fLastBin(lastBin), fDoFit(true)
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

  unsigned int ScaleFactor::NDim() const{
    if (!fDoFit) return 0;
    else return 1;
  }
  
  double ScaleFactor::GetVal(const double* parameters, int whichParam, int whichBin) const{
    double value = fFitHist->GetBinContent(whichBin)*((parameters+whichParam)[0]);
    if (fDoFit) value*((parameters+whichParam)[0]);
    return value;
  }

  int ScaleFactor::GetNBins() { return fFitHist->GetNbinsX(); }

  int ScaleFactor::GetFirstFitBin() { return fFirstBin; }
  int ScaleFactor::GetLastFitBin() { return fLastBin; }

  void ScaleFactor::SetFit(bool doFit){ fDoFit = doFit; }

}
