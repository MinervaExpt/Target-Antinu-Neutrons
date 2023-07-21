//File: ScaleFactor.cpp
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/ScaleFactor.h"

namespace fit{
  ScaleFactor::ScaleFactor(const TH1D* fitHist, int firstBin, int lastBin):Fit(fitHist, firstBin, lastBin)
  {
  }

  unsigned int ScaleFactor::NDim() const{
    if (!fDoFit) return 0;
    else return 1;
  }
  
  double ScaleFactor::GetVal(const double* parameters, int whichParam, int whichBin) const{
    double value = fFitHist->GetBinContent(whichBin);
    if (fDoFit) value = value*((parameters+whichParam)[0]);//Previously did nothing, fixed now to be the returned value only when fitting.
    return value;
  }

}
