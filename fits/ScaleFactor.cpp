//File: ScaleFactor.cpp
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/ScaleFactor.h"

namespace fit{
  ScaleFactor::ScaleFactor(const std::vector<TH1D*> fitHists, TString name, int firstBin, int lastBin):Fit(fitHists, name, firstBin, lastBin)
  {
  }

  unsigned int ScaleFactor::NDim() const{
    if (!fDoFit || fExtVal0) return 0;
    else return 1;
  }
  
  double ScaleFactor::GetFitVal(const double* parameters, int whichParam, int whichBin, double extVal0) const{
    double fitVal = -999.0;
    if (fDoFit) fitVal = (fExtVal0) ? extVal0 : (parameters+whichParam)[0];
    return fitVal;
  }

  double ScaleFactor::GetFitErr(const double* parameters, const double* errors, int whichParam, int whichBin, double extErr0) const{
    double fitErr = -999.0;
    if (fDoFit) fitErr = (fExtVal0) ? extErr0 : (errors+whichParam)[0];
    return fitErr;
  }

}
