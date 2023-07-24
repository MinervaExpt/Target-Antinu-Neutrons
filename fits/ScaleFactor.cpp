//File: ScaleFactor.cpp
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/ScaleFactor.h"

namespace fit{
  ScaleFactor::ScaleFactor(const std::vector<TH1D*> fitHists, int firstBin, int lastBin):Fit(fitHists, firstBin, lastBin)
  {
  }

  unsigned int ScaleFactor::NDim() const{
    if (!fDoFit) return 0;
    else return 1;
  }
  
  std::vector<double> ScaleFactor::GetVals(const double* parameters, int whichParam, int whichBin) const{
    std::vector<double> values;
    if (fDoFit){
      for (auto hist:fFitHists){
	double value = hist->GetBinContent(whichBin)*((parameters+whichParam)[0]);
	values.push_back(value);
      }
    }
    return values;
  }

}
