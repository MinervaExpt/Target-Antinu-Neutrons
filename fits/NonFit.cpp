//File: NonFit.cpp
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/NonFit.h"

namespace fit{
  NonFit::NonFit(const std::vector<TH1D*> fitHists, int firstBin, int lastBin):Fit(fitHists, firstBin, lastBin)
  {
  }

  unsigned int NonFit::NDim() const{
    return 0;
  }
  
  std::vector<double> NonFit::GetVals(const double* parameters, int whichParam, int whichBin) const{
    std::vector<double> values;
    if (fDoFit){
      for (auto hist:fFitHists){
	double value = hist->GetBinContent(whichBin);
	values.push_back(value);
      }
    }
    return values;
  }

}
