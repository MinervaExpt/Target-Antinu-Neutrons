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
  
  double NonFit::GetFitVal(const double* parameters, int whichParam, int whichBin) const{
    double fitVal = 0.0;
    if (fDoFit) fitVal = 1.0;
    return fitVal;
  }

}
