//File: Line.cpp
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/Line.h"

namespace fit{
  Line::Line(const std::vector<TH1D*> fitHists, int firstBin, int lastBin):Fit(fitHists, firstBin, lastBin)
  {
  }

  unsigned int Line::NDim() const{
    if (!fDoFit) return 0;
    else return 2;
  }
  
  double Line::GetFitVal(const double* parameters, int whichParam, int whichBin) const{
    double fitVal = 0.0;
    if (fDoFit){
      double binCenter = fFitHists.at(0)->GetBinCenter(whichBin);
      fitVal = (parameters+whichParam)[0]+(parameters+whichParam)[1]*(binCenter);
    }
    return fitVal;
  }

}
