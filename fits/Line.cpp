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
  
  std::vector<double> Line::GetVals(const double* parameters, int whichParam, int whichBin) const{
    std::vector<double> values;
    if (fDoFit){
      for (auto hist:fFitHists){
	double binCenter = hist->GetBinCenter(whichBin);
	double value = hist->GetBinContent(whichBin)*((parameters+whichParam)[0]+(parameters+whichParam)[1]*(binCenter));
	values.push_back(value);
      }
    }
    return values;
  }

}
