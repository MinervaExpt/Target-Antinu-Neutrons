//File: Line.cpp
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/Line.h"

namespace fit{
  Line::Line(const std::vector<TH1D*> fitHists, TString name, int firstBin, int lastBin):Fit(fitHists, name, firstBin, lastBin)
  {
  }

  unsigned int Line::NDim() const{
    if (!fDoFit) return 0;
    else if (fExtVal0) return 1;
    else return 2;
  }
  
  double Line::GetFitVal(const double* parameters, int whichParam, int whichBin, double extVal0) const{
    double fitVal = -999.0;
    if (fDoFit){
      double lowCenter = fFitHists.at(0)->GetBinCenter(fFirstBin);
      double totDiff = fFitHists.at(0)->GetBinCenter(fLastBin)-lowCenter;
      double binDiff = fFitHists.at(0)->GetBinCenter(whichBin)-lowCenter;
      double range = fFitHists.at(0)->GetBinCenter(fLastBin)-lowCenter;
      double val0 = (fExtVal0) ? extVal0 : (parameters+whichParam)[0];
      double val1 = (fExtVal0) ? (parameters+whichParam)[0] : (parameters+whichParam)[1];
      fitVal = val0 + (val1-val0)*(binDiff/totDiff);
    }
    return fitVal;
  }

  double Line::GetFitErr(const double* parameters, const double* errors, int whichParam, int whichBin, double extErr0) const{
    double fitErr = -999.0;
    if (fDoFit){
      double lowCenter = fFitHists.at(0)->GetBinCenter(fFirstBin);
      double totDiff = fFitHists.at(0)->GetBinCenter(fLastBin)-lowCenter;
      double binDiff = fFitHists.at(0)->GetBinCenter(whichBin)-lowCenter;
      double err0 = (fExtVal0) ? extErr0 : (errors+whichParam)[0];
      double err1 = (fExtVal0) ? (errors+whichParam)[0] : (errors+whichParam)[1];
      double var = (1-(binDiff/totDiff))*(1-(binDiff/totDiff))*err0*err0 + (binDiff/totDiff)*(binDiff/totDiff)*err1*err1;
      fitErr = std::pow(var, 0.5);
    }
    return fitErr;
  }
}
