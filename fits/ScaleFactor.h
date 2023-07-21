//File: ScaleFactor.h
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef SCALEFACTOR_H
#define SCALEFACTOR_H

#include "TH1D.h"
#include "TVector3.h"
#include "Math/IFunction.h"
#include "stdlib.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <bitset>

namespace fit{
  class ScaleFactor{
  private:
    //Members to hold the histos
    const TH1D* fFitHist;
    //Members for fit range
    int fFirstBin;
    int fLastBin;
    bool fDoFit;

  public:
    //CTOR
    ScaleFactor(const TH1D* fitHist, const int firstBin = 1, const int lastBin = -1);

    unsigned int NDim() const;

    //Function which the ROOT fitter will minimize
    double GetVal(const double* parameters, int whichParam, int whichBin) const;

    int GetNBins();
    int GetFirstFitBin();
    int GetLastFitBin();

    void SetFit(bool doFit);

    //DTOR
    virtual ~ScaleFactor() = default;
  };
}

#endif
