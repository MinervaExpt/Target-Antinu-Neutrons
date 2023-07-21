//File: Fit.h
//Info: Base class for defining a fit function to pass to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef FIT_H
#define FIT_H

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
  class Fit{
    //  private:
    //
  public:
    //Members to hold the histos
    const TH1D* fFitHist;
    //Members for fit range
    int fFirstBin;
    int fLastBin;
    bool fDoFit;

    //CTOR
    Fit(const TH1D* fitHist, const int firstBin = 1, const int lastBin = -1);

    virtual unsigned int NDim() const;

    //Function which the ROOT fitter will minimize
    virtual double GetVal(const double* parameters, int whichParam, int whichBin) const;

    int GetNBins();
    int GetFirstFitBin();
    int GetLastFitBin();

    bool CheckFit();
    void SetFit(bool doFit);

    //DTOR
    virtual ~Fit() = default;
  };
}

#endif
