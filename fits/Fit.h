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
    const std::vector<TH1D*> fFitHists;
    //Members for fit range
    int fNBins;
    int fFirstBin;
    int fLastBin;
    //Quick Access to number of regions
    int fNRegions;
    //Check if you are using this fit
    bool fDoFit;

    //CTOR
    Fit(const std::vector<TH1D*> fitHists, const int firstBin = 1, const int lastBin = -1);

    virtual unsigned int NDim() const;

    //Function which the ROOT fitter will minimize
    virtual double GetFitVal(const double* paramters, int whichParam, int whichBin) const;

    virtual double GetFitErr(const double* paramters, const double* errors, int whichParam, int whichBin) const;

    //Function which gets bin content and scales it by the fit function value from GetFitVal(...)
    std::vector<double> GetVals(const double* parameters, int whichParam, int whichBin) const;

    int GetNBins();
    int GetFirstFitBin();
    int GetLastFitBin();

    int GetNRegions();

    bool CheckFit();

    //DTOR
    virtual ~Fit() = default;
  };
}

#endif
