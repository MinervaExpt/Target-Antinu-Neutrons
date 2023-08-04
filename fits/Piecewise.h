//File: Piecewise.h
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef PIECEWISE_H
#define PIECEWISE_H

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
#include "fits/Fit.h"

namespace fit{
  class Piecewise: public Fit{
  private:
    std::vector<fit::Fit*> fFits;
    int fNDim;
  public:
    //CTOR
    Piecewise(const std::vector<fit::Fit*> fits, const std::vector<TH1D*> fitHists, TString name, const int firstBin = 1, const int lastBin = -1);

    unsigned int NDim() const override;

    //Function which the ROOT fitter will minimize
    double GetFitVal(const double* parameters, int whichParam, int whichBin, double extVal0 = -999.0) const override;

    double GetFitErr(const double* parameters, const double* errors, int whichParam, int whichBin, double extErr0 = -999.0) const override;

    double RecurseGetVal0(const double* parameters, int whichParam, int iFit) const;

    double RecurseGetErr0(const double* parameters, const double* errors, int whichParam, int iFit) const;

    int FindFit(int whichBin) const;
    
    //DTOR
    ~Piecewise() = default;
  };
}

#endif
