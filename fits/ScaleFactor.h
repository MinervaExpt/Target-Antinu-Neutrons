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
#include "fits/Fit.h"

namespace fit{
  class ScaleFactor: public Fit{
  public:
    //CTOR
    ScaleFactor(const TH1D* fitHist, const int firstBin = 1, const int lastBin = -1);

    virtual unsigned int NDim() const override;

    //Function which the ROOT fitter will minimize
    virtual double GetVal(const double* parameters, int whichParam, int whichBin) const override;

    //DTOR
    virtual ~ScaleFactor() = default;
  };
}

#endif
