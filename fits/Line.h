//File: Line.h
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef LINE_H
#define LINE_H

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
  class Line: public Fit{
  public:
    //CTOR
    Line(const std::vector<TH1D*> fitHists, TString name, const int firstBin = 1, const int lastBin = -1);

    unsigned int NDim() const override;

    //Function which the ROOT fitter will minimize
    double GetFitVal(const double* parameters, int whichParam, int whichBin, double extVal0 = -999.0) const override;

    double GetFitErr(const double* parameters, const double* errors, int whichParam, int whichBin, double extErr0 = -999.0) const override;

    //DTOR
    ~Line() = default;
  };
}

#endif
