//File: FitMgr.h
//Info: Class for handling varied fit functions to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 07/21/2023.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef FITMGR_H
#define FITMGR_H

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
#include "fits/ScaleFactor.h"

namespace fit{
  class FitMgr: public ROOT::Math::IBaseFunctionMultiDimTempl<double>{
  private:
    //Members to hold the histos
    std::vector<fit::ScaleFactor*> fFits;
    std::vector<TH1D*> fUnfitHists;
    const TH1D* fDataHist;
    //Members for fit range
    int fFirstBin;
    int fLastBin;
    bool fDoFit;

  public:
    //CTOR
    FitMgr(const std::vector<fit::ScaleFactor*> fits, const std::vector<TH1D*> unfitHists, const TH1D* dataHist);

    unsigned int NDim() const override;

    //Function which the ROOT fitter will minimize
    double DoEval(const double* parameters) const override;

    //Required for ROOT fittable function base class :( (indeed this is sad, Andrew)
    IBaseFunctionMultiDimTempl<double>* Clone() const override;

    //DTOR
    virtual ~FitMgr() = default;
  };
}

#endif
