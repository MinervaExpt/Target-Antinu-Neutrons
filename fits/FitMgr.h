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
#include "fits/Fit.h"

namespace fit{
  class FitMgr: public ROOT::Math::IBaseFunctionMultiDimTempl<double>{
  private:
    //Members to hold the histos
    const std::vector<fit::Fit*> fFits;
    //std::vector<TH1D*> fUnfitHists; Removing to create an Fit-derived class which isn't fit, but will clena up the structure to match up checks between the fit and non-fit sets of histos.
    const std::vector<TH1D*> fDataHists;
    //Members for fit range
    int fFirstBin;
    int fLastBin;
    //Track of the number of regions
    int fNRegions;
    //Check if doing the overall fit at all.
    bool fDoFit;

  public:
    //CTOR
    FitMgr(const std::vector<fit::Fit*> fits, const std::vector<TH1D*> dataHist);

    unsigned int NDim() const override;

    //Function which the ROOT fitter will minimize
    double DoEval(const double* parameters) const override;

    //Check if everything is good to fit.
    bool CheckFit();

    //Required for ROOT fittable function base class :( (indeed this is sad, Andrew)
    IBaseFunctionMultiDimTempl<double>* Clone() const override;

    //DTOR
    virtual ~FitMgr() = default;
  };
}

#endif
