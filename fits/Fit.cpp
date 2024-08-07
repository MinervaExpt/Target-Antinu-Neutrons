//File: Fit.cpp
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/Fit.h"

namespace fit{
  bool compFitStartBins(Fit* fit1, Fit* fit2){
    return (fit1->GetFirstFitBin() < fit2->GetFirstFitBin());
  }
  
  Fit::Fit(const std::vector<TH1D*> fitHists, TString name, int firstBin, int lastBin): fFitHists(fitHists), fNBins(0), fName(name), fFirstBin(firstBin), fLastBin(lastBin), fNRegions(0), fDoFit(true), fExtVal0(false)
  {
    if (fFitHists.size() == 0){
      std::cout << "Function has no histos to fit... Setting the fit condition to false." << std::endl;
      fDoFit = false;
    }

    int nBins = 0;
    if (fDoFit) nBins = fFitHists.at(0)->GetNbinsX();

    //Do I want to ensure that the histos are the same exact bins too? For now, I'll leave it as is.
    for(unsigned int iFit; iFit < fFitHists.size(); ++iFit){
      if (!fDoFit) break;
      if (fFitHists.at(iFit)->GetNbinsX() != nBins){
	std::cout << "No. of bins not consistent... Setting the fit condition to false." << std::endl;
	fDoFit = false;
      }
    }

    //Avoid fitting overflow, and set NBins value.
    if (fDoFit){
      fNBins = nBins;
      fNRegions = fFitHists.size();
      if (fLastBin < 0 || fLastBin > nBins) fLastBin=nBins;
    }
  }

  unsigned int Fit::NDim() const{
    return 0;//This base class should not be used without a specified fit function to be overwritten in derived classes, which also therefore define the number of relevant parameters.
  }
  
  double Fit::GetFitVal(const double* parameters, int whichParam, int whichBin, double extVal0) const{
    return -999.0;//This base class should not be used without a specified fit function to be overwritten in derived classes.
  }

  double Fit::GetFitErr(const double* parameters, const double* errors, int whichParam, int whichBin, double extErr0) const{
    return -999.0;//This base class should not be used without a specified fit function to be overwritten in derived classes.
  }

  std::vector<double> Fit::GetVals(const double* parameters, int whichParam, int whichBin) const{
    std::vector<double> values;
    if (fDoFit){
      double fitVal = GetFitVal(parameters,whichParam,whichBin);
      for (auto hist:fFitHists){
	double value = hist->GetBinContent(whichBin)*fitVal;
        values.push_back(value);
      }
    }
    return values;
  }

  int Fit::GetNBins() { return fNBins; }
  int Fit::GetFirstFitBin() { return fFirstBin; }
  int Fit::GetLastFitBin() { return fLastBin; }

  int Fit::GetNRegions() { return fNRegions; }

  bool Fit::CheckFit() { return fDoFit; }
  
  TString Fit::FitName() { return fName; }

  void Fit::SetExtVal0(bool setVal) { fExtVal0 = setVal; }
  bool Fit::CheckExtVal0() { return fExtVal0; }
}
