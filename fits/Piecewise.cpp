//File: Piecewise.cpp
//Info: Class for passing as a function to TMinuit for minimization.
//      Largely derived by example from Andrew Olivier's 
//      github.com/MinervaExpt/NucCCNeutrons/tree/develop/fits/Universe.h 
//      as of 03/16/2022.
//
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#include "fits/Piecewise.h"

namespace fit{
  Piecewise::Piecewise(const std::vector<fit::Fit*> fits, const std::vector<TH1D*> fitHists, TString name, int firstBin, int lastBin):Fit(fitHists, name, firstBin, lastBin), fFits(fits)
  {
    //TO-DO: Add a check on the histograms among the fits to ensure they are the same.
    int lowBin = 0;
    int hiBin = 0;
    int prevLastBin = 0;
    if (fFits.size()==0){
      std::cout << "Cannot fit a piecewise function of no pieces." << std::endl;
      fDoFit = false;
    }
    else{
      lowBin = fFits.at(0)->GetFirstFitBin();
      hiBin = fFits.at(0)->GetLastFitBin();
      prevLastBin = hiBin;
      for (unsigned int iFit=1; iFit < fFits.size(); ++iFit){
	int bin = fFits.at(iFit)->GetFirstFitBin();
	if (bin == prevLastBin){
	  fFits.at(iFit)->SetExtVal0(true);
	}
	else if (bin == prevLastBin+1){
	  fFits.at(iFit)->SetExtVal0(false);
	}
	else{
	  std::cout << "Not fitting piecewise because pieces do not overlap correctly." << std::endl;
	  fDoFit=false;
	  break;
	}
	hiBin = fFits.at(iFit)->GetLastFitBin();
	prevLastBin = hiBin;
      }
    }
    //std::cout << "At initialization of piecewise" << std::endl;
    //std::cout << lowBin << ", " << fFirstBin << std::endl;
    //std::cout << hiBin << ", " << fLastBin << std::endl;
    if (lowBin != fFirstBin || hiBin != fLastBin){
      std::cout << "Not fitting piecewise because combined fit range does not match declared fit range" << std::endl;
      fDoFit = false;
    }
  }

  unsigned int Piecewise::NDim() const{
    unsigned int nDimSum =0;
    if (fDoFit){
      for (auto fit:fFits){
	nDimSum += fit->NDim();
      }
    }
    return nDimSum;
  }
  
  double Piecewise::GetFitVal(const double* parameters, int whichParam, int whichBin, double extVal0) const{
    //std::cout << "In piecewise getfitval" << std::endl;
    double fitVal = -999.0;
    if (fDoFit){
      //std::cout << "whichBin: " << whichBin << std::endl;
      //std::cout << "whichParam: " << whichParam << std::endl;
      int whichFit = FindFit(whichBin); 
      //std::cout << "which Fit: " << iFit << std::endl;
      int theParam = whichParam;
      for (int iFit = 0; iFit < whichFit; ++iFit) theParam += fFits.at(iFit)->NDim();
      double val0 = RecurseGetVal0(parameters, theParam, whichFit);
      //std::cout << "0 val: " << val0 << std::endl;
      fitVal = fFits.at(whichFit)->GetFitVal(parameters, theParam, whichBin, val0);
      //std::cout << "fitVal: " << fitVal << std::endl;
    }
    //std::cout << "fitVal for fit: " << fName << ", " << fitVal << std::endl;
    return fitVal;
  }

  double Piecewise::GetFitErr(const double* parameters, const double* errors, int whichParam, int whichBin, double extErr0) const{
    double fitErr = -999.0;
    if (fDoFit){
      int whichFit = FindFit(whichBin);
      int theParam = whichParam;
      for (int iFit = 0; iFit < whichFit; ++iFit) theParam += fFits.at(iFit)->NDim();
      double err0 = RecurseGetErr0(parameters, errors, theParam, whichFit);
      fitErr = fFits.at(whichFit)->GetFitErr(parameters, errors, theParam, whichBin, err0);
    }
    return fitErr;
  }

  double Piecewise::RecurseGetVal0(const double* parameters, int whichParam, int iFit) const{
    if (iFit <= 0 || !fDoFit) return -999.0;
    double val0 = -999.0;
    if (fFits.at(iFit-1)->CheckExtVal0()){
      val0 = RecurseGetVal0(parameters, whichParam - fFits.at(iFit-1)->NDim(), iFit-1);
    }
    double val = fFits.at(iFit-1)->GetFitVal(parameters, whichParam - fFits.at(iFit-1)->NDim(), fFits.at(iFit)->GetFirstFitBin(), val0);
    return val;
  }

  double Piecewise::RecurseGetErr0(const double* parameters, const double* errors, int whichParam, int iFit) const{
    if (iFit <= 0 || !fDoFit) return -999.0;
    double err0 = -999.0;
    if (fFits.at(iFit-1)->CheckExtVal0()){
      err0 = RecurseGetErr0(parameters, errors, whichParam-fFits.at(iFit-1)->NDim(), iFit-1);
    }
    double err = fFits.at(iFit-1)->GetFitErr(parameters, errors, whichParam-fFits.at(iFit-1)->NDim(), fFits.at(iFit)->GetFirstFitBin(), err0);
    return err;
  }

  int Piecewise::FindFit(int whichBin) const{
    int whichFit = 0;
    for (unsigned int iFit=0; iFit < fFits.size(); ++iFit){
      //std::cout << "iFit: " << iFit << ", last bin: " << fFits.at(iFit)->GetLastFitBin() << std::endl;
      if (whichBin > fFits.at(iFit)->GetLastFitBin()) ++whichFit;
      //std::cout << "whichFit: " << whichFit << std::endl;
    }
    if (whichFit == fFits.size()) return -1;
    return whichFit;
  }
}
