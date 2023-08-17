//File: checkPossibleBinning.C
//Info: Takes input migration histo for a 1D variable and calculates bin edges with at least N% on the diagonal.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

//C++ includes
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <bitset>
#include <time.h>
#include <sys/stat.h>

//ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"
#include "TParameter.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH2D.h"

using namespace std;
using namespace PlotUtils;

pair<double,double> combineSums(MnvH2D* h, int currBin, int nextBin, double currSum, double currTotal){
  double retSum = currSum;
  double retTotal = currTotal;
  for (int iBin=1; iBin <= h->GetNbinsX(); ++iBin){
    retTotal += h->GetBinContent(iBin, nextBin);
    if (iBin < nextBin && iBin >= currBin){
      retSum += h->GetBinContent(iBin, nextBin);
    }
    else if (iBin == nextBin){
      for (int iBinY=currBin; iBinY <= nextBin; ++iBinY){
	retSum += h->GetBinContent(iBin, iBinY);
      }
    }
  }
  pair<double,double> retPair = make_pair(retSum, retTotal);
  return retPair;
}

void checkPossibleBinning(MnvH2D* hIN, vector<double> binUpEdges, double nPercent) {
  MnvH2D* h = (MnvH2D*)hIN->Clone();
  cout << "Total Y bins: " << h->GetNbinsY() << endl;
  cout << "Total X bins: " << h->GetNbinsX() << endl;
 
  pair<double,double> sums =  make_pair(0.0, 0.0);
  bool combineBins = false;
  int lastBin = 1;
  int extBin = 1;
  for (int iBinY = 1; iBinY <= h->GetNbinsY(); ++iBinY){
    if (!combineBins){
      sums = combineSums(h,iBinY,iBinY,0.0,0.0);
      lastBin = iBinY;
    }
    else{
      sums = combineSums(h, lastBin, iBinY, sums.first, sums.second);
    }
    //    cout << "Current diagonal bin sum: " << sums.first << ", current row Sum: " << sums.second << endl;
    double percent = (sums.second > 0.0) ? 100.0*(sums.first/sums.second) : 0.0;
    //cout << "Percent: " << percent << endl;
    if (h->GetYaxis()->GetBinUpEdge(iBinY) > binUpEdges.back()){
      cout << "Reached end of proposed binning at bin edge: " << h->GetYaxis()->GetBinUpEdge(iBinY) << ", Leaving script." << endl;
      break;
    }
    else if (fabs(h->GetYaxis()->GetBinUpEdge(iBinY)-binUpEdges.at(extBin)) < 1e-9){
      combineBins = false;
      ++extBin;
      cout << "Combine bin upper bound: " << h->GetYaxis()->GetBinUpEdge(iBinY) << ", with diagonal percent of row: " << percent <<", and total bin content: " << sums.first << endl;
    }
    else{
      combineBins = true;
      //cout << "Not done combining yet." << endl;
    }
  }

  double percent = (sums.second > 0.0) ? 100.0*(sums.first/sums.second) : 0.0;
  cout << "Final percent: " << percent << ", and total bin content: " << sums.first << endl;

  delete h;
}
