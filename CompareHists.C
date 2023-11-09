//File: CompareHists.C
//Info: Quick macro to load to compare any two MnvH1D's passed by pointer.
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

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"

using namespace std;
using namespace PlotUtils;

void CompareHists(const MnvH1D* m1, const MnvH1D* m2, double tolerance) {
  cout << "First check: See if they have the same number of bins" << endl;
  if (m1->GetNbinsX() != m2->GetNbinsX()){
    cout << "Not same bins. Not comparing." << endl;
    return;
  }

  cout << "Second check: See if all of the contents and total errors are within the tolerance." << endl;

  /*
  TH1D* h1 = (TH1D*)(m1->GetCVHistoWithError().Clone());
  TH1D* h2 = (TH1D*)(m2->GetCVHistoWithError().Clone());
  */

  int nBins = m1->GetNbinsX();
  for (int iBin=0; iBin <= nBins+1; ++iBin){
    if (abs(m1->GetBinContent(iBin) - m2->GetBinContent(iBin)) > tolerance){
      cout << "Bin: " << iBin << " has content: " << m1->GetBinContent(iBin) << " in h1, and content: " << m2->GetBinContent(iBin) << " in h2." << endl;
      return;
    }
    if (abs(m1->GetBinError(iBin) - m2->GetBinError(iBin)) > tolerance){
      cout << "Bin: " << iBin << " has error: " << m1->GetBinError(iBin) << " in h1, and error: " << m2->GetBinError(iBin) << " in h2." << endl;
      return;
    }
  }

  cout << "Seemingly the same." << endl;
  return;
}
