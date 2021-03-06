//File: GetBackgroundID.h
//Brief: Get the CCQENu Background ID from the universe at hand
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef UTIL_GETBKGID_H
#define UTIL_GETBKGID_H

#include "event/CVUniverse.h"

namespace util
{
  int GetBackgroundID(const CVUniverse& univ){
    int nPiC = 0;
    int nPi0 = 0;
    std::vector<int> FSPDGs = univ.GetFSPartPDG();
    for (int i=0; i < FSPDGs.size(); ++i){
      if (abs(FSPDGs.at(i)) == 211) nPiC++;
      else if (FSPDGs.at(i) == 111) nPi0++;
    }
    if ((nPiC + nPi0) > 1){
      return 3;
    }
    else if (nPiC == 1){
      return 1;
    }
    else if (nPi0 == 1){
      return 2;
    }
    else return 0;
  }
}

#endif //UTIL_GETINGREDIENT_H
