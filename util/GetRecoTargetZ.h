//File: GetRecoTargetZ.h
//Brief: Get the reco. target section from vertex location.
//TODO: Add in the material sections as well.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef UTIL_GETRECOTGTZ_H
#define UTIL_GETRECOTGTZ_H

#include "PlotUtils/TargetUtils.h"

namespace util
{
  int GetRecoTargetZ(const double vtx_x, const double vtx_y, const double vtx_z){
    //Thickness different for different materials, so fix this later...
    double Tgt1Lo = PlotUtils::TargetUtils::Get().GetTarget1CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt1::Pb/2;
    double Tgt2Lo = PlotUtils::TargetUtils::Get().GetTarget2CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt2::Pb/2;
    double Tgt3Lo = PlotUtils::TargetUtils::Get().GetTarget3CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt3::Pb/2;
    double Tgt4Lo = PlotUtils::TargetUtils::Get().GetTarget4CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt4::Pb/2;
    double Tgt5Lo = PlotUtils::TargetUtils::Get().GetTarget5CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt5::Pb/2;
    double TgtWLo = PlotUtils::TargetProp::WaterTarget::Face;

    if (!PlotUtils::TargetUtils::Get().IsInHexagon(vtx_x,vtx_y)) return -1;
    if (PlotUtils::TargetUtils::Get().InTarget1ZMC(vtx_z, 82)) return 1;
    else if (PlotUtils::TargetUtils::Get().InTarget2ZMC(vtx_z, 82)) return 2; 
    else if (PlotUtils::TargetUtils::Get().InTarget3ZMC(vtx_z, 6)) return 3;
    else if (PlotUtils::TargetUtils::Get().InTarget4ZMC(vtx_z, 82)) return 4;
    else if (PlotUtils::TargetUtils::Get().InTarget5ZMC(vtx_z, 82)) return 5;
    else if (PlotUtils::TargetUtils::Get().InWaterTargetZMC(vtx_z)) return 6;
    else if (vtx_z < Tgt1Lo) return 10;
    else if (vtx_z < Tgt2Lo) return 21;
    else if (vtx_z < Tgt3Lo) return 32;
    else if (vtx_z < TgtWLo) return 63;
    else if (vtx_z < Tgt4Lo) return 46;
    else if (vtx_z < Tgt5Lo) return 54;
    else return 0;
  }
}

#endif //UTIL_GETRECOTGTZ_H
