//File: GetRecoTargetZ.h
//Brief: Get the reco. target section from vertex location.
//       This has expanded to do all that I need for targets and could use a 
//       better name as a result.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef UTIL_GETRECOTGTZ_H
#define UTIL_GETRECOTGTZ_H

#include "PlotUtils/TargetUtils.h"

namespace util
{
  std::map<int,double> TgtDSCut = {{1,4523.0},
			      {2,4744.0},
			      {3,5009.0},
			      {4,5686.0},
			      {5,5819.0},
			      {6,5465.0}};

  std::map<int, double> TgtUSPlaneBack = {{1,4455.0},
					  {2,4676.0},
					  {3,4897.0},
					  {4,5619.0},
					  {5,5751.0},
					  {6,5162.0}};
  
  double step_big = 23.57;
  double step_small = 20.64;

  int GetRecoTargetZ(const double vtx_x, const double vtx_y, const double vtx_z){
    //Thickness different for different materials, so fix this later...
    double Tgt1Lo = PlotUtils::TargetUtils::Get().GetTarget1CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt1::Pb/2;
    double Tgt2Lo = PlotUtils::TargetUtils::Get().GetTarget2CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt2::Pb/2;
    double Tgt3Lo = PlotUtils::TargetUtils::Get().GetTarget3CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt3::Pb/2;
    double Tgt4Lo = PlotUtils::TargetUtils::Get().GetTarget4CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt4::Pb/2;
    double Tgt5Lo = PlotUtils::TargetUtils::Get().GetTarget5CenterZMC() - PlotUtils::TargetProp::ThicknessMC::Tgt5::Pb/2;
    double TgtWLo = PlotUtils::TargetProp::WaterTarget::Face;

    if (!PlotUtils::TargetUtils::Get().IsInHexagon(vtx_x,vtx_y)) return -1;
    else if (vtx_z < Tgt1Lo) return 10;
    else if (vtx_z < TgtDSCut[1]) return 1;
    else if (vtx_z < Tgt2Lo) return 21;
    else if (vtx_z < TgtDSCut[2]) return 2;
    else if (vtx_z < Tgt3Lo) return 32;
    else if (vtx_z < TgtDSCut[3]) return 3;
    else if (vtx_z < TgtWLo) return 63;
    else if (vtx_z < TgtDSCut[6]) return 6;
    else if (vtx_z < Tgt4Lo) return 46;
    else if (vtx_z < TgtDSCut[4]) return 4;
    else if (vtx_z < Tgt5Lo) return 54;
    else if (vtx_z < TgtDSCut[5]) return 5;
    else return 0;
  }

  std::vector<double> XYProjToTgt(int tgt, double vtx_x, double vtx_y, double vtx_z, std::vector<double> muonP){
    std::vector<double> newXY = {-999,-999};
    newXY[0] = vtx_x;
    newXY[1] = vtx_y;
    double zCenter;
    if (tgt == 1) zCenter = PlotUtils::TargetUtils::Get().GetTarget1CenterZMC();
    else if (tgt == 2) zCenter = PlotUtils::TargetUtils::Get().GetTarget2CenterZMC();
    else if (tgt == 3) zCenter = PlotUtils::TargetUtils::Get().GetTarget3CenterZMC();
    //else if (tgt == 36) zCenter = PlotUtils::TargetUtils::Get().GetTarget3CarbonCenterZMC();
    else if (tgt == 4) zCenter = PlotUtils::TargetUtils::Get().GetTarget4CenterZMC();
    else if (tgt == 5) zCenter = PlotUtils::TargetUtils::Get().GetTarget5CenterZMC();
    else if (tgt == 6) zCenter = 0.5*(PlotUtils::TargetProp::WaterTarget::Face+PlotUtils::TargetProp::WaterTarget::Back);
    else return newXY;
    newXY[0] = vtx_x + (muonP[0]/muonP[2])*(zCenter-vtx_z);
    newXY[1] = vtx_y + (muonP[1]/muonP[2])*(zCenter-vtx_z);
    return newXY;
  }

  int GetRecoTarget15(double vtx_x, double vtx_y){
    double u = PlotUtils::TargetUtils::Get().GetCoordU(vtx_x,vtx_y);
    double udist = u - PlotUtils::TargetProp::offset_pb_fe;
    if (fabs(udist) <= 25.0) return 99;
    if (udist < 0.0) return 26;
    return 82;
  }

  int GetRecoTarget2(double vtx_x, double vtx_y){
    double d = PlotUtils::TargetUtils::Get().GetCoordD(vtx_x,vtx_y);
    double ddist = d - PlotUtils::TargetProp::offset_pb_fe;
    if (fabs(ddist) <= 25.0) return 99;
    if (ddist < 0.0) return 26;
    return 82;
  }

  int GetRecoTarget3(double vtx_x, double vtx_y){
    double c = PlotUtils::TargetUtils::Get().GetCoordC(vtx_x,vtx_y);
    double cdist = c;
    if (fabs(cdist) <= 25.0) return 99;
    if (cdist >= 0.0) return 6;
    if (fabs(vtx_x) <= 25.0) return 99;
    if (vtx_x < 0.0) return 26;
    return 82;
  }

  int GetRecoTargetCode(double vtx_x, double vtx_y, double vtx_z, std::vector<double>muonP){
    int TgtByZ = GetRecoTargetZ(vtx_x,vtx_y,vtx_z);
    if (TgtByZ < 0) return -999;
    //Commented out to accommodate change which removed these as options from the TargetCode list so that I can treat the plastic as cut out, and let any actually weird events find their way into the Other.
    /*if (TgtByZ == 10) return 1000;
    if (TgtByZ == 21) return 2100;
    if (TgtByZ == 32) return 3200;
    if (TgtByZ == 63) return 6300;
    if (TgtByZ == 46) return 4600;
    if (TgtByZ == 54) return 5400;
    if (TgtByZ == 0) return 7500;
    */
    if (TgtByZ == 10 || TgtByZ == 21 || TgtByZ == 32 || TgtByZ == 63 || TgtByZ == 46 || TgtByZ == 54 || TgtByZ == 0) return -1;
    if (TgtByZ == 4){
      std::vector<double> newXY = XYProjToTgt(4, vtx_x, vtx_y, vtx_z, muonP);
      if (!PlotUtils::TargetUtils::Get().IsInHexagon(newXY[0],newXY[1])) return -999;
      return 4482;
    }
    if (TgtByZ == 6){
      std::vector<double> newXY = XYProjToTgt(6, vtx_x, vtx_y, vtx_z, muonP);
      if (!PlotUtils::TargetUtils::Get().IsInHexagon(newXY[0],newXY[1])) return -999;
      return 6666;
    }
    if (TgtByZ == 1){
      std::vector<double> newXY = XYProjToTgt(1, vtx_x, vtx_y, vtx_z, muonP);
      if (!PlotUtils::TargetUtils::Get().IsInHexagon(newXY[0],newXY[1])) return -999;
      int mat = GetRecoTarget15(newXY[0],newXY[1]);
      if (mat == 99) return 1199;
      if (mat == 26) return 1126;
      if (mat == 82) return 1182;
    }
    if (TgtByZ == 2){
      std::vector<double> newXY = XYProjToTgt(2, vtx_x, vtx_y, vtx_z, muonP);
      if (!PlotUtils::TargetUtils::Get().IsInHexagon(newXY[0],newXY[1])) return -999;
      int mat = GetRecoTarget2(newXY[0],newXY[1]);
      if (mat == 99) return 2299;
      if (mat == 26) return 2226;
      if (mat == 82) return 2282;
    }
    if (TgtByZ == 3){
      std::vector<double> newXY = XYProjToTgt(3, vtx_x, vtx_y, vtx_z, muonP);
      if (!PlotUtils::TargetUtils::Get().IsInHexagon(newXY[0],newXY[1])) return -999;
      int mat = GetRecoTarget3(newXY[0],newXY[1]);
      if (mat == 99) return 3399;
      if (mat == 6) return 3306;
      if (mat == 26) return 3326;
      if (mat == 82) return 3382;
    }
    if (TgtByZ == 5){
      std::vector<double> newXY = XYProjToTgt(5, vtx_x, vtx_y, vtx_z, muonP);
      if (!PlotUtils::TargetUtils::Get().IsInHexagon(newXY[0],newXY[1])) return -999;
      int mat = GetRecoTarget15(newXY[0],newXY[1]);
      if (mat == 99) return 5599;
      if (mat == 26) return 5526;
      if (mat == 82) return 5582;
    }
    return -999;
  }

  int GetUSTgtCode(int TgtByZ, double vtx_x, double vtx_y, double vtx_z, std::vector<double> muonP, int ReqTgtID){
    if ( TgtByZ < 10 ) return -999;
    std::vector<double> newXY = {-999, -999};
    double zCenter = 0.0;
    if (TgtByZ == 10){
      if (ReqTgtID != 1 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(1,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = PlotUtils::TargetUtils::Get().GetTarget1CenterZMC();
    }
    else if (TgtByZ == 21){
      if (ReqTgtID != 2 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(2,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = PlotUtils::TargetUtils::Get().GetTarget2CenterZMC();
    }
    else if (TgtByZ == 32){ 
      if (ReqTgtID != 3 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(3,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = PlotUtils::TargetUtils::Get().GetTarget3CenterZMC();
    }
    else if (TgtByZ == 63){
      if (ReqTgtID != 6 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(6,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = 0.5*(PlotUtils::TargetProp::WaterTarget::Face+PlotUtils::TargetProp::WaterTarget::Back);
    }
    else if (TgtByZ == 46){ 
      if (ReqTgtID != 4 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(4,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = PlotUtils::TargetUtils::Get().GetTarget4CenterZMC();
    }
    else if (TgtByZ == 54){ 
      if (ReqTgtID != 5 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(5,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = PlotUtils::TargetUtils::Get().GetTarget5CenterZMC();
    }
    return GetRecoTargetCode(newXY[0], newXY[1], zCenter, muonP);
  }

  int GetDSTgtCode(int TgtByZ, double vtx_x, double vtx_y, double vtx_z, std::vector<double> muonP, int ReqTgtID){
    if ( TgtByZ !=0 && TgtByZ < 11 ) return -999;
    std::vector<double> newXY = {-999, -999};
    double zCenter = 0.0;
    if (TgtByZ == 21){
      if (ReqTgtID != 1 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(1,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = PlotUtils::TargetUtils::Get().GetTarget1CenterZMC();
    }
    else if (TgtByZ == 32){
      if (ReqTgtID != 2 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(2,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = PlotUtils::TargetUtils::Get().GetTarget2CenterZMC();
    }
    else if (TgtByZ == 63){
      if (ReqTgtID != 3 && ReqTgtID != -1) return -999; 
      newXY = XYProjToTgt(3,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = PlotUtils::TargetUtils::Get().GetTarget3CenterZMC();
    }
    else if (TgtByZ == 46){
      if (ReqTgtID != 6 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(6,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = 0.5*(PlotUtils::TargetProp::WaterTarget::Face+PlotUtils::TargetProp::WaterTarget::Back);
    }
    else if (TgtByZ == 54){ 
      if (ReqTgtID != 4 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(4,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = PlotUtils::TargetUtils::Get().GetTarget4CenterZMC();
    }
    else if (TgtByZ == 0){ 
      if (ReqTgtID != 5 && ReqTgtID != -1) return -999;
      newXY = XYProjToTgt(5,vtx_x,vtx_y,vtx_z,muonP);
      zCenter = PlotUtils::TargetUtils::Get().GetTarget5CenterZMC();
    }
    return GetRecoTargetCode(newXY[0], newXY[1], zCenter, muonP);
  }

  bool CorrectTargetMaterial(int tgtCode, int tgtType){
    if (tgtCode < 1000) return false;
    else if (tgtCode == 4482){
      if (tgtType == 82) return true;
      else return false;
    }
    else if (tgtCode == 6666){
      if (tgtType == 8) return true;
      else return false;
    }
    else if (tgtCode == 1126){
      if (tgtType == 26) return true;
      else return false;
    }
    else if (tgtCode == 1182){
      if (tgtType == 82) return true;
      else return false;
    }
    else if (tgtCode == 1199){
      if (tgtType == 26 || tgtType == 82) return true;
      else return false;
    }
    else if (tgtCode == 2226){
      if (tgtType == 26) return true;
      else return false;
    }
    else if (tgtCode == 2282){
      if (tgtType == 82) return true;
      else return false;
    }
    else if (tgtCode == 2299){
      if (tgtType == 26 || tgtType == 82) return true;
      else return false;
    }
    else if (tgtCode == 5526){
      if (tgtType == 26) return true;
      else return false;
    }
    else if (tgtCode == 5582){
      if (tgtType == 82) return true;
      else return false;
    }
    else if (tgtCode == 5599){
      if (tgtType == 26 || tgtType == 82) return true;
      else return false;
    }
    else if (tgtCode == 3306){
      if (tgtType == 6) return true;
      else return false;
    }
    else if (tgtCode == 3326){
      if (tgtType == 26) return true;
      else return false;
    }
    else if (tgtCode == 3382){
      if (tgtType == 82) return true;
      else return false;
    }
    else if (tgtCode == 3399){
      if (tgtType == 26 || tgtType == 82 || tgtType == 6) return true;
      else return false;
    }
    else return true;
  }

  double GetNPlanesDSOfTarget(int tgtCode, double vtx_z){
    int tgtID = tgtCode/1000;
    //std::cout << "Calculated Tgt ID: " << tgtID << std::endl;
    double nPlanes = -999.0;
    double upperBound = TgtDSCut[tgtID];
    int tmpN = 1;
    int small = 0;
    int nMax = 11;
    if (tgtID == 4) nMax = 4;
    if (vtx_z < upperBound) return nPlanes;
    upperBound += step_small;
    while (tmpN < nMax){
      if (vtx_z < upperBound) break;
      if (small == 0){
	upperBound += step_big;
	small = 1;
      }
      else{
	upperBound += step_small;
	small = 0;
      }
      tmpN++;
    }
    if (tmpN < nMax) nPlanes = (double)tmpN;
    return nPlanes;
  }

  double GetNPlanesUSOfTarget(int tgtCode, double vtx_z){
    int tgtID = tgtCode/1000;
    //std::cout << "Calculated Tgt ID: " << tgtID << std::endl;
    double nPlanes = -999.0;
    double lowerBound = TgtUSPlaneBack[tgtID];
    int tmpN = 1;
    int small = 1;
    int nMax = 11;
    if (tgtID == 5) nMax = 4;
    if (vtx_z > lowerBound) return nPlanes;
    lowerBound -= step_small;
    while (tmpN < nMax){
      if (vtx_z > lowerBound) break;
      if (small == 0){
	lowerBound -= step_big;
	small = 1;
      }
      else{
	lowerBound -= step_small;
	small = 0;
      }
      tmpN++;
    }
    if (tmpN < nMax) nPlanes = (double)(-1.0*tmpN);
    return nPlanes;
  }
}

#endif //UTIL_GETRECOTGTZ_H
