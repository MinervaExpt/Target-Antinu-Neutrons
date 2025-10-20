#ifndef UNIVSTUDY_H
#define UNIVSTUDY_H

#include "event/CVUniverse.h"

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;

// What things to study
/*
QE normalization (shape?)
FS Neutron energy weighting?
Proton Energy is undershot according to what I understand from Joel (how to upweight appropriately?)
*/

// Weights can be applied with GetWeightRatioToCV()!!!!!!!! This code doesn't see the event though, and perhaps not the index of the neutron candidate in question...
// For the case of the neutron energy distribution from GENIE, have to do something to keep the total probability in principle... Not sure how to weight or not.

class ProtonTrueEnergy: public CVUniverse{
public:
  ProtonTrueEnergy(PlotUtils::ChainWrapper* chw): CVUniverse(chw)
  {
  }
  
  virtual ~ProtonTrueEnergy() = default;
  
  std::string ShortName() const override
  {
    return "ProtonTrueEnergy";
  }
  
  std::string LatexName() const override
  {
    return "Shift Proton Energy from Neutron Interactions";
  }

  double GetWeightRatioToCV() const override
  {
    double weight = 1.0;
    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchPIDName = "_BlobMCPID";
      std::string branchParentEName = "_BlobMCTrackE";
      
      if (GetVecElem((toolName+branchPIDName).c_str(), m_LeadNeutIndex) == 2212){
	double protKE = (GetVecElem((toolName+branchParentEName).c_str(), m_LeadNeutIndex) - M_p);
	if (protKE < 5) weight = 0.75;
	else if (protKE < 50) weight = 1.25;
      }
    }
        
    return weight;
  }
	
};

class UntrackedUpweight: public CVUniverse{
public:
  double nSigma;
  
  UntrackedUpweight(PlotUtils::ChainWrapper* chw, double nSig=1.0): CVUniverse(chw), nSigma(nSig)
  {
  }
  
  virtual ~UntrackedUpweight() = default;
  
  std::string ShortName() const override
  {
    return "UntrackedUpweight";
  }
  
  std::string LatexName() const override
  {
    return "Upweight Trackable Parents of neutron candidates";
  }

  double GetWeightRatioToCV() const override
  {
    double weight = 1.0;
    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchPIDName = "_BlobTopMCPID";
      std::string branchParentEName = "_BlobMCTrackE";

      int pid = GetVecElem((toolName+branchPIDName).c_str(), m_LeadNeutIndex);

      //Just protons for a bit... Curious if the fraction of upweight in the backgrounds accidentally overpowered that in the signal...
      if (pid == 2212){/* || fabs(pid) == 211){*/
	weight = 1.0 + 0.5 * nSigma;
      }
    }
        
    return weight;
  }
	
};

class MuonUpweight: public CVUniverse{
public:
  double nSigma;
  
  MuonUpweight(PlotUtils::ChainWrapper* chw, double nSig=1.0): CVUniverse(chw), nSigma(nSig)
  {
  }
  
  virtual ~MuonUpweight() = default;
  
  std::string ShortName() const override
  {
    return "MuonUpweight";
  }
  
  std::string LatexName() const override
  {
    return "Upweight Muon-Induced Neutron Candidates";
  }

  double GetWeightRatioToCV() const override
  {
    double weight = 1.0;
    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchPIDName = "_BlobTopMCPID";

      int pid = GetVecElem((toolName+branchPIDName).c_str(), m_LeadNeutIndex);

      //Just protons for a bit... Curious if the fraction of upweight in the backgrounds accidentally overpowered that in the signal...
      if (fabs(pid) == 13){/* || fabs(pid) == 211){*/
	weight = 1.0 + 1.0 * nSigma;
      }
    }
        
    return weight;
  }
	
};

class Neutron3DDiff: public CVUniverse{
public:
  double prob;
  TRandom3* fRandom;
  
  Neutron3DDiff(PlotUtils::ChainWrapper* chw, double p=0.5): CVUniverse(chw), prob(p)
  {
    fRandom = new TRandom3(0);
  }
  
  virtual ~Neutron3DDiff() = default;
  
  std::string ShortName() const override
  {
    return "Neutron3DDiff";
  }
  
  std::string LatexName() const override
  {
    return "Weight neutron candidate is 3D for neutron-induced candidates.";
  }

  NeutronCandidates::NeutCand GetNeutCand(int index) override
  {
    std::vector<double> vtx = GetVtx();
    TVector3 EvtVtx;
    EvtVtx.SetXYZ(vtx.at(0),vtx.at(1),vtx.at(2));
    NeutronCandidates::intCandData intData;
    NeutronCandidates::doubleCandData doubleData;
    std::string toolName = GetAnaToolName();
    for (const auto& intMember: NeutronCandidates::GetBranchIntMap()){
      intData[intMember.first]={};
      for (const auto& branchName: intMember.second){
        intData[intMember.first].push_back(GetVecElemInt((toolName+branchName).c_str(),index));
      }
    }
    for (const auto& doubleMember: NeutronCandidates::GetBranchDoubleMap()){
      doubleData[doubleMember.first]={};
      for (const auto& branchName: doubleMember.second){
        doubleData[doubleMember.first].push_back(GetVecElem((toolName+branchName).c_str(),index));
      }
    }

    NeutronCandidates::NeutCand ret = NeutronCandidates::NeutCand(intData,doubleData,EvtVtx);
    
    if (ret.GetMCPID() == 2112 && ret.GetIs3D()==1 && fRandom->Binomial(1,prob)){
      std::vector<int> is3D;
      is3D.push_back(0);
      ret.SetIs3D(is3D);
    }
    
    return ret;
  }
  
};

class NeutronFSEnergy: public CVUniverse{
public:
  NeutronFSEnergy(PlotUtils::ChainWrapper* chw): CVUniverse(chw)
  {
  }
  
  virtual ~NeutronFSEnergy() = default;
  
  std::string ShortName() const override
  {
    return "NeutronFSEnergy";
  }
  
  std::string LatexName() const override
  {
    return "Shift FS Neutron Energy";
  }

  double GetWeightRatioToCV() const override
  {
    double weight = 1.0;
    double maxNeutKE = GetMaxFSNeutronKE();
    if (maxNeutKE < 30) weight = 0.8;
    else if (maxNeutKE < 60) weight = 1.2;
    else weight = 1.0;
    
    return weight;
  }
	
};

class ProtonEnergy: public CVUniverse{
public:
  ProtonEnergy(PlotUtils::ChainWrapper* chw): CVUniverse(chw)
  {
  }
  
  virtual ~ProtonEnergy() = default;
  
  std::string ShortName() const override
  {
    return "ProtonEnergy";
  }
  
  std::string LatexName() const override
  {
    return "Shift Proton Reco Energy from Neutron Interactions";
  }

  double GetLeadNeutCandE() const override
  {
    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchEName = "_BlobTotalE";
      std::string branchPIDName = "_BlobMCPID";

      double factor = 1.0;

      if (GetVecElem((toolName+branchPIDName).c_str(), m_LeadNeutIndex) == 2212) factor = 1.0-0.035;
      
      return factor*GetVecElem((toolName+branchEName).c_str(), m_LeadNeutIndex);
    }

    return -999;
  }
  
};

class QENorm: public CVUniverse{
public:
  QENorm(PlotUtils::ChainWrapper* chw): CVUniverse(chw)
  {
  }

  virtual ~QENorm() = default;

  std::string ShortName() const override
  {
    return "QENorm";
  }

  std::string LatexName() const override
  {
    return "QE normalization shape effect";
  }

  std::vector<double> SplinePtsY = {1, 1, 1, 0.993, 0.987, 0.987, 0.987, 0.987, 0.987, 0.987, 1, 1.013, 1.04, 1.06, 1.09, 1.13, 1.18, 1.23, 1.29, 1.36, 1.44, 1.52, 1.61, 1.71, 1.83, 1.97, 2.10, 2.24, 2.36, 2.50};
  std::vector<double> SplinePtsX = {0.0, 0.0185, 0.0244, 0.0313, 0.0395, 0.0507, 0.0650, 0.0809, 0.102, 0.129, 0.171, 0.222, 0.293, 0.381, 0.496, 0.627, 0.803, 1.00, 1.21, 1.46, 1.72, 1.99, 2.34, 2.75, 3.09, 3.47, 3.79, 4.13, 4.51, 5};

  int FindMinPtIndex(double Q2) const
  {
    int point = -1;
    for (int iPt=0; iPt < SplinePtsX.size(); ++iPt){
      double pt = SplinePtsX.at(iPt);
      if (pt <= Q2){
	point = iPt;
      }
      else break;
    }
    return point;
  }
  
  double RoughWeightingToTejin(double Q2) const
  {
    double weight = 1.0;
    int idx = FindMinPtIndex(Q2);
    if (idx >= (SplinePtsX.size() - 1)){
      weight = SplinePtsY.at(SplinePtsX.size()-1);
    }
    else if (idx >= 0){
      double loX = SplinePtsX.at(idx);
      double loY = SplinePtsY.at(idx);
      double hiX = SplinePtsX.at(idx+1);
      double hiY = SplinePtsY.at(idx+1);
      weight = loY + ((Q2-loX)/(hiX-loX))*(hiY-loY); //Linear extrapolation between spline pts
    }
    
    return weight;
  }
    
  double GetWeightRatioToCV() const override
  {
    double weight = 1.0;
    double trueQ2 = GetQ2True()*MeVGeV*MeVGeV;
    weight = RoughWeightingToTejin(trueQ2);
    return weight;
  }

};

UniverseMap GetStudyUnivs(PlotUtils::ChainWrapper* chw)
{
  UniverseMap error_bands;
  //error_bands["ProtonTrueEnergy"].push_back(new ProtonTrueEnergy(chw));
  //error_bands["MuonUpweight"].push_back(new MuonUpweight(chw,1.0));
  //error_bands["MuonUpweight"].push_back(new MuonUpweight(chw,-1.0));
  //error_bands["AllowHighAngleTracks"].push_back(new AllowHighAngleTracks(chw));
  //error_bands["ProtonEnergy"].push_back(new ProtonEnergy(chw));
  //error_bands["NeutronFSEnergy"].push_back(new NeutronFSEnergy(chw));
  //error_bands["QENorm"].push_back(new QENorm(chw));

  error_bands["MuonUpweight"].push_back(new MuonUpweight(chw,1.0));
  error_bands["UntrackedUpweight"].push_back(new UntrackedUpweight(chw,1.0));
  error_bands["Neutron3DDiff"].push_back(new Neutron3DDiff(chw));
  return error_bands;
};

#endif
