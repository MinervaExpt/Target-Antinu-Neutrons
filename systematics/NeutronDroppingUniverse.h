//This is meant to hold what it needs to run a universe which sometimes drops neutrons as per the prescription from Miranda Elkin's LE work.

#ifndef NEUTDROPUNIV_H
#define NEUTDROPUNIV_H

#include "event/CVUniverse.h"

#include "TRandom3.h"

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;

class DropGEANTNeutrons: public CVUniverse{
public:
  DropGEANTNeutrons(PlotUtils::ChainWrapper* chw, double prob=0.5, double thresh=10.0): CVUniverse(chw), m_Prob(prob), m_Thresh(thresh)
  {
    m_Rnd = new TRandom3(0);
    if (m_Prob > 1.0) m_Prob = 1.0;
    else if (m_Prob < 0.0) m_Prob = 0.0;
  }
  
  virtual ~DropGEANTNeutrons() = default;
  
  std::string ShortName() const override
  {
    return "DropGEANTNeutrons";
  }
  
  std::string LatexName() const override
  {
    return "Drop GEANT Neutrons";
  }

  //This function isn't called on the truth tree so I'm not worried about needing to do anything special to make this make sense in the truth systematics.
  NeutronCandidates::NeutCands GetLeadNeutCandOnly() override
  {    
    std::vector<NeutronCandidates::NeutCand> cands = {};    
    int nBlobs = GetNNeutBlobs();
    
    if (nBlobs > 0){
      std::vector<double> Es = GetNeutCandEs();
      std::string toolName = GetAnaToolName();
      std::string branchNameParent = "_BlobParentMCPID";
      std::string branchNamePID = "_BlobMCPID";
      int leadNeutEIndex = -999;
      double maxE = -999;
      for (unsigned int idx = 0; idx < Es.size(); ++idx){
	if (Es.at(idx) > maxE){
	  if (Es.at(idx) < m_Thresh){
	      int parentID = GetVecElemInt((toolName+branchNameParent).c_str(), idx);
	      int ID = GetVecElemInt((toolName+branchNamePID).c_str(), idx);
	      //std::cout << "parent: " << parentID << ", self: " << ID << std::endl;
	      if ((parentID==2112 || ID==2112) && !m_Rnd->Binomial(1,m_Prob)) continue;
	  }
	  maxE = Es.at(idx);
	  leadNeutEIndex = idx;
	}
      }

      if (leadNeutEIndex >= 0){
	cands.push_back(GetNeutCand(leadNeutEIndex));
      }

      m_LeadNeutIndex = leadNeutEIndex;
    }

    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  }

  private:
  double m_Prob;
  double m_Thresh;
  TRandom3* m_Rnd;
};

/*
class DropGENIENeutrons: public CVUniverse{
public:
  DropGENIENeutrons(PlotUtils::ChainWrapper* chw): CVUniverse(chw)
  {
  }
  
  virtual ~DropGENIENeutrons() = default;
  
  std::string ShortName() const override
  {
    return "DropGENIENeutrons";
  }
  
  std::string LatexName() const override
  {
    return "Drop GENIE Neutrons";
  }

  //Need to write a truth function? If so, how do I ensure that the neutron removed is removed in true and reco both...

  NeutronCandidates::NeutCands GetLeadNeutCandOnly() override
  {
    std::vector<NeutronCandidates::NeutCand> cands = {};
    int nBlobs = GetNNeutBlobs();
    
    if (nBlobs > 0){
      std::vector<double> Es = GetNeutCandEs();
      int leadNeutEIndex = std::max_element(Es.begin(),Es.end()) - Es.begin();

      //This lives as a test where it shouldn't affect the total number of events.
      cands.push_back(GetNeutCand(leadNeutEIndex));
      
      if (Es.at(leadNeutEIndex) > 100.0){
	m_LeadNeutIndex = leadNeutEIndex;
      }
    }

    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  }
  
};
*/

UniverseMap GetNeutronDroppingUnivs(PlotUtils::ChainWrapper* chw, double prob=0.5)
{
  UniverseMap error_bands;
  //error_bands["DropGENIENeutrons"].push_back(new DropGENIENeutrons(chw, prob));
  error_bands["DropGEANTNeutrons"].push_back(new DropGEANTNeutrons(chw, prob));

  return error_bands;
};

#endif
