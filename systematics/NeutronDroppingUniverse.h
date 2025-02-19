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
	      if ((parentID==2112 || ID==2112) && m_Rnd->Binomial(1,m_Prob)) continue;
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

class DropGENIENeutrons: public CVUniverse{
public:
  DropGENIENeutrons(PlotUtils::ChainWrapper* chw, double prob=0.5, double thresh=50.0): CVUniverse(chw), m_Prob(prob), m_Thresh(thresh)
  {
    m_Rnd = new TRandom3(0);
    if (m_Prob > 1.0) m_Prob = 1.0;
    else if (m_Prob < 0.0) m_Prob = 0.0;
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

  //This only handles the reco side of things. Not trying to get the truth side at this current juncture.
  NeutronCandidates::NeutCands GetLeadNeutCandOnly() override
  {
    std::vector<double> energiesToDrop;
    std::vector<int> FS_PDGs = GetFSPartPDG();
    std::vector<double> FS_Es = GetFSPartE();

    unsigned int nFS = FS_PDGs.size();
    for (unsigned int iFS = 0; iFS < nFS; ++iFS){
      int PDG = FS_PDGs.at(iFS);
      double E = FS_Es.at(iFS);
      if (PDG==2112 && (E-M_n) < m_Thresh && m_Rnd->Binomial(1,m_Prob)) energiesToDrop.push_back(E);
    }

    unsigned int nSkip = energiesToDrop.size();
    
    std::vector<NeutronCandidates::NeutCand> cands = {}; 
    int nBlobs = GetNNeutBlobs();
    
    if (nBlobs > 0){
      std::vector<double> Es = GetNeutCandEs();
      std::string toolName = GetAnaToolName();
      std::string branchNameTopPID = "_BlobTopMCPID";
      std::string branchNameTopE = "_BlobMCTopTrackE";
      int leadNeutEIndex = -999;
      double maxE = -999;
      for (unsigned int idx = 0; idx < Es.size(); ++idx){
	if (Es.at(idx) > maxE){
	  double TopE = GetVecElem((toolName+branchNameTopE).c_str(), idx);
	  int TopPID = GetVecElemInt((toolName+branchNameTopPID).c_str(), idx);

	  //Only check the kinetic energy as one to skip if the energy is below the threshold for possibly being skipped.
	  if(nSkip > 0 && TopPID==2112 && (TopE-M_n) < m_Thresh){
	    bool skip = false;
	    for (unsigned int iEn=0; iEn < nSkip; ++iEn){
	      if (fabs(TopE-energiesToDrop.at(iEn)) < 0.1){
		skip = true;
		break;
	      }
	    }
	    if (skip) continue;
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

UniverseMap GetNeutronDroppingUnivs(PlotUtils::ChainWrapper* chw, double probGENIE=0.5, double probGEANT=0.5)
{
  UniverseMap error_bands;
  error_bands["DropGENIENeutrons"].push_back(new DropGENIENeutrons(chw, probGENIE));
  error_bands["DropGEANTNeutrons"].push_back(new DropGEANTNeutrons(chw, probGEANT));

  return error_bands;
};

#endif
