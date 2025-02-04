//This is meant to hold what it needs to run a universe which sometimes drops neutrons as per the prescription from Miranda Elkin's LE work.

#include "event/CVUniverse.h"

#include "TRandom.h"

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;

/*
class DropGEANTNeutrons: public CVUniverse{
public:
DropGEANTNeutrons(/*Whatever is needed to configure the correct handling,*//*PlotUtils::ChainWrapper* chw, double prob=0.5, double thresh=10.0): CVUniverse(chw), fProb(prob), fThresh(thresh), fLeadNeutIndex(-1354), fRnd(new TRandom())
  {
    if (fProb > 1.0) fProb = 1.0;
    else if (fProb < 0.0) fProb = 0.0;
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
  
  NeutronCandidates::NeutCands GetLeadNeutCandOnly() override
  {
    TString toolName = GetAnaToolName();
    
    std::vector<NeutronCandidates::NeutCand> cands = {};    
    int nBlobs = GetNNeutBlobs();
    
    if (nBlobs > 0){
      std::vector<double> Es = GetNeutCandEs();
      std::vector<double> parentIDs = GetVec<double>(toolName+"_BlobParentMCPID");
      int leadNeutEIndex = -999;
      double maxE = -999;
      for (unsigned int idx = 0; idx < Es.size(); ++idx){
	if (Es.at(idx) > maxE){
	  if ((Es.at(idx) < fThresh) && !fRnd->Binomial(1,fProb)) continue;
	  maxE = Es.at(idx);
	  leadNeutEIndex = idx;
	}
      }

      if (leadNeutEIndex >= 0){
	cands.push_back(GetNeutCand(leadNeutEIndex));
      }

      fLeadNeutIndex = leadNeutEIndex;
    }

    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  }

private:
  double fProb;
  double fThresh;
  TRandom* fRnd;
};
*/

class DropGENIENeutrons: public CVUniverse{
public:
  DropGENIENeutrons(/*Whatever is needed to configure the correct handling,*/PlotUtils::ChainWrapper* chw): CVUniverse(chw)
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

      if (Es.at(leadNeutEIndex) > 100.0){
	cands.push_back(GetNeutCand(leadNeutEIndex));
      }
    }

    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  }
  
};

UniverseMap GetNeutronDroppingUnivs(PlotUtils::ChainWrapper* chw/*,prob*/)
{
  UniverseMap error_bands;
  error_bands["DropNeutrons"].push_back(new DropGENIENeutrons(chw/*,prob*/));
  //error_bands["DropGEANTNeutrons"].push_back(new DropGEANTNeutrons(chw, prob));

  return error_bands;
};
