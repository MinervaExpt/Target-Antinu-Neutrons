//This is meant to hold what it needs to run a universe which sometimes drops neutrons as per the prescription from Miranda Elkin's LE work.

#ifndef FSIREPLACEUNIV_H
#define FSIREPLACEUNIV_H

#include "event/CVUniverse.h"

//Something to make the functions available for the

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;

class FSIReplace: public CVUniverse{
public:
  FSIReplace(PlotUtils::ChainWrapper* chw): CVUniverse(chw)
  {
    m_indices = {};//Initializing the indices of the final state particles
  }
  
  virtual ~FSIReplace() = default;
  
  std::string ShortName() const override
  {
    return "FSIReplace";
  }
  
  std::string LatexName() const override
  {
    return "Replace BAD FSI";
  }

  //Borrowed directly from weight_fsi
  double GetGenieBEinMeV(int A) const{
    // From UserPhysicsOptions.xml file
    // these are hard coded, so we can get an exact match
    // but beware if you try to port this code to any new GENIE.
    
    if (A == 1) return 0; // hydrogen
    if (A == 6) return 17.0; // lithium
    if (A == 12) return 25.0; // carbon
    if (A == 16) return 27.0; // oxygen
    if (A == 24) return 32.0; // magnesium
    if (A == 40) return 29.5; // argon
    if (A == 48) return 30.0; // Ti48
    if (A == 56) return 36.0; // 56 iron
    if (A == 58) return 36.0; // 58 nickel
    if (A >= 206 && A <= 208) return 44.0; // 208 lead
    
    // Specialty in MINERvA,picked off numbers by hand.
    if (A == 28) return 8.219751;  // silicon
    if (A == 27) return 8.115287;  // aluminum
    if (A == 14) return 7.185166;  // nitrogen
    if (A == 55) return 8.653063;  // iron55 or manganese55
    if (A == 35) return 8.347164;  // chlorine  
    
    if (A == 4) return 5.0;  // this is rough, 
    
    //else
    // this is a problem, all other numbers come from a semi-empirical binding energy formula
    // need to back off the QE precision for those.
    return 8.0;   // GENIE defaults to something like this.  Needs to be exact.  Check it.
  }

  bool DoesParentage14(int mother, std::vector<int> mothers, std::vector<int> statuses) const{
    int newMother = mother;
    while (newMother >= 0){
      if (statuses.at(newMother) == 14) return true;
      else if (newMother == mothers.at(newMother)) break;//Avoids recursive loop that I didn't expect to be possible... scattering off of electron seems to be the issue...
      else newMother = mothers.at(newMother);
    }
    return false;
  }

  std::vector<int> RecursedIndices(int first, int last, std::vector<int> statuses, std::vector<int> firsts, std::vector<int> lasts) const{
    std::vector<int> ret;
    for (int idx=first; idx<=last; ++idx){
      if (statuses.at(idx) == 1) ret.push_back(idx);
      else{
	std::vector<int> tmp = RecursedIndices(firsts.at(idx), lasts.at(idx), statuses, firsts, lasts);
	for (auto index: tmp) ret.push_back(index);
      }
    }
    return ret;
  }
  
  //Copies the logic of weight_fsi to decide if elastic FSI called. Saves the preFSI index if elastic, otherwise saves what should be all the post FSI stuff.
  std::vector<int> FillIndices() const{
    std::vector<int> indices;
    std::vector<int> statuses = GetVec<int>("mc_er_status");
    std::vector<int> PDGs = GetVec<int>("mc_er_ID");
    std::vector<int> mothers = GetVec<int>("mc_er_mother");
    std::vector<int> firsts = GetVec<int>("mc_er_FD");
    std::vector<int> lasts = GetVec<int>("mc_er_LD");
    std::vector<double> Es = GetVec<double>("mc_er_E");
    for (int idx=0; idx < statuses.size(); ++idx){
      int status = statuses.at(idx);
      int ID = PDGs.at(idx);
      int fd = firsts.at(idx);
      int ld = lasts.at(idx);
      int mother = mothers.at(idx);
      double E = Es.at(idx);
      //std::cout << "Particle with index: " << idx << ", PDG: " << ID << ", Status: " << status << ", energy: " << E << ", Mother: " << mother << ", First Daughter Index: " << fd << ", Last Daughter Index: " << ld << std::endl;
	
      if (status == 14){
	//int fd = firsts.at(idx);
	//int ld = lasts.at(idx);
	//int ID = PDGs.at(idx);
	bool keepPostFSI = true;
	
	if (fd==ld && PDGs.at(fd)==ID && (ID==2112 || ID==2212 || ID==111 || fabs(ID)==211)){
	  double offset = GetGenieBEinMeV(GetTargetA());
	  double tolerance = 0.0001;
	  if (fabs(offset - 8.0) < 0.1) tolerance = 1.2;
	  if (offset < 6.0) tolerance = 1.2;
	  if (GetInteractionType() != 1) offset = 0.0;
	  if (fabs(Es.at(idx) - Es.at(fd) - offset) >= tolerance) keepPostFSI = false; //Replace the postFSI with the preFSI for elastics
	}
	
	if (keepPostFSI){
	  std::vector<int> recursedIndices = RecursedIndices(fd, ld, statuses, firsts, lasts);
	  for (auto index: recursedIndices) indices.push_back(index);
	}
	else{
	  indices.push_back(idx);
	}
      }
      else if (status == 1 && !DoesParentage14(idx, mothers, statuses)){
	indices.push_back(idx);
      }
    }
    return indices;
  }
  
  void OnNewEntry() override{
    m_LeadNeutIndex = -999;//Resetting to avoid any possible mishaps with the indexing of an array.
    m_indices.clear();
    m_indices = FillIndices();
  }
  
  std::vector<int> GetFSPartPDG() const override {
    std::vector<int> indices;
    if (m_indices.size() <= 0){
      indices = FillIndices();
    }
    else indices = m_indices;
    std::vector<int> FSPDGs;
    std::vector<int> MCPDGs = GetVec<int>("mc_er_ID");

    for (auto idx:indices){
      FSPDGs.push_back(MCPDGs.at(idx));
    }
    
    return FSPDGs;
  }

  std::vector<double> GetFSPartE() const override {
    std::vector<int> indices;
    if (m_indices.size() <= 0){
      indices = FillIndices();
    }
    else indices = m_indices;
    std::vector<double> FSEs;
    std::vector<double> MCEs = GetVec<double>("mc_er_E");

    for (auto idx:indices){
      FSEs.push_back(MCEs.at(idx));
    }
    
    return FSEs;
  }

  std::vector<int> m_indices;
};

UniverseMap GetFSIReplaceUnivs(PlotUtils::ChainWrapper* chw)
{
  UniverseMap error_bands;
  error_bands["FSIReplace"].push_back(new FSIReplace(chw));

  return error_bands;
};

#endif
