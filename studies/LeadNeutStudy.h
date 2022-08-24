//studies includes
#include "studies/Study.h"

//David's includes
#include "event/NeutCands.h"
#include "util/NeutronVariable.h"
#include "util/Variable.h"
#include "util/GetBackgroundID.h"
#include "event/CVUniverse.h"

//TODO: ONLY FILL THESE VARIABLES WHEN THE EVENT IS SELECTED?

class LeadNeutStudy: public Study
{
  private:
    Variable* fMatchesLead;

  public:
    LeadNeutStudy(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		  std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
		  std::vector<CVUniverse*>& data_error_bands): Study()
    {
      std::vector<double> myBins = {-1.5,-0.5,0.5,1.5};

      fMatchesLead = new Variable("LeadCandMatches");

      fMatchesLead->InitializeMCHists(mc_error_bands, truth_error_bands);
      fMatchesLead->InitializeDATAHists(data_error_bands);
    }

    void SaveOrDraw(TDirectory& outDir){
      return;
    }

    void SaveOrDrawMC(TFile& outFile)
    {
      fMatchesLead->WriteMC(outFile);
    }

    void SaveOrDrawData(TFile& outFile)
    {
      fMatchesLead->WriteData(outFile);
    }

  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      NeutronCandidates::NeutCand leadCand = evt.GetLeadingNeutCand();
      int intType = evt.GetIntType();
      int tgtType = evt.GetTgtZ();
      int leadBlobType = leadCand.GetPDGBin();
      TLorentzVector maxFSNeutMom = univ.GetMaxFSNeutronMom();
      TLorentzVector leadFSMom = leadCand.GetTopMCMom();
      int matchValue = leadCand.MatchesFSNeutron(maxFSNeutMom,0.01);
      std::cout << "DEBUGGING" << std::endl;
      std::cout << "Lead Blob Type: " << leadBlobType << std::endl;
      std::cout << "Lead Blob Type: " << leadBlobType << std::endl;
      std::cout << "Max FS Neutron 4-Momentum:" << std::endl;
      std::cout << "Px: " << maxFSNeutMom.X() << ", Py: " << maxFSNeutMom.Y() << ", Pz: " << maxFSNeutMom.Z() << ", E: " << maxFSNeutMom.T() << std::endl;
      std::cout << "Lead Blob Type FS Particle 4-Momentum:" << std::endl;
      std::cout << "Px: " << leadFSMom.X() << ", Py: " << leadFSMom.Y() << ", Pz: " << leadFSMom.Z() << ", E: " << leadFSMom.T() << std::endl;
      std::cout << "Reported Match Value: " << matchValue << std::endl;
      std::cout << std::endl;

      if (evt.IsMC()){
	if (evt.IsSignal()){
	  fMatchesLead->selectedMCReco->FillUniverse(&univ, matchValue, weight);
	  fMatchesLead->dataHist->FillUniverse(&univ, matchValue, weight);
	  fMatchesLead->selectedSignalReco->FillUniverse(&univ, matchValue, weight);
	  (*fMatchesLead->m_SigIntTypeHists)[intType].FillUniverse(&univ, matchValue, weight);
	  (*fMatchesLead->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, matchValue, weight);
	  (*fMatchesLead->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, matchValue, weight);
	}
	
	else {
	  int bkgd_ID = -1;
	  bkgd_ID = util::GetBackgroundID(univ);

	  fMatchesLead->selectedMCReco->FillUniverse(&univ, matchValue, weight);
	  fMatchesLead->dataHist->FillUniverse(&univ, matchValue, weight);
	  (*fMatchesLead->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, matchValue, weight);
	  (*fMatchesLead->m_BkgIntTypeHists)[intType].FillUniverse(&univ, matchValue, weight);
	  (*fMatchesLead->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, matchValue, weight);
	  (*fMatchesLead->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, matchValue, weight);
	}
      }
      
      return;
    }
    
    //All of your plots happen here so far.
    void fillSelectedSignal(const CVUniverse& univ, const NeutronEvent& evt, const double weight)
    {
      return;
      //(*m_VarToGENIELabel)[univ.GetInteractionType()].FillUniverse(&univ, fReco(univ, evt), weight);
    }

    //Do nothing for now...  Good place for efficiency denominators in the future.
    void fillTruthSignal(const CVUniverse& univ, const NeutronEvent& evt, const double weight) { return; }
};
