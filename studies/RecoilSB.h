//studies includes
#include "studies/Study.h"

//David's includes
#include "event/NeutCands.h"
#include "util/NeutronVariable.h"
#include "util/EventVariable.h"
#include "util/Variable.h"
#include "util/Variable2D.h"
#include "util/GetBackgroundID.h"
#include "event/CVUniverse.h"

//Renamed to PreRecoil
class PreRecoil: public Study
{
  private:
    std::vector<Variable*> fVars;
    std::vector<Variable2D*> fVars2D;
    std::map<int,Variable*> fRecoilBinned;
    bool fSplitRecoil;

  public:
    PreRecoil(std::vector<Variable*> vars,
		     std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		     std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
	     std::vector<CVUniverse*>& data_error_bands, bool splitRecoil): Study(), fSplitRecoil(splitRecoil)
    {

      if (fSplitRecoil){
	//Constructing recoil variable in bins of pt/pz. Just make a map of it all.
	std::vector<double> myRecoilBins;
	const double myRecoilBinWidth = 1.0/50.;
	for(int whichBin = 0; whichBin < 51; ++whichBin) myRecoilBins.push_back(myRecoilBinWidth * whichBin);

	int nBins = 14;//Hard-coded from bins from Amit's analysis... Not appropriate binning given my cuts, but going to get the structure going first.
	int lowBin = -1;//Default return value. hard-coded as well.

	for (int iBin=lowBin; iBin < nBins; ++iBin){
	  std::string binName = std::to_string(iBin);
	  if (iBin == lowBin) binName = "lost";
	  fRecoilBinned[iBin] = new Variable(("recoilE_bin_"+binName).c_str(), "Recoil E [GeV]", myRecoilBins, &CVUniverse::GetDANRecoilEnergyGeV);
	  fRecoilBinned[iBin]->InitializeMCHists(mc_error_bands, truth_error_bands);
	  fRecoilBinned[iBin]->InitializeDATAHists(data_error_bands);
	}
      }

      for (auto& var : vars){
	fVars.push_back(new Variable((var->GetName()+"_PreRecoilCut").c_str(), var->GetAxisLabel(), var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
      }

      fVars2D = {
	new Variable2D(*fVars[4],*fVars[3]),//recoil v. Q2     
	new Variable2D(*fVars[fVars.size()-1],*fVars[0]),//pT v. recoilQ2Bin     
      };

      for(auto& var: fVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars) var->InitializeDATAHists(data_error_bands);

      for(auto& var: fVars2D) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars2D) var->InitializeDATAHists(data_error_bands);
    }

    void SaveOrDraw(TDirectory& outDir){
      return;
    }

    void SaveOrDrawMC(TFile& outFile)
    {
      for (auto& recoil : fRecoilBinned) recoil.second->WriteMC(outFile);
      for (auto& var : fVars) var->WriteMC(outFile);
      for (auto& var : fVars2D) var->Write(outFile);
    }

    void SaveOrDrawData(TFile& outFile)
    {
      for (auto& recoil : fRecoilBinned) recoil.second->WriteData(outFile);
      for (auto& var : fVars) var->WriteData(outFile);
    }
    
  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      std::bitset<64> SBStat = evt.GetSideBandStat();
      int intType = evt.GetIntType();
      int tgtType = evt.GetTgtZ();
      int leadBlobType = evt.GetLeadingNeutCand().GetPDGBin();
      int iBin = evt.GetBinPTPZ();
      
      if (evt.IsMC()){
	if (evt.IsSignal()){
	  if (fSplitRecoil){
	    fRecoilBinned[iBin]->selectedMCReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight); //"Fake data" for closure
	    fRecoilBinned[iBin]->selectedSignalReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_SigIntTypeHists)[intType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	  }

	  for (auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }

	  for (auto& var: fVars2D){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	  }
	}
	
	else{
	  int bkgd_ID = -1;	  
	  bkgd_ID = util::GetBackgroundID(univ);
	  
	  if (fSplitRecoil){
	    fRecoilBinned[iBin]->selectedMCReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight); //"Fake data" for closure
	    (*fRecoilBinned[iBin]->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_BkgIntTypeHists)[intType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	  }

	  for(auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }

	  for(auto& var: fVars2D){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	  }
	}
      }
      
      else{
	if (fSplitRecoil) fRecoilBinned[iBin]->dataHist->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), 1);

	for (auto& var : fVars){
	  var->dataHist->FillUniverse(&univ, var->GetRecoValue(univ), 1);
	}

	for (auto& var : fVars2D){
	  var->dataHist->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), 1);
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

//Actually requires the sideband to be satisfied. Written as though the recoil cut is the 0th sideband cut
class RecoilSB: public Study
{
  private:
    std::vector<Variable*> fVars;
    std::vector<Variable2D*> fVars2D;
    std::map<int,Variable*> fRecoilBinned;
    bool fSplitRecoil;

  public:
    RecoilSB(std::vector<Variable*> vars,
		     std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		     std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
	     std::vector<CVUniverse*>& data_error_bands, bool splitRecoil): Study(), fSplitRecoil(splitRecoil)
    {

      if (fSplitRecoil){
	//Constructing recoil variable in bins of pt/pz. Just make a map of it all.
	std::vector<double> myRecoilBins;
	const double myRecoilBinWidth = 1.0/50.;
	for(int whichBin = 0; whichBin < 51; ++whichBin) myRecoilBins.push_back(myRecoilBinWidth * whichBin);

	int nBins = 14;//Hard-coded from bins from Amit's analysis... Not appropriate binning given my cuts, but going to get the structure going first.
	int lowBin = -1;//Default return value. hard-coded as well.

	for (int iBin=lowBin; iBin < nBins; ++iBin){
	  std::string binName = std::to_string(iBin);
	  if (iBin == lowBin) binName = "lost";
	  fRecoilBinned[iBin] = new Variable(("recoilE_RecoilSB_bin_"+binName).c_str(), "Recoil E [GeV]", myRecoilBins, &CVUniverse::GetDANRecoilEnergyGeV);
	  fRecoilBinned[iBin]->InitializeMCHists(mc_error_bands, truth_error_bands);
	  fRecoilBinned[iBin]->InitializeDATAHists(data_error_bands);
	}
      }

      for (auto& var : vars){
	fVars.push_back(new Variable((var->GetName()+"_RecoilSB").c_str(), var->GetAxisLabel(), var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
      }

      fVars2D = {
	new Variable2D(*fVars[4],*fVars[3]),//recoil v. Q2     
	new Variable2D(*fVars[fVars.size()-1],*fVars[0]),//pT v. recoilQ2Bin     
      };

      for(auto& var: fVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars) var->InitializeDATAHists(data_error_bands);

      for(auto& var: fVars2D) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars2D) var->InitializeDATAHists(data_error_bands);
    }

    void SaveOrDraw(TDirectory& outDir){
      return;
    }

    void SaveOrDrawMC(TFile& outFile)
    {
      for (auto& recoil : fRecoilBinned) recoil.second->WriteMC(outFile);
      for (auto& var : fVars) var->WriteMC(outFile);
      for (auto& var : fVars2D) var->Write(outFile);
    }

    void SaveOrDrawData(TFile& outFile)
    {
      for (auto& recoil : fRecoilBinned) recoil.second->WriteData(outFile);
      for (auto& var : fVars) var->WriteData(outFile);
    }
    
  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      std::bitset<64> SBStat = evt.GetSideBandStat();
      int intType = evt.GetIntType();
      int tgtType = evt.GetTgtZ();
      int leadBlobType = evt.GetLeadingNeutCand().GetPDGBin();
      int iBin = evt.GetBinPTPZ();

      if (SBStat[0] == 1) return;
      
      if (evt.IsMC()){
	if (evt.IsSignal()){
	  if (fSplitRecoil){
	    fRecoilBinned[iBin]->selectedMCReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight); //"Fake data" for closure
	    fRecoilBinned[iBin]->selectedSignalReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_SigIntTypeHists)[intType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	  }

	  for (auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }

	  for (auto& var: fVars2D){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	  }
	}
	
	else{
	  int bkgd_ID = -1;	  
	  bkgd_ID = util::GetBackgroundID(univ);
	  
	  if (fSplitRecoil){
	    fRecoilBinned[iBin]->selectedMCReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight); //"Fake data" for closure
	    (*fRecoilBinned[iBin]->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_BkgIntTypeHists)[intType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	  }

	  for(auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }

	  for(auto& var: fVars2D){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	  }
	}
      }
      
      else{
	if (fSplitRecoil) fRecoilBinned[iBin]->dataHist->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), 1);

	for (auto& var : fVars){
	  var->dataHist->FillUniverse(&univ, var->GetRecoValue(univ), 1);
	}

	for (auto& var : fVars2D){
	  var->dataHist->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), 1);
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
