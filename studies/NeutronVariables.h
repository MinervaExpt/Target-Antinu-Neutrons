//studies includes
#include "studies/Study.h"

//David's includes
#include "event/NeutCands.h"
#include "util/NeutronVariable.h"
#include "util/Variable.h"
#include "util/GetBackgroundID.h"
#include "event/CVUniverse.h"
#include "TString.h"

//TODO: ONLY FILL THESE VARIABLES WHEN THE EVENT IS SELECTED?

class NeutronVariables: public Study
{
  private:
    std::vector<NeutronVariable*> fLeadVars;
    std::vector<NeutronVariable*> fAllVars;
    std::vector<NeutronVariable*> fTgtVars;
    std::vector<NeutronVariable*> fTrackVars;
    double fBound;
    double fMinZ;
    bool fAllCuts;

  public:
    NeutronVariables(double tgtBoundary,
		     double minZ,
		     bool allCuts,
		     std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		     std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
		     std::vector<CVUniverse*>& data_error_bands): Study(), fBound(tgtBoundary), fMinZ(minZ), fAllCuts(allCuts)
    {
      std::vector<double> myBlobEBins;
      const double myBlobEBinWidth = 3.;
      for (int whichBin=0; whichBin < 51; ++whichBin) myBlobEBins.push_back(myBlobEBinWidth * whichBin);

      std::vector<double> mySelfAngleBins;
      const double mySelfAngleBinWidth = 3.2/180.;
      for (int whichBin=0; whichBin < 181; ++whichBin) mySelfAngleBins.push_back(mySelfAngleBinWidth * whichBin);

      std::vector<double> myPDGBins;
      const double myPDGBinWidth = 1.;
      for (int whichBin=0; whichBin < 11; ++whichBin) myPDGBins.push_back(myPDGBinWidth * whichBin);

      std::vector<double> myLenBins;
      const double myLenBinWidth = 10.;
      for (int whichBin=0; whichBin < 51; ++whichBin) myLenBins.push_back(myLenBinWidth * whichBin);

      std::vector<double> myDEDXBins;
      const double myDEDXBinWidth = 2.;
      for (int whichBin=0; whichBin < 26; ++whichBin) myDEDXBins.push_back(myDEDXBinWidth * whichBin);

      std::vector<double> myVtxDistBins;
      const double myVtxDistBinWidth = 10.;
      for (int whichBin=0; whichBin < 101; ++whichBin) myVtxDistBins.push_back(myVtxDistBinWidth * whichBin);

      std::vector<double> myZPosBins;
      const double myZPosBinWidth = 10.;
      for (int whichBin=0; whichBin < ((9300-fMinZ)/10)+1; ++whichBin) myZPosBins.push_back(myZPosBinWidth * whichBin + fMinZ);

      std::string nameTag = (fAllCuts) ? "" : "_PreRecoilCut" ;

      fLeadVars = {
	new NeutronVariable(("leadBlob"+nameTag+"_blobE").c_str(),"E [MeV]", myBlobEBins,&NeutronCandidates::NeutCand::GetTotalE),
	new NeutronVariable(("leadBlob"+nameTag+"_SelfAngle").c_str(),"Angle [rad]", mySelfAngleBins,&NeutronCandidates::NeutCand::GetAngleToFP),
	new NeutronVariable(("leadBlob"+nameTag+"_primary_parent").c_str(),"", myPDGBins,&NeutronCandidates::NeutCand::GetPDGBin),
	new NeutronVariable(("leadBlob"+nameTag+"_length").c_str(),"len. [mm]", myLenBins,&NeutronCandidates::NeutCand::GetLength),
	new NeutronVariable(("leadBlob"+nameTag+"_avg_dEdx").c_str(),"dE/dx [MeV/mm]", myDEDXBins,&NeutronCandidates::NeutCand::GetDEDX),
	new NeutronVariable(("leadBlob"+nameTag+"_dist").c_str(),"dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxDist),
	new NeutronVariable(("leadBlob"+nameTag+"_Zdist").c_str(),"Z dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxZDist),
	new NeutronVariable(("leadBlob"+nameTag+"_ZPos").c_str(),"Z [mm]", myZPosBins,&NeutronCandidates::NeutCand::GetZPos),
      };
      fAllVars = {
	new NeutronVariable(("All"+nameTag+"_blobE").c_str(),"E [MeV]", myBlobEBins,&NeutronCandidates::NeutCand::GetTotalE),
	new NeutronVariable(("All"+nameTag+"_primary_parent").c_str(),"", myPDGBins,&NeutronCandidates::NeutCand::GetPDGBin),
	new NeutronVariable(("All"+nameTag+"_length").c_str(),"len. [mm]", myLenBins,&NeutronCandidates::NeutCand::GetLength),
	new NeutronVariable(("All"+nameTag+"_avg_dEdx").c_str(),"dE/dx [MeV/mm]", myDEDXBins,&NeutronCandidates::NeutCand::GetDEDX),
	new NeutronVariable(("All"+nameTag+"_dist").c_str(),"dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxDist),
	new NeutronVariable(("All"+nameTag+"_Zdist").c_str(),"Z dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxZDist),
      };
      fTgtVars = {
	new NeutronVariable(("target"+nameTag+"_blobE").c_str(),"E [MeV]", myBlobEBins,&NeutronCandidates::NeutCand::GetTotalE),
	new NeutronVariable(("target"+nameTag+"_primary_parent").c_str(),"", myPDGBins,&NeutronCandidates::NeutCand::GetPDGBin),
	new NeutronVariable(("target"+nameTag+"_length").c_str(),"len. [mm]", myLenBins,&NeutronCandidates::NeutCand::GetLength),
	new NeutronVariable(("target"+nameTag+"_avg_dEdx").c_str(),"dE/dx [MeV/mm]", myDEDXBins,&NeutronCandidates::NeutCand::GetDEDX),
	new NeutronVariable(("target"+nameTag+"_dist").c_str(),"dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxDist),
	new NeutronVariable(("target"+nameTag+"_Zdist").c_str(),"Z dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxZDist),
      };
      fTrackVars = {
	new NeutronVariable(("tracker"+nameTag+"_blobE").c_str(),"E [MeV]", myBlobEBins,&NeutronCandidates::NeutCand::GetTotalE),
	new NeutronVariable(("tracker"+nameTag+"_primary_parent").c_str(),"", myPDGBins,&NeutronCandidates::NeutCand::GetPDGBin),
	new NeutronVariable(("tracker"+nameTag+"_length").c_str(),"len. [mm]", myLenBins,&NeutronCandidates::NeutCand::GetLength),
	new NeutronVariable(("tracker"+nameTag+"_avg_dEdx").c_str(),"dE/dx [MeV/mm]", myDEDXBins,&NeutronCandidates::NeutCand::GetDEDX),
	new NeutronVariable(("tracker"+nameTag+"_dist").c_str(),"dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxDist),
	new NeutronVariable(("tracker"+nameTag+"_Zdist").c_str(),"Z dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxZDist),
      };

      for(auto& var: fLeadVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fLeadVars) var->InitializeDATAHists(data_error_bands);
      for(auto& var: fAllVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fAllVars) var->InitializeDATAHists(data_error_bands);
      for(auto& var: fTgtVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fTgtVars) var->InitializeDATAHists(data_error_bands);
      for(auto& var: fTrackVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fTrackVars) var->InitializeDATAHists(data_error_bands);
    }

    void SaveOrDraw(TDirectory& outDir){
      return;
    }

    void SaveOrDrawMC(TFile& outFile)
    {
      for (auto& var : fLeadVars) var->WriteMC(outFile);
      for (auto& var : fAllVars) var->WriteMC(outFile);
      for (auto& var : fTgtVars) var->WriteMC(outFile);
      for (auto& var : fTrackVars) var->WriteMC(outFile);
    }

    void SaveOrDrawData(TFile& outFile)
    {
      for (auto& var : fLeadVars) var->WriteData(outFile);
      for (auto& var : fAllVars) var->WriteData(outFile);
      for (auto& var : fTgtVars) var->WriteData(outFile);
      for (auto& var : fTrackVars) var->WriteData(outFile);
    }

  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      NeutronCandidates::NeutCand leadCand = evt.GetLeadingNeutCand();
      NeutronCandidates::NeutCands neutCands = evt.GetNeutCands();
      int intType = evt.GetIntType();
      int tgtType = evt.GetTgtZ();
      int leadBlobType = leadCand.GetPDGBin();

      std::bitset<64> SBStat = evt.GetSideBandStat();

      if (fAllCuts && !SBStat.all()) return;

      if (evt.IsMC()){
	if (evt.IsSignal()){
	  for (auto& var : fLeadVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	    (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	  }
	  for (auto& cand: neutCands.GetCandidates()){
	    for (auto& var : fAllVars){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	    }
	    if (cand.second.GetBegPos().Z() <= fBound){
	      for (auto& var : fTgtVars){
		var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      }
	    }
	    else {
	      for (auto& var : fTrackVars){
		var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      }
	    }
	  }
	}
	
	else {
	  int bkgd_ID = -1;
	  bkgd_ID = util::GetBackgroundID(univ);

	  for (auto& var : fLeadVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	    (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	  }
	  for (auto& cand: neutCands.GetCandidates()){
	    for (auto& var : fAllVars){
	      var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	    }
	    if (cand.second.GetBegPos().Z() <= fBound){
	      for (auto& var : fTgtVars){
		var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      }
	    }
	    else {
	      for (auto& var : fTrackVars){
		var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
		(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      }
	    }
	  }
	}
      }
      else{
	for (auto& var : fLeadVars){
	  var->dataHist->FillUniverse(&univ, var->GetRecoValue(leadCand), 1);
	}
	for (auto& cand: neutCands.GetCandidates()){
	  for (auto& var : fAllVars){
	    var->dataHist->FillUniverse(&univ, var->GetRecoValue(cand.second), 1);
	  }
	  if (cand.second.GetBegPos().Z() <= fBound){
	    for (auto& var : fTgtVars){
	      var->dataHist->FillUniverse(&univ, var->GetRecoValue(cand.second), 1);	     
	    }
	  }
	  else {
	    for (auto& var : fTrackVars){
	      var->dataHist->FillUniverse(&univ, var->GetRecoValue(cand.second), 1);
	    }
	  }
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
