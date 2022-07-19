//studies includes
#include "studies/Study.h"

//David's includes
#include "event/NeutCands.h"
#include "util/NeutronVariable.h"
#include "util/EventVariable.h"
#include "util/Variable.h"
#include "util/Variable2D.h"
#include "util/GetBackgroundID.h"
#include "util/Categorized.h"
#include "event/CVUniverse.h"

//Renamed to PreRecoil, includes all target breakdowns needed to fit plastic backgrounds etc.
class PreRecoil: public Study
{
  private:
    std::vector<Variable*> fVars;
    std::vector<util::Categorized<Variable, int>*> fVars_US_ByTgt;//US Plastic inside target
    std::vector<util::Categorized<Variable, int>*> fVars_DS_ByTgt;//DS Plastic inside target
    std::vector<util::Categorized<Variable, int>*> fVars_US_Post_ByTgt;//See above US defined object, but before recoil cut
    std::vector<util::Categorized<Variable, int>*> fVars_DS_Post_ByTgt;//See immediately above, but DS Plastic
    std::vector<util::Categorized<Variable, int>*> fVars_ByTgt;
    std::map<int, util::Categorized<Variable,int>*> fVtx_US_ByTgt;
    std::map<int, util::Categorized<Variable,int>*> fVtx_DS_ByTgt;
    std::vector<Variable2D*> fVars2D;
    std::map<int,Variable*> fRecoilBinned;
    TString fFVregionName;
    bool fSplitRecoil;

  public:
    PreRecoil(std::vector<Variable*> vars,
		     std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		     std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
	      std::vector<CVUniverse*>& data_error_bands, bool splitRecoil, TString FVregionName): Study(), fSplitRecoil(splitRecoil), fFVregionName(FVregionName)
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
	  fRecoilBinned[iBin] = new Variable(false,("recoilE_bin_"+binName).c_str(), "Recoil E [GeV]", myRecoilBins, &CVUniverse::GetDANRecoilEnergyGeV);
	  fRecoilBinned[iBin]->InitializeMCHists(mc_error_bands, truth_error_bands);
	  fRecoilBinned[iBin]->InitializeDATAHists(data_error_bands);
	}
      }

      std::vector<double> coarseVtxBins;//Modified to be a number of planes instead of vertex Z values.
      for (int iBin=0; iBin<22;++iBin){
	coarseVtxBins.push_back((double)iBin-10.5);
      }

      fVars_ByTgt = {};
      fVtx_US_ByTgt = {};
      fVtx_DS_ByTgt = {};
      for (auto& var : vars){
	fVars.push_back(new Variable(false, (var->GetName()+"_PreRecoilCut").c_str(), var->GetAxisLabel(), var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	if (fFVregionName == "Targets"){
	  TString nameCheck = var->GetName();
	  if (nameCheck=="vtxZ"){
	    for (auto& tgt: util::TgtCodeList){
	      fVtx_US_ByTgt[tgt.first]=new util::Categorized<Variable, int>(("ByTgt_"+tgt.second).c_str(), "US_ByType", false, (var->GetName()+"_ByTgt_"+tgt.second+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtTypeList,coarseVtxBins,var->GetRecoFunc(),var->GetTrueFunc());
	      fVtx_DS_ByTgt[tgt.first]=new util::Categorized<Variable, int>(("ByTgt_"+tgt.second).c_str(), "DS_ByType", false, (var->GetName()+"_ByTgt_"+tgt.second+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtTypeList,coarseVtxBins,var->GetRecoFunc(),var->GetTrueFunc());
	    }
	  }
	  else if (nameCheck == "pTmu" || nameCheck == "recoilE"){
	    fVars_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList,var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	    fVars_US_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_InnerUSPlastic_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList,var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	    fVars_DS_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_InnerDSPlastic_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList,var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	    fVars_US_Post_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_InnerUSPlastic").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList,var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	    fVars_DS_Post_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_InnerDSPlastic").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList,var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	  }
	  else{
	    fVars_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", false, (var->GetName()+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList,var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	  }
	}
      }
      
      //Commented out in the main code so commented here as well.
      fVars2D = {
	//new Variable2D(*fVars[4],*fVars[3]),//recoil v. Q2     
	//new Variable2D(*fVars[fVars.size()-1],*fVars[0]),//pT v. recoilQ2Bin     
      };

      for(auto& var: fVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars) var->InitializeDATAHists(data_error_bands);

      for(auto& tgt: fVars_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars_ByTgt){
	tgt->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVars_US_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars_US_ByTgt){
	tgt->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVars_DS_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars_DS_ByTgt){
	tgt->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVars_US_Post_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars_US_Post_ByTgt){
	tgt->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVars_DS_Post_ByTgt){
	tgt->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVars_DS_Post_ByTgt){
	tgt->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVtx_US_ByTgt){
	tgt.second->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVtx_US_ByTgt){
	tgt.second->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fVtx_DS_ByTgt){
	tgt.second->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fVtx_DS_ByTgt){
	tgt.second->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

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
      for(auto& tgt: fVars_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars_US_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars_DS_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars_US_Post_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVars_DS_Post_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVtx_US_ByTgt){
	tgt.second->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fVtx_DS_ByTgt){
	tgt.second->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for (auto& var : fVars2D) var->Write(outFile);
    }

    void SaveOrDrawData(TFile& outFile)
    {
      for (auto& recoil : fRecoilBinned) recoil.second->WriteData(outFile);
      for (auto& var : fVars) var->WriteData(outFile);
      for(auto& tgt: fVars_ByTgt){
	tgt->visit([&outFile](Variable& var)
		   {
		     var.WriteData(outFile);
		   });
      }
      for(auto& tgt: fVtx_US_ByTgt){
	tgt.second->visit([&outFile](Variable& var)
		   {
		     var.WriteData(outFile);
		   });
      }
      for(auto& tgt: fVtx_DS_ByTgt){
	tgt.second->visit([&outFile](Variable& var)
		   {
		     var.WriteData(outFile);
		   });
      }
    }
    
  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      std::bitset<64> SBStat = evt.GetSideBandStat();
      int intType = evt.GetIntType();
      int intCode = evt.GetIntCode();
      int tgtType = evt.GetTgtZ();
      int leadBlobType = evt.GetLeadingNeutCand().GetPDGBin();
      int iBin = evt.GetBinPTPZ();
      const bool isTgts = (fVars_ByTgt.size() > 0) ? true: false;
      std::vector<double> muonMom = {univ.GetMuon4V().X(),univ.GetMuon4V().Y(),univ.GetMuon4V().Z()};

      // At some point just put this into the event structure itself.
      std::vector<double> vtx = univ.GetVtx();
      double vtx_x = vtx.at(0);
      double vtx_y = vtx.at(1);
      double vtx_z = vtx.at(2);

      int tgtID = util::GetRecoTargetZ(vtx_x,vtx_y,vtx_z);
      int tgtCode = util::GetRecoTargetCode(vtx_x,vtx_y,vtx_z,muonMom);
      int tgtZ = univ.GetTargetZ();
      int USTgt = -1;
      int DSTgt = -1;
      if (tgtCode == -1){
	USTgt = util::GetUSTgtCode(tgtID,vtx_x,vtx_y,vtx_z,muonMom);
	DSTgt = util::GetDSTgtCode(tgtID,vtx_x,vtx_y,vtx_z,muonMom);
	/*Debugging removed
	if (tgtID >= 10 && USTgt == -999 && tgtType == 1){
	  std::cout << "TGT ID: " << tgtID << std::endl;
	  std::cout << "US Tgt Code: " << USTgt << std::endl;
	  std::cout << "DS Tgt Code: " << DSTgt << std::endl;
	  std::cout << "Event Weight: " << weight << std::endl;
	  if(evt.IsSignal()) std::cout << "Signal" << std::endl;
	  else std::cout << "BKG" << std::endl;
	  if (tgtID == 10) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(1111,vtx_z) << std::endl;
	  else if (tgtID == 21) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(3333,vtx_z) << std::endl;
	  else if (tgtID == 32) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(3333,vtx_z) << std::endl;
	  else if (tgtID == 63) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(6666,vtx_z) << std::endl;
	  else if (tgtID == 46) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(4444,vtx_z) << std::endl;
	  else if (tgtID == 54) std::cout << "N Planes US: " << util::GetNPlanesUSOfTarget(5555,vtx_z) << std::endl;
	  std::cout << "" << std::endl;
	}
	else if ((tgtID == 0 || tgtID > 11) && DSTgt == -999 && tgtType == 1){
	  std::cout << "TGT ID: " << tgtID << std::endl;
	  std::cout << "US Tgt Code: " << USTgt << std::endl;
	  std::cout << "DS Tgt Code: " << DSTgt << std::endl;
	  std::cout << "Event Weight: " << weight << std::endl;
	  if(evt.IsSignal()) std::cout << "Signal" << std::endl;
	  else std::cout << "BKG" << std::endl;

	  if(tgtID==21) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(1111,vtx_z) << std::endl;
	  else if(tgtID==32) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(2222,vtx_z) << std::endl;
	  else if(tgtID==63) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(3333,vtx_z) << std::endl;
	  else if(tgtID==46) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(6666,vtx_z) << std::endl;
	  else if(tgtID==54) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(4444,vtx_z) << std::endl;
	  else if(tgtID==0) std::cout << "N Planes DS: " << util::GetNPlanesDSOfTarget(5555,vtx_z) << std::endl;
	  std::cout << "" << std::endl;
	}
	*/
      }
      
      if (evt.IsMC()){
	if (evt.IsSignal()){
	  if (fSplitRecoil){
	    fRecoilBinned[iBin]->selectedMCReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight); //"Fake data" for closure
	    fRecoilBinned[iBin]->selectedSignalReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_SigIntTypeHists)[intType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    //(*fRecoilBinned[iBin]->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	  }

	  for (auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }

	  //Only fill for non-interstitial plastics
	  if(tgtCode != -1){
	    for (auto& var: fVars_ByTgt){
	      (*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
	      (*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
	      (*(*var)[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
	      (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
	      //(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
	    }
	  }
	  else{
	    if(USTgt > 0){
	      /* Debugging removed
	      std::cout << "IsSignal" << std::endl;
	      std::cout << "Material: " << tgtType << std::endl;
	      std::cout << "USTgt: " << USTgt << std::endl;
	      std::cout << "Vertex Z: " << vtx_z << std::endl;
	      std::cout << "Variable Vtx Z: " << (*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ) << std::endl;
	      std::cout << "NPlanes US for Target: " << util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
	      std::cout << "" << std::endl;
	      */
	      (*(fVtx_US_ByTgt[USTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
	      (*(fVtx_US_ByTgt[USTgt]))[tgtType].selectedSignalReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
	      (*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_SigIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
	      (*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
	      //(*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);	    
	    }
	    if(DSTgt > 0){
	      /*Debugging removed
	      std::cout << "IsSignal" << std::endl;
	      std::cout << "Material: " << tgtType << std::endl;
	      std::cout << "DSTgt: " << DSTgt << std::endl;
	      std::cout << "Vertex Z: " << vtx_z << std::endl;
	      std::cout << "Variable Vtx Z: " << (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ) << std::endl;
	      std::cout << "NPlanes DS for Target: " << util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
	      std::cout << "" << std::endl;
	      */
	      (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
	      (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].selectedSignalReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
	      (*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_SigIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
	      (*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
	      //(*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);	    
	    }
	  }

	  for (auto& var: fVars2D){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	  }
	}
	
	else{
	  int bkgd_ID = -1;
	  if (isTgts) bkgd_ID = 44;
	  if (tgtType == 200 || tgtType == 300){
	    if (evt.IsFSSignal()){
	      if(tgtCode != -1){
		if (tgtType == 200){
		  for (auto& var: fVars_US_ByTgt){
		    (*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
		    (*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    (*(*var)[tgtCode].m_SigIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    //(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		  }
		  if (SBStat.all()){
		    for (auto& var: fVars_US_Post_ByTgt){
		      (*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
		      (*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      (*(*var)[tgtCode].m_SigIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      //(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    }
		  }
		}
		else{
		  for (auto& var: fVars_DS_ByTgt){
		    (*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
		    (*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    (*(*var)[tgtCode].m_SigIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    //(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		  }
		  if (SBStat.all()){
		    for (auto& var: fVars_DS_Post_ByTgt){
		      (*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
		      (*var)[tgtCode].selectedSignalReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      (*(*var)[tgtCode].m_SigIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      //(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    }
		  }
		}
	      }
	    }
	    else{
	      bkgd_ID = util::GetBackgroundID(univ);
	      if(tgtCode != -1){
		if (tgtType == 200){
		  for(auto& var: fVars_US_ByTgt){
		    (*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
		    (*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    (*(*var)[tgtCode].m_BkgIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    //(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		  }
		  if (SBStat.all()){
		    for(auto& var: fVars_US_Post_ByTgt){
		      (*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
		      (*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      (*(*var)[tgtCode].m_BkgIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      //(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    }
		  }
		}
		else{
		  for(auto& var: fVars_DS_ByTgt){
		    (*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
		    (*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    (*(*var)[tgtCode].m_BkgIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    //(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		  }
		  if (SBStat.all()){
		    for(auto& var: fVars_DS_Post_ByTgt){
		      (*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
		      (*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      (*(*var)[tgtCode].m_BkgIntTypeHists)[intCode].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		      //(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
		    }
		  }
		}
	      }
	    }
	    bkgd_ID = intType;
	  }
	  else if (util::CorrectTargetMaterial(tgtCode,tgtZ)) bkgd_ID = util::GetBackgroundID(univ);
	  
	  if (fSplitRecoil){
	    fRecoilBinned[iBin]->selectedMCReco->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight); //"Fake data" for closure
	    (*fRecoilBinned[iBin]->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_BkgIntTypeHists)[intType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    (*fRecoilBinned[iBin]->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	    //(*fRecoilBinned[iBin]->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	  }

	  for(auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }

	  //Only fill for non-interstitial plastics
	  if(tgtCode != -1){
	    for(auto& var: fVars_ByTgt){
	      (*var)[tgtCode].selectedMCReco->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight); //"Fake data" for closure
	      (*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
	      (*(*var)[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
	      (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
	      //(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), weight);
	    }
	  }
	  else{
	    if(USTgt > 0){
	      /* Debugging removed
	      std::cout << "IsBKG" << std::endl;
	      std::cout << "Material: " << tgtType << std::endl;
	      std::cout << "USTgt: " << USTgt << std::endl;
	      std::cout << "Vertex Z: " << vtx_z << std::endl;
	      std::cout << "Variable Vtx Z: " << (*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ) << std::endl;
	      std::cout << "NPlanes US for Target: " << util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
	      std::cout << "" << std::endl;
	      */
	      (*(fVtx_US_ByTgt[USTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
	      (*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
	      (*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_BkgIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
	      (*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
	      //(*(*(fVtx_US_ByTgt[USTgt]))[tgtType].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), weight);
	    }
	    if(DSTgt > 0){
	      /*Debugging removed
	      std::cout << "IsBKG" << std::endl;
	      std::cout << "Material: " << tgtType << std::endl;
	      std::cout << "DSTgt: " << DSTgt << std::endl;
	      std::cout << "Vertex Z: " << vtx_z << std::endl;
	      std::cout << "Variable Vtx Z: " << (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ) << std::endl;
	      std::cout << "NPlanes DS for Target: " << util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
	      std::cout << "" << std::endl;
	      */
	      (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight); //"Fake data" for closure
	      (*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
	      (*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_BkgIntTypeHists)[intType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
	      (*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
	      //(*(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), weight);
	    }
	  }

	  for(auto& var: fVars2D){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	  }
	}
      }
      
      else{
	if (fSplitRecoil) fRecoilBinned[iBin]->dataHist->FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), 1);

	for (auto& var : fVars){
	  var->dataHist->FillUniverse(&univ, var->GetRecoValue(univ), 1);
	}

	//Only fill for non-interstitial plastics
	if(tgtCode != -1){
	  for (auto& var : fVars_ByTgt){
	    (*var)[tgtCode].dataHist->FillUniverse(&univ, (*var)[tgtCode].GetRecoValue(univ), 1);
	  }
	}

	/* REMOVE FILLING THE PLASTIC DATA. It's useless.
	else{
	  if(USTgt > 0){
	    std::cout << "IsData" << std::endl;
	    std::cout << "Material: " << tgtType << std::endl;
	    std::cout << "USTgt: " << USTgt << std::endl;
	    std::cout << "Vertex Z: " << vtx_z << std::endl;
	    std::cout << "Variable Vtx Z: " << (*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ) << std::endl;
	    std::cout << "NPlanes US for Target: " << util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
	    std::cout << "" << std::endl;
	    (*(fVtx_US_ByTgt[USTgt]))[tgtType].dataHist->FillUniverse(&univ, util::GetNPlanesUSOfTarget(USTgt,(*(fVtx_US_ByTgt[USTgt]))[tgtType].GetRecoValue(univ)), 1);
	  }
	  if(DSTgt > 0){
	    std::cout << "IsData" << std::endl;
	    std::cout << "Material: " << tgtType << std::endl;
	    std::cout << "DSTgt: " << DSTgt << std::endl;
	    std::cout << "Vertex Z: " << vtx_z << std::endl;
	    std::cout << "Variable Vtx Z: " << (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ) << std::endl;
	    std::cout << "NPlanes DS for Target: " << util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)) << std::endl;
	    std::cout << "" << std::endl;
	    (*(fVtx_DS_ByTgt[DSTgt]))[tgtType].dataHist->FillUniverse(&univ, util::GetNPlanesDSOfTarget(DSTgt,(*(fVtx_DS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ)), 1);
	  }
	}
	*/

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
	  fRecoilBinned[iBin] = new Variable(false,("recoilE_RecoilSB_bin_"+binName).c_str(), "Recoil E [GeV]", myRecoilBins, &CVUniverse::GetDANRecoilEnergyGeV);
	  fRecoilBinned[iBin]->InitializeMCHists(mc_error_bands, truth_error_bands);
	  fRecoilBinned[iBin]->InitializeDATAHists(data_error_bands);
	}
      }

      for (auto& var : vars){
	fVars.push_back(new Variable(false,(var->GetName()+"_RecoilSB").c_str(), var->GetAxisLabel(), var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
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
	    //(*fRecoilBinned[iBin]->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	  }

	  for (auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }

	  for (auto& var: fVars2D){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    //(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
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
	    //(*fRecoilBinned[iBin]->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, fRecoilBinned[iBin]->GetRecoValue(univ), weight);
	  }

	  for(auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }

	  for(auto& var: fVars2D){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
	    //(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValueX(univ), var->GetRecoValueY(univ), weight);
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
