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

//Renamed to PreRecoil
class PreRecoil: public Study
{
  private:
    std::vector<Variable*> fVars;
    std::vector<util::Categorized<Variable, int>*> fVars_ByTgt;
    std::map<int, util::Categorized<Variable,int>*> fUS_ByTgt;
    std::map<int, util::Categorized<Variable,int>*> fDS_ByTgt;
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

      std::vector<double> coarseVtxBins;
      for (int iBin=0; iBin<201;++iBin){
	coarseVtxBins.push_back((double)iBin*2000.0/200.0+4000.0);
      }

      fVars_ByTgt = {};
      fUS_ByTgt = {};
      fDS_ByTgt = {};
      for (auto& var : vars){
	fVars.push_back(new Variable(false, (var->GetName()+"_PreRecoilCut").c_str(), var->GetAxisLabel(), var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
	if (fFVregionName == "Targets"){
	  TString nameCheck = var->GetName();
	  if (nameCheck=="vtxZ"){
	    for (auto& tgt: util::TgtCodeList){
	      fUS_ByTgt[tgt.first]=new util::Categorized<Variable, int>(("ByTgt_"+tgt.second).c_str(), "US_ByType", var->IsAnaVar(), (var->GetName()+"_ByTgt_"+tgt.second+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtTypeList,coarseVtxBins,var->GetRecoFunc(),var->GetTrueFunc());
	      fDS_ByTgt[tgt.first]=new util::Categorized<Variable, int>(("ByTgt_"+tgt.second).c_str(), "DS_ByType", var->IsAnaVar(), (var->GetName()+"_ByTgt_"+tgt.second+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtTypeList,coarseVtxBins,var->GetRecoFunc(),var->GetTrueFunc());
	    }
	  }
	  else{
	    fVars_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", var->IsAnaVar(), (var->GetName()+"_PreRecoilCut").c_str(),var->GetAxisLabel().c_str(), util::TgtCodeList,var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
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

      for(auto& tgt: fUS_ByTgt){
	tgt.second->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fUS_ByTgt){
	tgt.second->visit([&data_error_bands](Variable& var)
		   {
		     var.InitializeDATAHists(data_error_bands);
		   });
      }

      for(auto& tgt: fDS_ByTgt){
	tgt.second->visit([&mc_error_bands, &truth_error_bands](Variable& var)
		   {
		     var.InitializeMCHists(mc_error_bands, truth_error_bands);
		   });
      }
      for(auto& tgt: fDS_ByTgt){
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
      for(auto& tgt: fUS_ByTgt){
	tgt.second->visit([&outFile](Variable& var)
		   {
		     var.WriteMC(outFile);
		   });
      }
      for(auto& tgt: fDS_ByTgt){
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
      for(auto& tgt: fUS_ByTgt){
	tgt.second->visit([&outFile](Variable& var)
		   {
		     var.WriteData(outFile);
		   });
      }
      for(auto& tgt: fDS_ByTgt){
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
      int tgtType = evt.GetTgtZ();
      int leadBlobType = evt.GetLeadingNeutCand().GetPDGBin();
      int iBin = evt.GetBinPTPZ();

      // At some point just put this into the event structure itself.
      std::vector<double> vtx = univ.GetVtx();
      double vtx_x = vtx.at(0);
      double vtx_y = vtx.at(1);
      double vtx_z = vtx.at(2);

      int tgtID = util::GetRecoTargetZ(vtx_x,vtx_y,vtx_z);
      int tgtCode = util::GetRecoTargetCode(vtx_x,vtx_y,vtx_z);
      int USTgt = -1;
      int DSTgt = -1;
      if (tgtCode == -1){
	USTgt = util::USTgtMap[tgtID];
	DSTgt = util::DSTgtMap[tgtID];
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
	    if(USTgt != -1){
	      (*(fUS_ByTgt[USTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), weight); //"Fake data" for closure
	      (*(fUS_ByTgt[USTgt]))[tgtType].selectedSignalReco->FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), weight);
	      (*(*(fUS_ByTgt[USTgt]))[tgtType].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), weight);
	      (*(*(fUS_ByTgt[USTgt]))[tgtType].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), weight);
	      //(*(*(fUS_ByTgt[USTgt]))[tgtType].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), weight);	    
	    }
	    if(DSTgt != -1){
	      (*(fDS_ByTgt[DSTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), weight); //"Fake data" for closure
	      (*(fDS_ByTgt[DSTgt]))[tgtType].selectedSignalReco->FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), weight);
	      (*(*(fDS_ByTgt[DSTgt]))[tgtType].m_SigIntTypeHists)[intType].FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), weight);
	      (*(*(fDS_ByTgt[DSTgt]))[tgtType].m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), weight);
	      //(*(*(fDS_ByTgt[DSTgt]))[tgtType].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), weight);	    
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
	    if(USTgt != -1){
	      (*(fUS_ByTgt[USTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), weight); //"Fake data" for closure
	      (*(*(fUS_ByTgt[USTgt]))[tgtType].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), weight);
	      (*(*(fUS_ByTgt[USTgt]))[tgtType].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), weight);
	      (*(*(fUS_ByTgt[USTgt]))[tgtType].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), weight);
	      //(*(*(fUS_ByTgt[USTgt]))[tgtType].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), weight);
	    }
	    if(DSTgt != -1){
	      (*(fDS_ByTgt[DSTgt]))[tgtType].selectedMCReco->FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), weight); //"Fake data" for closure
	      (*(*(fDS_ByTgt[DSTgt]))[tgtType].m_backgroundHists)[bkgd_ID].FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), weight);
	      (*(*(fDS_ByTgt[DSTgt]))[tgtType].m_BkgIntTypeHists)[intType].FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), weight);
	      (*(*(fDS_ByTgt[DSTgt]))[tgtType].m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), weight);
	      //(*(*(fDS_ByTgt[DSTgt]))[tgtType].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), weight);
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
	else{
	  if(USTgt != -1){
	    (*(fUS_ByTgt[USTgt]))[tgtType].dataHist->FillUniverse(&univ, (*(fUS_ByTgt[USTgt]))[tgtType].GetRecoValue(univ), 1);
	  }
	  if(DSTgt != -1){
	    (*(fDS_ByTgt[DSTgt]))[tgtType].dataHist->FillUniverse(&univ, (*(fDS_ByTgt[DSTgt]))[tgtType].GetRecoValue(univ), 1);
	  }
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
