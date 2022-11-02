#define MC_OUT_FILE_NAME "runEventLoopMC"
#define DATA_OUT_FILE_NAME "runEventLoopData"

#define USAGE \
"\n*** USAGE ***\n"\
"runEventLoop <dataPlaylist.txt> <mcPlaylist.txt> <skip_systematics> <MnvTune_v> <FV> <TgtNum> optional: <fileName_tag> <NeutKE>\n\n"\
"*** Explanation ***\n"\
"Reduce MasterAnaDev AnaTuples to event selection histograms to extract a\n"\
"single-differential inclusive cross section for the 2021 MINERvA 101 tutorial.\n\n"\
"*** The Input Files ***\n"\
"Playlist files are plaintext files with 1 file name per line.  Filenames may be\n"\
"xrootd URLs or refer to the local filesystem.  The first playlist file's\n"\
"entries will be treated like data, and the second playlist's entries must\n"\
"have the \"Truth\" tree to use for calculating the efficiency denominator.\n\n"\
"*** Output ***\n"\
"Produces a two files starting with names, " MC_OUT_FILE_NAME " and " DATA_OUT_FILE_NAME ", with\n"\
"all histograms needed for the ExtractCrossSection program also built by this\n"\
"package.  You'll need a .rootlogon.C that loads ROOT object definitions from\n"\
"PlotUtils to access systematics information from these files.\n\n"\
"*** Environment Variables and Inputs ***\n"\
"Setting up this package appends to PATH and LD_LIBRARY_PATH.  PLOTUTILSROOT,\n"\
"MPARAMFILESROOT, and MPARAMFILES must be set according to the setup scripts in\n"\
"those packages for systematics and flux reweighters to function.\n"\
"If the <skip_systematics> argument is nonzero at all, output histograms will have no error bands.\n"\
"This is useful for debugging the CV and running warping studies.\n"\
"<MnvTune_v> has to be 1 or 2 and is the only current input control for the arguments to the function.\n"\
"<FV> parameter needs to be set to Tracker or Targets to select vtx location.\n"\
"<TgtNum> parameter needs to be set -1 for all targets or tracker, or 1-6 for those targets (6 is water tgt).\n"\
"This parameter also controls the plotting range of the vtxZ variable.\n\n"\
"*** Optional Inputs ***\n"\
"<fileName_tag> lives as a naming tag used for extra identification of file contents if desired.\n"\
"<NeutKE> If <=0, skips neutron cuts for deliverables that may not use them.\n"\
"If >0, then this is the requirement in MeV for a signal neutron.\n\n"\
"*** Return Codes ***\n"\
"0 indicates success.  All histograms are valid only in this case.  Any other\n"\
"return code indicates that histograms should not be used.  Error messages\n"\
"about what went wrong will be printed to stderr.  So, they'll end up in your\n"\
"terminal, but you can separate them from everything else with something like:\n"\
"\"runEventLoop data.txt mc.txt 2> errors.txt\"\n"

enum ErrorCodes
{
  success = 0,
  badCmdLine = 1,
  badInputFile = 2,
  badFileRead = 3,
  badOutputFile = 4
};

//PlotUtils includes
//No junk from PlotUtils please!  I already
//know that MnvH1D does horrible horrible things.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"

//Includes from this package
#include "event/CVUniverse.h"
#include "event/MichelEvent.h"
#include "event/NeutCands.h"
#include "systematics/Systematics.h"
#include "cuts/MaxPzMu.h"
#include "cuts/CCQECuts.h"
#include "cuts/NeutCuts.h"
#include "util/Variable.h"
#include "util/NeutronVariable.h"
#include "util/Variable2D.h"
#include "util/GetFluxIntegral.h"
#include "util/GetPlaylist.h"
#include "util/GetBackgroundID.h"
#include "util/GetRecoTargetZ.h"
#include "util/Categorized.h"
#include "cuts/SignalDefinition.h"
#include "cuts/q3RecoCut.h"
#include "studies/Study.h"
#include "studies/NeutronVariables.h"
#include "studies/EMSideBands.h"
#include "studies/MichelAndNBlobSB.h"
#include "studies/RecoilSB.h"
#include "studies/LeadNeutStudy.h"
//#include "Binning.h" //TODO: Fix me

//PlotUtils includes
#include "PlotUtils/makeChainWrapper.h"
#include "PlotUtils/HistWrapper.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MacroUtil.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/CCInclusiveCuts.h"
#include "PlotUtils/CCInclusiveSignal.h"
#include "PlotUtils/CrashOnROOTMessage.h" //Sets up ROOT's debug callbacks by itself
#include "PlotUtils/Cutter.h"
#include "PlotUtils/Model.h"
#include "PlotUtils/FluxAndCVReweighter.h"
#include "PlotUtils/GENIEReweighter.h"
#include "PlotUtils/GeantNeutronCVReweighter.h"
#include "PlotUtils/LowRecoil2p2hReweighter.h"
#include "PlotUtils/LowQ2PiReweighter.h"
#include "PlotUtils/RPAReweighter.h"
#include "PlotUtils/MINOSEfficiencyReweighter.h"
#include "PlotUtils/TargetUtils.h"
#pragma GCC diagnostic pop

//ROOT includes
#include "TParameter.h"
#include "TMath.h"

//c++ includes
#include <iostream>
#include <cstdlib> //getenv()

//==============================================================================
// Loop and Fill
//==============================================================================
void LoopAndFillEventSelection(
    PlotUtils::ChainWrapper* chain,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    std::vector<Variable*> vars,
    std::vector<util::Categorized<Variable, int>*> vars_ByTgt,
    std::vector<Variable2D*> vars2D,
    std::vector<util::Categorized<Variable2D, int>*> vars2D_ByTgt,
    std::vector<Study*> studies,
    PlotUtils::Cutter<CVUniverse, NeutronEvent>& michelcuts,
    PlotUtils::Model<CVUniverse, NeutronEvent>& model,
    bool doNeutron = true)
{
  assert(!error_bands["cv"].empty() && "\"cv\" error band is empty!  Can't set Model weight.");
  auto& cvUniv = error_bands["cv"].front();

  std::cout << "Starting MC reco loop...\n";
  const int nEntries = chain->GetEntries();
  for (int i=0; i<nEntries; ++i)
  {
    if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::endl;

    cvUniv->SetEntry(i);
    //NeutronEvent cvEvent(cvUniv->GetLeadNeutCandOnly());
    NeutronEvent cvEvent = doNeutron ? NeutronEvent(cvUniv->GetLeadNeutCandOnly()) : NeutronEvent();
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);
    //For testing.
    //const double cvWeight = 1.0;

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : error_bands)
    {
      std::vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes)
      {    
        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);
        
        // This is where you would Access/create a Michel
        //NeutronEvent myevent(universe->GetLeadNeutCandOnly()); // make sure your event is inside the error band loop. 
	NeutronEvent myevent = doNeutron ? NeutronEvent(universe->GetLeadNeutCandOnly()) : NeutronEvent();
	myevent.SetIsMC();

	myevent.SetEMBlobInfo(universe->GetEMNBlobsTotalEnergyTotalNHits());
	std::bitset<64> SBStat = michelcuts.isMCSelected(*universe, myevent, cvWeight);
	myevent.SetSideBandStat(SBStat);

	if (SBStat.none()) continue;

        //weight is ignored in isMCSelected() for all but the CV Universe.
        //if (!michelcuts.isMCSelected(*universe, myevent, cvWeight).all()) continue; //all is another function that will later help me with sidebands
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the per-universe weight for events that will actually use it.
        //const double weight = 1.0; //Dummy weight for testing/validation pre-weight
        const bool isFSSignal = michelcuts.isSignal(*universe, weight);
	const bool isTgts = (vars_ByTgt.size() > 0 || vars2D_ByTgt.size()) ? true : false; 
	bool tmpIsSignal = isFSSignal;
	int intType = universe->GetInteractionType();
	int tgtZ = universe->GetTargetZ();
	//util::GetTargetType(tgtZ,);
	myevent.SetIntCode(intType);

	//get muon momentum to project into targets
	std::vector<double> muonMom = {universe->GetMuon4V().X(),universe->GetMuon4V().Y(),universe->GetMuon4V().Z()};

	std::vector<double> vtx = universe->GetVtx();
	double vtx_x = vtx.at(0);
	double vtx_y = vtx.at(1);
	double vtx_z = vtx.at(2);

	/*
	if ((TString)(universe->ShortName()) == "cv" && vtx_z > 5650 && vtx_z < 5800){
	  std::cout << "Event Vertex Z: " << vtx_z << std::endl;
	  }*/

	std::vector<double> mc_vtx = universe->GetTrueVtx();
	double mc_vtx_x = mc_vtx.at(0);
	double mc_vtx_y = mc_vtx.at(1);
	double mc_vtx_z = mc_vtx.at(2);

	int tgtID = util::GetRecoTargetZ(vtx_x,vtx_y,vtx_z);
	int tgtCode = util::GetRecoTargetCode(vtx_x,vtx_y,vtx_z,muonMom);
	//int trueTgtID = util::GetTightTargetZ(mc_vtx_x,mc_vtx_y,mc_vtx_z);//Get Target By Z for the true vtx. Tight constraints to not include US or DS plastic.
	int trueTgtCode = util::GetTrueTgtCode(tgtZ, mc_vtx_x, mc_vtx_y, mc_vtx_z);

	int tgtType = tgtZ;
	if (tgtZ==6 && trueTgtCode != 3306) tgtType = 1; //Carbon not in target 3 is in plastic... Maybe need to require trueTgtCode to be carbon... but that can wait on the assumption there's no carbon outside of the carbon target.
	if (tgtZ==1 && trueTgtCode == 6666) tgtType = 8; //hydrogen in target 6 (water target) is water
	if (tgtZ==8 && trueTgtCode != 6666) tgtType = -1; //oxygen not in the water target is ill-defined

	if (tgtType==1 && tgtID > 0 && tgtID < 7){
	  if (mc_vtx_z < vtx_z){
	    tgtType = 200; //US Plastic has true vertex upstream
	    intType = 62;
	  }
	  if (mc_vtx_z > vtx_z){
	    tgtType = 300; //DS Plastic has true vertex downstream
	    intType = 63;
	  }
	}

	if (tgtID > 0 && tgtID < 7){
	  tmpIsSignal = tmpIsSignal ? util::CorrectTargetMaterial(tgtCode,trueTgtCode) : tmpIsSignal; //Only call it signal if the true material for that section is correct.
	}

	const bool isSignal = tmpIsSignal;

	/* Just for checking the output of GetRecoTargetZ
	std::cout << "Checking Target Breakdown" << std::endl;
	std::cout << "Target Z " << tgtType << std::endl;
	std::cout << "Reco Vtx.: " << vtx.at(0) << ", " << vtx.at(1) << ", " << vtx.at(2) << std::endl;
	std::cout << "Reco Vtx. X: " << vtx_x << std::endl;
	std::cout << "Reco Vtx. Y: " << vtx_y << std::endl;
	std::cout << "Reco Vtx. Z: " << vtx_z << std::endl;
	std::cout << "Returned Tgt. ID: " << tgtID << std::endl;
	std::cout << "" << std::endl;
	*/

	myevent.SetFSSignal(isFSSignal);
        myevent.SetSignal(isSignal);
	myevent.SetIntType(intType);
	myevent.SetTgtZ(tgtType);

	myevent.SetBinPTPZ(universe->GetPTPZBin());

	myevent.SetMaxFSNeutronKE(universe->GetMaxFSNeutronKE());
	//myevent.SetBKGID(bkgd_ID); Maybe add later...

	std::cout << std::setprecision(16);

	/*
	if ((TString)(universe->ShortName()) == "cv"){
	  std::cout << "Event: " << universe->GetDouble("eventID") << " has last bit: " << SBStat[0] << " and second-to-last bit: " << SBStat[1] << std::endl;
	  std::cout << "Has: " << myevent.GetEMNBlobs() << "EM Blobs" << std::endl;
	  std::cout << "With a ratio E/NHit of: " << myevent.GetEMBlobENHitRatio() << "EM Blobs" << std::endl;
	  std::cout << "And No. Michel of: " << universe->GetNImprovedMichel() << std::endl;
	  std::cout << "Int Type: " << intType << std::endl;
	  if (isSignal) std::cout << "Is Signal with lead true neutron KE: " << myevent.GetMaxFSNeutronKE() << std::endl;
	  else std::cout << "Is Not Signal." << std::endl;
	}
	*/

	int leadBlobType = myevent.GetLeadingNeutCand().GetPDGBin();
	/*
	if ((TString)(universe->ShortName()) == "cv"){
	  std::cout << "Event in MC CV: " << universe->GetDouble("eventID") << std::endl;
	  if (SBStat.all()) std::cout << "Selected" << std::endl;
	  else std::cout << "NOT Selected" << std::endl;
	  std::cout << "Int Type: " << intType << std::endl;
	  std::cout << "Lead Blob Type: " << leadBlobType << std::endl;
	  if (isSignal) std::cout << "Is Signal with lead true neutron KE: " << myevent.GetMaxFSNeutronKE() << std::endl;
	  else std::cout << "Is Not Signal w/ bkgd. ID: " << util::GetBackgroundID(*universe) << std::endl;
	  for (auto& var:vars) std::cout << "Var " << var->GetName() << ": " << var->GetRecoValue(*universe) << std::endl;
	  std::cout << "" << std::endl;
	}
	*/

	//if((TString)(universe->ShortName()) == "cv"){ ONLY FOR DEBUGGING
	for(auto& study: studies) study->Selected(*universe, myevent, weight);

	if (!SBStat.all()) continue;

        for(auto& var: vars) if(var->IsFill()) var->selectedMCReco->FillUniverse(universe, var->GetRecoValue(*universe), weight); //"Fake data" for closure

	for(auto& var: vars2D) if(var->IsFill()) var->selectedMCReco->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight); //"Fake data" for closure

	//Only fill by tgt for non-interstitial plastic events
	if (tgtCode != -1){
	  for(auto& var: vars_ByTgt) if((*var)[tgtCode].IsFill()) (*var)[tgtCode].selectedMCReco->FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), weight); //"Fake data" for closure
	  for(auto& var: vars2D_ByTgt) if((*var)[tgtCode].IsFill()) (*var)[tgtCode].selectedMCReco->FillUniverse(universe, (*var)[tgtCode].GetRecoValueX(*universe), (*var)[tgtCode].GetRecoValueY(*universe), weight); //"Fake data" for closure
	}

        if(isSignal)
        {
          //for(auto& study: studies) study->SelectedSignal(*universe, myevent, weight);

          for(auto& var: vars)
          {
            //Cross section components
	    if (var->IsFill()){
	      if (var->IsAnaVar()) var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValue(*universe), weight);
	      if (var->IsAnaVar()) var->migration->FillUniverse(universe, var->GetRecoValue(*universe), var->GetTrueValue(*universe), weight);
	      var->selectedSignalReco->FillUniverse(universe, var->GetRecoValue(*universe), weight); //Efficiency numerator in reco variables.  Useful for warping studies.

	      if (var->IsBroken()){
		//Various breakdowns of selected signal reco
		(*var->m_SigIntTypeHists)[intType].FillUniverse(universe, var->GetRecoValue(*universe), weight);
		(*var->m_SigTargetTypeHists)[tgtType].FillUniverse(universe, var->GetRecoValue(*universe), weight);
		(*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(universe, var->GetRecoValue(*universe), weight);
	      }
	    }
          }
	  
          for(auto& var: vars2D)
          {
	    if (var->IsFill()){
	      //Cross section components
	      var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
	      var->selectedSignalReco->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);; //Efficiency numerator in reco variables.  Useful for warping studies.

	      //Various breakdowns of selected signal reco
	      (*var->m_SigIntTypeHists)[intType].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);;
	      (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);;
	      (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
	      
	      //var->efficiencyNumerator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
	    }
	  }

	  //Only fill by tgt for non-interstitial plastic events
	  if (tgtCode != -1){
	    for(auto& var: vars_ByTgt){
	      if ((*var)[tgtCode].IsFill()){
		//Cross section components
		if ((*var)[tgtCode].IsAnaVar()) (*var)[tgtCode].efficiencyNumerator->FillUniverse(universe, (*var)[tgtCode].GetTrueValue(*universe), weight);
		if ((*var)[tgtCode].IsAnaVar()) (*var)[tgtCode].migration->FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), (*var)[tgtCode].GetTrueValue(*universe), weight);
		(*var)[tgtCode].selectedSignalReco->FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), weight); //Efficiency numerator in reco variables.  Useful for warping studies.

		if ((*var)[tgtCode].IsBroken()){
		  //Various breakdowns of selected signal reco
		  (*(*var)[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), weight);
		  (*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), weight);
		  (*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), weight);
		}
	      }
	    }
	    
	    for(auto& var: vars2D_ByTgt){
	      if ((*var)[tgtCode].IsFill()){
		//Cross section components
		(*var)[tgtCode].efficiencyNumerator->FillUniverse(universe, (*var)[tgtCode].GetTrueValueX(*universe), (*var)[tgtCode].GetTrueValueY(*universe), weight);
		(*var)[tgtCode].selectedSignalReco->FillUniverse(universe, (*var)[tgtCode].GetRecoValueX(*universe), (*var)[tgtCode].GetRecoValueY(*universe), weight);; //Efficiency numerator in reco variables.  Useful for warping studies.
		//Various breakdowns of selected signal reco
		(*(*var)[tgtCode].m_SigIntTypeHists)[intType].FillUniverse(universe, (*var)[tgtCode].GetRecoValueX(*universe), (*var)[tgtCode].GetRecoValueY(*universe), weight);;
		(*(*var)[tgtCode].m_SigTargetTypeHists)[tgtType].FillUniverse(universe, (*var)[tgtCode].GetRecoValueX(*universe), (*var)[tgtCode].GetRecoValueY(*universe), weight);;
		(*(*var)[tgtCode].m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(universe, (*var)[tgtCode].GetRecoValueX(*universe), (*var)[tgtCode].GetRecoValueY(*universe), weight);
	      
		//(*var)[tgtCode].efficiencyNumerator->FillUniverse(universe, (*var)[tgtCode].GetTrueValueX(*universe), (*var)[tgtCode].GetTrueValueY(*universe), weight);
	      }
	    }
	  }
        }

        else
        {
          int bkgd_ID = -1;
	  if (isTgts){
	    bkgd_ID = 44;
	    //bkgd_ID = util::GetBackgroundID(*universe);
	    if (tgtType == 200 || tgtType == 300){ //Separating out the plastic backgrounds
	      bkgd_ID = intType;
	    }
	    else if (util::CorrectTargetMaterial(tgtCode,tgtZ) || tgtCode == -1) bkgd_ID = util::GetBackgroundID(*universe);
	    if (bkgd_ID == 44) intType = bkgd_ID;
	  }
	  else bkgd_ID = util::GetBackgroundID(*universe);
          //for(auto& study: studies) study->SelectedSignal(*universe, myevent, weight);

          for(auto& var: vars){
	    if (var->IsFill()){
	      (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValue(*universe), weight);
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(universe, var->GetRecoValue(*universe), weight);
	      if (var->IsBroken()){
		//Various breakdowns of selected backgrounds
		(*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(universe, var->GetRecoValue(*universe), weight);
		(*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(universe, var->GetRecoValue(*universe), weight);
	      }
	    }
	  }

	  for(auto& var: vars2D){
	    if(var->IsFill()){
	      (*var->m_backgroundHists)[bkgd_ID].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
	      //Various breakdowns of selected backgrounds
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
	      (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
	      (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), weight);
	    }
	  }

	  //Only fill by tgt for non-interstitial plastic events
	  if (tgtCode != -1){
	    for(auto& var: vars_ByTgt){
	      if ((*var)[tgtCode].IsFill()){
		(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), weight);
		//Various breakdowns of selected backgrounds
		(*(*var)[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), weight);
		if ((*var)[tgtCode].IsBroken()){
		  (*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), weight);
		  (*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), weight);
		}
	      }
	    }

	    for(auto& var: vars2D_ByTgt){
	      if ((*var)[tgtCode].IsFill()){ 
		(*(*var)[tgtCode].m_backgroundHists)[bkgd_ID].FillUniverse(universe, (*var)[tgtCode].GetRecoValueX(*universe), (*var)[tgtCode].GetRecoValueY(*universe), weight);
		//Various breakdowns of selected backgrounds
		(*(*var)[tgtCode].m_BkgIntTypeHists)[intType].FillUniverse(universe, (*var)[tgtCode].GetRecoValueX(*universe), (*var)[tgtCode].GetRecoValueY(*universe), weight);
		(*(*var)[tgtCode].m_BkgTargetTypeHists)[tgtType].FillUniverse(universe, (*var)[tgtCode].GetRecoValueX(*universe), (*var)[tgtCode].GetRecoValueY(*universe), weight);
		(*(*var)[tgtCode].m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(universe, (*var)[tgtCode].GetRecoValueX(*universe), (*var)[tgtCode].GetRecoValueY(*universe), weight);
	      }
	    }
	  }
        }
      } // End band's universe loop
    } // End Band loop
  } //End entries loop
  std::cout << "Finished MC reco loop.\n";
}

void LoopAndFillData( PlotUtils::ChainWrapper* data,
			        std::vector<CVUniverse*> data_band,
				std::vector<Variable*> vars,
		                std::vector<util::Categorized<Variable, int>*> vars_ByTgt,
                                std::vector<Variable2D*> vars2D,
		                std::vector<util::Categorized<Variable2D, int>*> vars2D_ByTgt,
                                std::vector<Study*> studies,
		                PlotUtils::Cutter<CVUniverse, NeutronEvent>& michelcuts,
		                bool doNeutron = true)

{
  std::cout << "Starting data loop...\n";
  const int nEntries = data->GetEntries();
  for (int i=0; i<data->GetEntries(); ++i) {
    for (auto universe : data_band) {
      universe->SetEntry(i);
      if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::endl;
      NeutronEvent myevent = doNeutron ? NeutronEvent(universe->GetLeadNeutCandOnly()) : NeutronEvent();

      myevent.SetEMBlobInfo(universe->GetEMNBlobsTotalEnergyTotalNHits());
      std::bitset<64> SBStat = michelcuts.isDataSelected(*universe, myevent);
      myevent.SetSideBandStat(SBStat);

      if (SBStat.none()) continue;

      std::vector<double> muonMom = {universe->GetMuon4V().X(),universe->GetMuon4V().Y(),universe->GetMuon4V().Z()};

      std::vector<double> vtx = universe->GetVtx();
      double vtx_x = vtx.at(0);
      double vtx_y = vtx.at(1);
      double vtx_z = vtx.at(2);
      
      //int tgtID = util::GetRecoTargetZ(vtx_x,vtx_y,vtx_z);
      int tgtCode = util::GetRecoTargetCode(vtx_x,vtx_y,vtx_z,muonMom);

      std::cout << std::setprecision(16);

      /*
      if ((TString)(universe->ShortName()) == "cv"){
	std::cout << "Event in Data CV: " << universe->GetDouble("eventID") << std::endl;
	if (SBStat.all()) std::cout << "Selected" << std::endl;
	else std::cout << "NOT Selected" << std::endl;
	for (auto& var:vars) std::cout << "Var " << var->GetName() << ": " << var->GetRecoValue(*universe) << std::endl;
	std::cout << "" << std::endl;
      }
      */

      myevent.SetBinPTPZ(universe->GetPTPZBin());

      for(auto& study: studies) study->Selected(*universe, myevent, 1); 

      if (!SBStat.all()) continue;

      //      if (!michelcuts.isDataSelected(*universe, myevent).all()) continue;

      

      for(auto& var: vars)
      {
        if(var->IsFill()) var->dataHist->FillUniverse(universe, var->GetRecoValue(*universe), 1);
      }
      
      for(auto& var: vars2D){
        if(var->IsFill()) var->dataHist->FillUniverse(universe, var->GetRecoValueX(*universe), var->GetRecoValueY(*universe), 1);
      }

      //Only fill by tgt for non-interstitial plastic events
      if (tgtCode != -1){
	for(auto& var: vars_ByTgt){
	  if ((*var)[tgtCode].IsFill()) (*var)[tgtCode].dataHist->FillUniverse(universe, (*var)[tgtCode].GetRecoValue(*universe), 1);
	}

	for(auto& var: vars2D_ByTgt){
	  if ((*var)[tgtCode].IsFill()) (*var)[tgtCode].dataHist->FillUniverse(universe, (*var)[tgtCode].GetRecoValueX(*universe), (*var)[tgtCode].GetRecoValueY(*universe), 1);
	}
      }
    }
  }
  std::cout << "Finished data loop.\n";
}

void LoopAndFillEffDenom( PlotUtils::ChainWrapper* truth,
    				std::map<std::string, std::vector<CVUniverse*> > truth_bands,
			        std::vector<Variable*> vars,
			        std::vector<util::Categorized<Variable, int>*> vars_ByTgt,
                                std::vector<Variable2D*> vars2D,
			        std::vector<util::Categorized<Variable2D, int>*> vars2D_ByTgt,
    				PlotUtils::Cutter<CVUniverse, NeutronEvent>& michelcuts,
                                PlotUtils::Model<CVUniverse, NeutronEvent>& model)
{
  assert(!truth_bands["cv"].empty() && "\"cv\" error band is empty!  Could not set Model entry.");
  auto& cvUniv = truth_bands["cv"].front();

  std::cout << "Starting efficiency denominator loop...\n";
  const int nEntries = truth->GetEntries();
  for (int i=0; i<nEntries; ++i)
  {
    if(i%1000==0) std::cout << i << " / " << nEntries << "\r" << std::endl;

    NeutronEvent cvEvent;
    cvUniv->SetEntry(i);
    model.SetEntry(*cvUniv, cvEvent);
    const double cvWeight = model.GetWeight(*cvUniv, cvEvent);

    //=========================================
    // Systematics loop(s)
    //=========================================
    for (auto band : truth_bands)
    {
      std::vector<CVUniverse*> truth_band_universes = band.second;
      for (auto universe : truth_band_universes)
      {
        NeutronEvent myevent; //Only used to keep the Model happy

        // Tell the Event which entry in the TChain it's looking at
        universe->SetEntry(i);

        if (!michelcuts.isEfficiencyDenom(*universe, cvWeight)) continue; //Weight is ignored for isEfficiencyDenom() in all but the CV universe 
        const double weight = model.GetWeight(*universe, myevent); //Only calculate the weight for events that will use it

	std::vector<double> mc_vtx = universe->GetTrueVtx();
	double mc_vtx_x = mc_vtx.at(0);
	double mc_vtx_y = mc_vtx.at(1);
	double mc_vtx_z = mc_vtx.at(2);

	int tgtZ = universe->GetTargetZ();

	int trueTgtCode = util::GetTrueTgtCode(tgtZ, mc_vtx_x, mc_vtx_y, mc_vtx_z);

        //Fill efficiency denominator now: 
        for(auto var: vars)
        {
          if (var->IsAnaVar() && var->IsFill()) var->efficiencyDenominator->FillUniverse(universe, var->GetTrueValue(*universe), weight);
        }
	
	//I think I just need to check the true target material here... I'm not certain about that, but I think that's right... Need to sort out the right way to check that as well...
	for (auto var: vars_ByTgt){
	  if ((*var)[trueTgtCode].IsAnaVar() && (*var)[trueTgtCode].IsFill()) (*var)[trueTgtCode].efficiencyDenominator->FillUniverse(universe, (*var)[trueTgtCode].GetTrueValue(*universe), weight);
	}
	
        for(auto var: vars2D)
        {
          if(var->IsFill()) var->efficiencyDenominator->FillUniverse(universe, var->GetTrueValueX(*universe), var->GetTrueValueY(*universe), weight);
        }

	//I think I just need to check the true target material here... I'm not certain about that, but I think that's right...
	for (auto var: vars2D_ByTgt){
          if((*var)[trueTgtCode].IsFill()) (*var)[trueTgtCode].efficiencyDenominator->FillUniverse(universe, (*var)[trueTgtCode].GetTrueValueX(*universe), (*var)[trueTgtCode].GetTrueValueY(*universe), weight);
	}

      }
    }
  }
  std::cout << "Finished efficiency denominator loop.\n";
}

//Returns false if recoTreeName could not be inferred
bool inferRecoTreeNameAndCheckTreeNames(const std::string& mcPlaylistName, const std::string& dataPlaylistName, std::string& recoTreeName)
{
  const std::vector<std::string> knownTreeNames = {"Truth", "Meta"};
  bool areFilesOK = false;

  std::ifstream playlist(mcPlaylistName);
  std::string firstFile = "";
  playlist >> firstFile;
  auto testFile = TFile::Open(firstFile.c_str());
  if(!testFile)
  {
    std::cerr << "Failed to open the first MC file at " << firstFile << "\n";
    return false;
  }

  //Does the MC playlist have the Truth tree?  This is needed for the efficiency denominator.
  const auto truthTree = testFile->Get("Truth");
  if(truthTree == nullptr || !truthTree->IsA()->InheritsFrom(TClass::GetClass("TTree")))
  {
    std::cerr << "Could not find the \"Truth\" tree in MC file named " << firstFile << "\n";
    return false;
  }

  //Figure out what the reco tree name is
  for(auto key: *testFile->GetListOfKeys())
  {
    if(static_cast<TKey*>(key)->ReadObj()->IsA()->InheritsFrom(TClass::GetClass("TTree"))
       && std::find(knownTreeNames.begin(), knownTreeNames.end(), key->GetName()) == knownTreeNames.end())
    {
      recoTreeName = key->GetName();
      areFilesOK = true;
    }
  }
  delete testFile;
  testFile = nullptr;

  //Make sure the data playlist's first file has the same reco tree
  playlist.open(dataPlaylistName);
  playlist >> firstFile;
  testFile = TFile::Open(firstFile.c_str());
  if(!testFile)
  {
    std::cerr << "Failed to open the first data file at " << firstFile << "\n";
    return false;
  }

  const auto recoTree = testFile->Get(recoTreeName.c_str());
  if(recoTree == nullptr || !recoTree->IsA()->InheritsFrom(TClass::GetClass("TTree")))
  {
    std::cerr << "Could not find the \"" << recoTreeName << "\" tree in data file named " << firstFile << "\n";
    return false;
  }

  return areFilesOK;
}

//==============================================================================
// Main
//==============================================================================
int main(const int argc, const char** argv)
{
  TH1::AddDirectory(false);

  //Validate input.
  const int nArgsMandatory = 7;
  const int nArgsOptional = 2;
  const int nArgsTotal = nArgsMandatory + nArgsOptional;
  if(argc < nArgsMandatory + 1 || argc > nArgsTotal + 1) //argc is the size of argv.  I check for number of arguments + 1 because
                                //argv[0] is always the path to the executable.
  {
    std::cerr << "Expected at least " << nArgsMandatory << " and no more than " << nArgsTotal << " arguments, but got " << argc - 1 << "\n" << USAGE << "\n";
    return badCmdLine;
  }

  //One playlist must contain only MC files, and the other must contain only data files.
  //Only checking the first file in each playlist because opening each file an extra time
  //remotely (e.g. through xrootd) can get expensive.
  const std::string FVregion = argv[5],
                    mc_file_list = argv[2],
                    data_file_list = argv[1];

  const int skipSystInt = atoi(argv[3]);

  const TString tuneVer = (TString)(argv[4]);

  const int TgtNum = atoi(argv[6]);

  TString FVregionNameTmp = (TString)FVregion;

  if (FVregionNameTmp == "Targets" && TgtNum != -1){
    FVregionNameTmp = ("SingleTarget_"+util::TgtList[TgtNum]).c_str();
  }
  
  const TString FVregionName = (TString)FVregionNameTmp;

  const bool doVtx = (atoi(argv[7]) != 0);

  TString nameExt = ".root";

  bool doNeutronCuts = true;
  bool splitRecoil = false;
  double neutKESig = 10.0;

  if (argc > nArgsMandatory + 1){
    if ((TString)(argv[nArgsMandatory+1]) != "" && (TString)(argv[nArgsMandatory+1]) != "NONE") nameExt = "_"+(TString)(argv[nArgsMandatory+1])+nameExt;
    if (argc == nArgsTotal + 1){
      neutKESig = atof(argv[nArgsTotal]);
      doNeutronCuts = (neutKESig > 0);
    }
  }

  if (doNeutronCuts) nameExt = "_wNeutCuts_neutKE_"+std::to_string(neutKESig)+nameExt;
  else splitRecoil = true;
 
  if (tuneVer != "1" && tuneVer != "2"){
    std::cerr << "Must choose between 1 and 2 for the <MnvTune_v> argument. Check usage printed below. \n" << USAGE << "\n";
    return badCmdLine;
  }

  if (FVregionName != "Tracker" && !FVregionName.Contains("Target")){
    std::cerr << "<FV> argument invalid. Check usage printed below. \n" << USAGE << "\n";
    return badCmdLine;
  }
  if (FVregionName.Contains("SingleTarget") && (TgtNum < 1 || TgtNum > 6)){
    std::cerr << "<TgtNum> argument invalid. Check usage printed below. \n" << USAGE << "\n";
    return badCmdLine;
  }
  if (FVregionName == "Tracker" && TgtNum != -1){
    std::cerr << "<TgtNum> argument invalid. Check usage printed below. \n" << USAGE << "\n";
    return badCmdLine;
  }

  nameExt = "_MnvTuneV"+tuneVer+"_FVregion_"+FVregionName+nameExt;

  //Check that necessary TTrees exist in the first file of mc_file_list and data_file_list
  std::string reco_tree_name;
  if(!inferRecoTreeNameAndCheckTreeNames(mc_file_list, data_file_list, reco_tree_name))
  {
    std::cerr << "Failed to find required trees in MC playlist " << mc_file_list << " and/or data playlist " << data_file_list << ".\n" << USAGE << "\n";
    return badInputFile;
  }

  //const bool doCCQENuValidation = (reco_tree_name == "CCQENu"); //Enables extra histograms and might influence which systematics I use.
  const bool doCCQENuValidation = false;

  //const bool is_grid = false;
  PlotUtils::MacroUtil options(reco_tree_name, mc_file_list, data_file_list, "minervame1A", true);
  options.m_plist_string = util::GetPlaylist(*options.m_mc, true); //TODO: Put GetPlaylist into PlotUtils::MacroUtil

  // You're required to make some decisions
  PlotUtils::MinervaUniverse::SetNuEConstraint(true);
  PlotUtils::MinervaUniverse::SetPlaylist(options.m_plist_string); //TODO: Infer this from the files somehow?
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(-14);
  PlotUtils::MinervaUniverse::SetNFluxUniverses(100);
  PlotUtils::MinervaUniverse::SetZExpansionFaReweight(false);

  //Now that we've defined what a cross section is, decide which sample and model
  //we're extracting a cross section for.
  PlotUtils::Cutter<CVUniverse, NeutronEvent>::reco_t sidebands, preCuts;
  PlotUtils::Cutter<CVUniverse, NeutronEvent>::truth_t signalDefinition, phaseSpace;

  double minZtmp=-1, maxZtmp=-1;
  if (FVregionName == "Tracker"){
    minZtmp = 5980;
    maxZtmp = 8422;
  }
  else if (FVregionName == "Targets"){
    minZtmp = 4200;
    maxZtmp = 5980;
  }
  else if (FVregionName.Contains("Single")){
    minZtmp = (TgtNum == 1) ? 4200 : util::GetLoZ(TgtNum);
    maxZtmp = (TgtNum == 5) ? 5980 : util::GetHiZ(TgtNum);
  }
  else{
    std::cerr << "<FV> argument invalid, but passed initial check. Check usage printed below, and try and understand why the checks aren't consistent. \n" << USAGE << "\n";
    return badCmdLine;
  }

  const double minZ = minZtmp, maxZ = maxZtmp, apothem = 850; //All in mm
  preCuts.emplace_back(new reco::ZRange<CVUniverse, NeutronEvent>(FVregion, minZ, maxZ));
  preCuts.emplace_back(new reco::Apothem<CVUniverse, NeutronEvent>(apothem));
  preCuts.emplace_back(new reco::MaxMuonAngle<CVUniverse, NeutronEvent>(17.0));
  //preCuts.emplace_back(new reco::MaxMuonAngle<CVUniverse, NeutronEvent>(20.0));
  preCuts.emplace_back(new reco::HasMINOSMatch<CVUniverse, NeutronEvent>());
  preCuts.emplace_back(new reco::NoDeadtime<CVUniverse, NeutronEvent>(1, "Deadtime"));
  preCuts.emplace_back(new MyCCQECuts::PMuRange<CVUniverse, NeutronEvent>("1.5 <= Pmu <= 20",1.5,20.0));
  preCuts.emplace_back(new MyCCQECuts::IsAntiNu<CVUniverse, NeutronEvent>());
  preCuts.emplace_back(new MyCCQECuts::IsSingleTrack<CVUniverse, NeutronEvent>());
  //preCuts.emplace_back(new MyCCQECuts::LooseRecoilCut<CVUniverse, NeutronEvent>());
  //preCuts.emplace_back(new MyCCQECuts::RecoilCut<CVUniverse, NeutronEvent>());
  if (doNeutronCuts){
    preCuts.emplace_back(new MyNeutCuts::LeadNeutIs3D<CVUniverse, NeutronEvent>());
    preCuts.emplace_back(new MyNeutCuts::LeadNeutIsFarFromMuon<CVUniverse, NeutronEvent>());
    preCuts.emplace_back(new MyNeutCuts::LeadNeutZDistMin<CVUniverse, NeutronEvent>());
  }
  //preCuts.emplace_back(new MyNeutCuts::LeadNeutInTracker<CVUniverse, NeutronEvent>(maxZ));
  //preCuts.emplace_back(new reco::IsNeutrino<CVUniverse, NeutronEvent>());

  sidebands.emplace_back(new MyCCQECuts::RecoilCut<CVUniverse, NeutronEvent>());
  //sidebands.emplace_back(new MyCCQECuts::AllEMBlobsCuts<CVUniverse, NeutronEvent>(true));
  //sidebands.emplace_back(new MyCCQECuts::NoMichels<CVUniverse, NeutronEvent>());
  //sidebands.emplace_back(new MyCCQECuts::EMNBlobsCut<CVUniverse, NeutronEvent>(true));
  //sidebands.emplace_back(new MyCCQECuts::EMBlobENHitRatioCut<CVUniverse, NeutronEvent>(true));
  
  //signalDefinition.emplace_back(new truth::IsNeutrino<CVUniverse>());
  signalDefinition.emplace_back(new MySignal::IsAntiNu<CVUniverse>());
  signalDefinition.emplace_back(new truth::IsCC<CVUniverse>());
  signalDefinition.emplace_back(new MySignal::IsCorrectFS<CVUniverse>(doNeutronCuts,neutKESig));

  phaseSpace.emplace_back(new truth::ZRange<CVUniverse>(FVregion, minZ, maxZ));
  phaseSpace.emplace_back(new truth::Apothem<CVUniverse>(apothem));
  //phaseSpace.emplace_back(new truth::MuonAngle<CVUniverse>(20.0));
  phaseSpace.emplace_back(new truth::MuonAngle<CVUniverse>(17.0));
  phaseSpace.emplace_back(new MySignal::TrueMuonPRange<CVUniverse>(1.5,20.));

  PlotUtils::Cutter<CVUniverse, NeutronEvent> mycuts(std::move(preCuts), std::move(sidebands) , std::move(signalDefinition),std::move(phaseSpace));

  std::vector<std::unique_ptr<PlotUtils::Reweighter<CVUniverse, NeutronEvent>>> MnvTune;
  MnvTune.emplace_back(new PlotUtils::FluxAndCVReweighter<CVUniverse, NeutronEvent>());
  MnvTune.emplace_back(new PlotUtils::GENIEReweighter<CVUniverse, NeutronEvent>(true, false));
  MnvTune.emplace_back(new PlotUtils::LowRecoil2p2hReweighter<CVUniverse, NeutronEvent>());
  MnvTune.emplace_back(new PlotUtils::MINOSEfficiencyReweighter<CVUniverse, NeutronEvent>());
  MnvTune.emplace_back(new PlotUtils::RPAReweighter<CVUniverse, NeutronEvent>());
  if (tuneVer == "2") MnvTune.emplace_back(new PlotUtils::LowQ2PiReweighter<CVUniverse, NeutronEvent>("JOINT"));
  MnvTune.emplace_back(new PlotUtils::GeantNeutronCVReweighter<CVUniverse, NeutronEvent>());

  PlotUtils::Model<CVUniverse, NeutronEvent> model(std::move(MnvTune));

  // Make a map of systematic universes
  // Leave out systematics when making validation histograms
  const bool doSystematics = (skipSystInt == 0);
  if(!doSystematics){
    nameExt = "_SkippedSyst"+nameExt;
    std::cout << "Skipping systematics (except 1 flux universe) because <systematics> argument is non-zero.\n";
    PlotUtils::MinervaUniverse::SetNFluxUniverses(2); //Necessary to get Flux integral later...  Doesn't work with just 1 flux universe though because _that_ triggers "spread errors".
  }

  //Line to help check the new naming control. Remove once tested.
  std::cout << nameExt << std::endl;

  std::map< std::string, std::vector<CVUniverse*> > error_bands;
  if(doSystematics) error_bands = GetStandardSystematics(options.m_mc);
  else{
    std::map<std::string, std::vector<CVUniverse*> > band_flux = PlotUtils::GetFluxSystematicsMap<CVUniverse>(options.m_mc, CVUniverse::GetNFluxUniverses());
    error_bands.insert(band_flux.begin(), band_flux.end()); //Necessary to get flux integral later...
  }
  error_bands["cv"] = {new CVUniverse(options.m_mc)};
  std::map< std::string, std::vector<CVUniverse*> > truth_bands;
  if(doSystematics) truth_bands = GetStandardSystematics(options.m_truth);
  truth_bands["cv"] = {new CVUniverse(options.m_truth)};

  double radianCorr = TMath::Pi()/180.;

  //Same as Amit's seemingly. Bin normalized these are smoother.
  std::vector<double> dansPTBins = {0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5},
                      //dansPTBins = {0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1, 1.25, 1.5, 2.5, 4.5},
                      myPTBins = {0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.625, 0.7, 0.775, 0.85, 0.925, 1, 1.125, 1.25, 1.5, 2.5, 4.5},
		      myQ2QEBins = {0,0.00625,0.0125,0.025,0.0375,0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1.0,1.2,2.0},
		      //dansPzBins = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15, 20, 40, 60},
		      dansPzBins = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 15},
                      robsEmuBins = {0,1,2,3,4,5,7,9,12,15,18,22,36,50,75,100,120},
		      tejinPmuBins = {1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10, 15, 20}, //include 0, can probably handle with plotting...
		      reactionBins = {-180*radianCorr,-125*radianCorr,-80*radianCorr,-55*radianCorr,-40*radianCorr,-30*radianCorr,-20*radianCorr,-10*radianCorr,0,10*radianCorr,20*radianCorr,30*radianCorr,40*radianCorr,55*radianCorr,80*radianCorr,125*radianCorr,180*radianCorr},
		      robsRecoilBins,
		      nBlobsBins,
		      n5Bins,
		      myRecoilBins,
		      myRecoilQ2Bins,
		      myPmuBins,
		      myVtxZBins,
		      myNeutAngleBins//,
		      //myBlobEBins
		      ;

  const double robsRecoilBinWidth = 50; //MeV
  for(int whichBin = 0; whichBin < 100 + 1; ++whichBin) robsRecoilBins.push_back(robsRecoilBinWidth * whichBin);

  const double nBlobsBinWidth = 1;
  for(int whichBin = 0; whichBin < 21; ++whichBin) nBlobsBins.push_back(nBlobsBinWidth * whichBin);

  const double n5BinWidth = 1;
  for(int whichBin = 0; whichBin < 6; ++whichBin) n5Bins.push_back(n5BinWidth * whichBin);

  const double myRecoilBinWidth = 1.0/50.;
  for(int whichBin = 0; whichBin < 51; ++whichBin) myRecoilBins.push_back(myRecoilBinWidth * whichBin);

  const double myRecoilQ2BinWidth = 1.0;
  for(int whichBin = 0; whichBin < 15*50+1; ++whichBin) myRecoilQ2Bins.push_back(myRecoilQ2BinWidth * whichBin);

  const double myPmuBinWidth = 0.5;
  for(int whichBin = 0; whichBin < 41; ++whichBin) myPmuBins.push_back(myPmuBinWidth * whichBin);

  const double nBinsNeutAngle = 180./15.;
  for(int whichBin = 0; whichBin <= nBinsNeutAngle+1;++whichBin){
    myNeutAngleBins.push_back((double)(whichBin)*(15.)*radianCorr);
  }

  const double myVtxZBinWidth = 1.;
  //const double myVtxZBase = 5800.; //Tracker for plot testing
  const double myVtxZBase = minZ; //Targets for later!!!
  const int nVtxZBins = ceil((maxZ-minZ)/myVtxZBinWidth);
  for(int whichBin = 0; whichBin < nVtxZBins; ++whichBin) myVtxZBins.push_back(myVtxZBinWidth * whichBin + myVtxZBase);

  //const double myBlobEBinWdith = 3.;
  //for(int whichBin = 0; whichBin < 51; ++whichBin) myBlobEBins.push_back(myBlobEBinWidth * whichBin);

  std::vector<Variable*> vars = {
    new Variable(true, "pTmu", "p_{T, #mu} [GeV/c]", dansPTBins, &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue),
    //new Variable(false, "pzmu", "p_{||, #mu} [GeV/c]", dansPzBins, &CVUniverse::GetMuonPz, &CVUniverse::GetMuonPzTrue),
    new Variable(true, "neutAngleToMuon", "Neutron Angle To Muon[radian]", myNeutAngleBins, &CVUniverse::GetLeadNeutAngleToMuon, &CVUniverse::GetMaxFSNeutronAngleToMuon),
    new Variable(true, "neutAngle", "#theta_{n}[radian]", myNeutAngleBins, &CVUniverse::GetLeadNeutAngleToBeam, &CVUniverse::GetMaxFSNeutronAngleToBeam),
    new Variable(true, "neutDeltaPhiT", "#delta#phi_{T}[radian]", myNeutAngleBins, &CVUniverse::GetLeadNeutDeltaPhiT, &CVUniverse::GetMaxFSNeutronDeltaPhiT),
    new Variable(true, "neutReactionAngle", "#delta#theta_{r}[radian]", reactionBins, &CVUniverse::GetLeadNeutReactionAngle, &CVUniverse::GetMaxFSNeutronReactionAngle),
    new Variable(true, "neutCoplanarityAngle", "#delta#theta_{p}[radian]", reactionBins, &CVUniverse::GetLeadNeutCoplanarityAngle, &CVUniverse::GetMaxFSNeutronCoplanarityAngle),
    //new Variable((TString)("MyBins"),"pTmu_MYBins", "p_{T, #mu} [GeV/c]", myPTBins, &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue),
    //new Variable(false, "pTmu_MYBins", "p_{T, #mu} [GeV/c]", myPTBins, &CVUniverse::GetMuonPT, &CVUniverse::GetMuonPTTrue),
    //new Variable(false, "nBlobs", "No.", nBlobsBins, &CVUniverse::GetNNeutBlobs),//Don't need GetDummyTrue perhaps...
    //new Variable(false, "recoilE", "Recoil E [GeV]", myRecoilBins, &CVUniverse::GetDANRecoilEnergyGeV),//Don't need GetDummyTrue perhaps...
    //new Variable(false, "Q2QE", "Q^{2}_{QE} [GeV^{2}]", myQ2QEBins, &CVUniverse::GetQ2QEPickledGeV),
    //new Variable("nMichel","No.", n5Bins, &CVUniverse::GetNImprovedMichel),
    //new Variable("nTrack","No.", n5Bins, &CVUniverse::GetNTracks),
    //new Variable("pmu", "p_{#mu} [GeV/c]", myPmuBins, &CVUniverse::GetMuonP, &CVUniverse::GetMuonPTrue),//Don't need GetDummyTrue perhaps...
    //new Variable(false,"vtxZ", "Z [mm]", myVtxZBins, &CVUniverse::GetVtxZ, &CVUniverse::GetTrueVtxZ),//Don't need GetDummyTrue perhaps...
    //new Variable("recQ2Bin","No.",myRecoilQ2Bins, &CVUniverse::GetRecoilQ2Bin),
  };

  /*
  std::map<int, std::string> TgtList = {{1,"Tgt1"},
					{2,"Tgt2"},
					{3,"Tgt3"},
					{4,"Tgt4"},
					{5,"Tgt5"},
					{6,"WaterTgt"},
					{10,"Plastic_US1"},
					{21,"Plastic_US2_DS1"},
					{32,"Plastic_US3_DS2"},
					{63,"Plastic_USWater_DS3"},
					{46,"Plastic_US4_DSWater"},
					{54,"Plastic_US5_DS4"},
					{0,"Plastic_DS5"}};
  */

  std::vector<util::Categorized<Variable, int>*> vars_ByTgt = {};
  if (FVregionName.Contains("Target")){
    for(auto& var: vars){ 
      var->SetFillVar(false);
      if (!doNeutronCuts) continue;
      TString nameCheck = var->GetName();
      if (nameCheck == "vtxZ") continue;
      vars_ByTgt.push_back(new util::Categorized<Variable, int>(var->GetDirectoryName(), "ByTgt", var->IsAnaVar(), var->GetName().c_str(),var->GetAxisLabel().c_str(),util::TgtCodeList[TgtNum],var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
    }
  }

  std::vector<Variable2D*> vars2D = {};
  /*
    new Variable2D(true,"pmu2D",*vars[1],*vars[0]),//pT v. p"z"
    //new Variable2D(false, *vars[5],*vars[4]),//recoil v. Q2
    //new Variable2D(false, *vars[vars.size()-1],*vars[0]),//pT v. recoilQ2Bin
    };*/
  std::vector<util::Categorized<Variable2D, int>*> vars2D_ByTgt = {};

  if (!doNeutronCuts){
    if (FVregionName.Contains("Target")){
      if (!doVtx){
	vars2D_ByTgt.push_back(new util::Categorized<Variable2D, int>("", "ByTgt", true, "pmu2D", util::TgtCodeList[TgtNum], *vars[1], *vars[0]));
      }
    }
    else{
      for (auto& var: vars) var->SetFillVar(false);
      vars2D.push_back(new Variable2D(true,"pmu2D",*vars[1],*vars[0]));
    }
  }
  
  //Commented out during testing of analysis variable flag
  /*
  if(doCCQENuValidation)
  {
  std::cerr << "Detected that tree name is CCQENu.  Making validation histograms.\n";
    vars.push_back(new Variable("Emu", "E_{#mu} [GeV]", robsEmuBins, &CVUniverse::GetEmuGeV, &CVUniverse::GetElepTrueGeV));
    vars.push_back(new Variable("Erecoil", "E_{recoil}", robsRecoilBins, &CVUniverse::GetRecoilE, &CVUniverse::Getq0True)); //TODO: q0 is not the same as recoil energy without a spline correction
    vars2D.push_back(new Variable2D(*vars[1], *vars[0]));
  }
  */
    
  CVUniverse* data_universe = new CVUniverse(options.m_data);
  std::vector<CVUniverse*> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse*> > data_error_bands;
  data_error_bands["cv"] = data_band;
  
  std::vector<Study*> data_studies;

  std::vector<Study*> studies = {
    //new EMSideBands(vars, error_bands, truth_bands, data_band),
    //new MichelAndNBlobSB(vars, error_bands, truth_bands, data_band),
    //new NeutronVariables(maxZ, minZ, error_bands, truth_bands, data_band),
    //new RecoilSB(vars, error_bands, truth_bands, data_band, splitRecoil),
    //new LeadNeutStudy(error_bands, truth_bands, data_band),
    new PreRecoil(vars, error_bands, truth_bands, data_band, splitRecoil, doNeutronCuts, FVregionName, TgtNum, doVtx),
  };

  for(auto& var: vars) if(var->IsFill()) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: vars) if(var->IsFill()) var->InitializeDATAHists(data_band);

  for(auto& tgt: vars_ByTgt){ 
    tgt->visit([&error_bands, &truth_bands](Variable& var)
	       {
		 if (var.IsFill()) var.InitializeMCHists(error_bands, truth_bands);
	       });
  }
  for(auto& tgt: vars_ByTgt){ 
    tgt->visit([&data_band](Variable& var)
	       {
		 if (var.IsFill()) var.InitializeDATAHists(data_band);
	       });
  }

  for(auto& var: vars2D) if(var->IsFill()) var->InitializeMCHists(error_bands, truth_bands);
  for(auto& var: vars2D) if(var->IsFill()) var->InitializeDATAHists(data_band);

  for(auto& tgt: vars2D_ByTgt){ 
    tgt->visit([&error_bands, &truth_bands](Variable2D& var)
	       {
		 if (var.IsFill()) var.InitializeMCHists(error_bands, truth_bands);
	       });
  }
  for(auto& tgt: vars2D_ByTgt){ 
    tgt->visit([&data_band](Variable2D& var)
	       {
		 if (var.IsFill()) var.InitializeDATAHists(data_band);
	       });
  }

  // Loop entries and fill
  try
  {
    CVUniverse::SetTruth(false);
    LoopAndFillEventSelection(options.m_mc, error_bands, vars, vars_ByTgt, vars2D, vars2D_ByTgt, studies, mycuts, model, doNeutronCuts);
    //LoopAndFillEventSelection(options.m_mc, error_bands, vars, vars_ByTgt, vars2D, vars2D_ByTgt, studies, mycuts, model);
    CVUniverse::SetTruth(true);
    LoopAndFillEffDenom(options.m_truth, truth_bands, vars, vars_ByTgt, vars2D, vars2D_ByTgt, mycuts, model);
    options.PrintMacroConfiguration(argv[0]);
    std::cout << "MC cut summary:\n" << mycuts << "\n";
    mycuts.resetStats();

    CVUniverse::SetTruth(false);
    LoopAndFillData(options.m_data, data_band, vars, vars_ByTgt, vars2D, vars2D_ByTgt, studies, mycuts, doNeutronCuts);
    //LoopAndFillData(options.m_data, data_band, vars, vars_ByTgt, vars2D, vars2D_ByTgt, studies, mycuts);
    std::cout << "Data cut summary:\n" << mycuts << "\n";

    std::cout << "Writing MC Output File." << std::endl;

    //Write MC results
    TFile* mcOutDir = TFile::Open((TString)(MC_OUT_FILE_NAME)+nameExt, "RECREATE");
    if(!mcOutDir)
    {
      std::cerr << "Failed to open a file named " << MC_OUT_FILE_NAME << nameExt << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    std::cout << "Actually Setting Directories in the File" << std::endl;

    for(auto& study: studies) study->SaveOrDrawMC(*mcOutDir);
    for(auto& var: vars) if(var->IsFill()) var->WriteMC(*mcOutDir);
    for(auto& tgt: vars_ByTgt){ 
      tgt->visit([mcOutDir](Variable& var)
		 {
		   if (var.IsFill()) var.WriteMC(*mcOutDir);
		 });
    }

    for(auto& var: vars2D) if(var->IsFill()) var->WriteMC(*mcOutDir);
    for(auto& tgt: vars2D_ByTgt){ 
      tgt->visit([mcOutDir](Variable2D& var)
		 {
		   if (var.IsFill()) var.WriteMC(*mcOutDir);
		 });
    }

    std::cout << "Writing" << std::endl;

    mcOutDir->Write();

    //Protons On Target
    auto mcPOT = new TParameter<double>("POTUsed", options.m_mc_pot);
    mcPOT->Write();

    PlotUtils::TargetUtils targetInfo;
    assert(error_bands["cv"].size() == 1 && "List of error bands must contain a universe named \"cv\" for the flux integral.");

    //Removed for antineutrino and targets since these functions don't cover those cases as written yet.

    /*
    for(const auto& var: vars)
    {
      //Flux integral only if systematics are being done (temporary solution)
      util::GetFluxIntegral(*error_bands["cv"].front(), var->efficiencyNumerator->hist)->Write((var->GetName() + "_reweightedflux_integrated").c_str());
      //Always use MC number of nucleons for cross section
      auto nNucleons = new TParameter<double>((var->GetName() + "_fiducial_nucleons").c_str(), targetInfo.GetTrackerNNucleons(minZ, maxZ, true, apothem));
      nNucleons->Write();
    }
    */

    std::cout << "Writing Data Output File" << std::endl;

    //Write data results
    TFile* dataOutDir = TFile::Open((TString)(DATA_OUT_FILE_NAME)+nameExt, "RECREATE");
    if(!dataOutDir)
    {
      std::cerr << "Failed to open a file named " << DATA_OUT_FILE_NAME << nameExt << " in the current directory for writing histograms.\n";
      return badOutputFile;
    }

    for(auto& study: studies) study->SaveOrDrawData(*dataOutDir);
    for(auto& var: vars) if(var->IsFill()) var->WriteData(*dataOutDir);
    for(auto& tgt: vars_ByTgt){
      tgt->visit([dataOutDir](Variable& var)
		 {
		   if (var.IsFill()) var.WriteData(*dataOutDir);
		 });
    }
    for(auto& var: vars2D) if(var->IsFill()) var->WriteData(*dataOutDir);
    for(auto& tgt: vars2D_ByTgt){ 
      tgt->visit([dataOutDir](Variable2D& var)
		 {
		   if (var.IsFill()) var.WriteData(*dataOutDir);
		 });
    }

    dataOutDir->Write();

    //Protons On Target
    auto dataPOT = new TParameter<double>("POTUsed", options.m_data_pot);
    dataPOT->Write();

    std::cout << "Success" << std::endl;
  }
  catch(const ROOT::exception& e)
  {
    std::cerr << "Ending on a ROOT error message.  No histograms will be produced.\n"
              << "If the message talks about \"TNetXNGFile\", this could be a problem with dCache.  The message is:\n"
              << e.what() << "\n" << USAGE << "\n";
    return badFileRead;
  }

  return success;
}
