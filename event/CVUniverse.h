// =============================================================================
// David Last adaptation of the CVUniverse for his analysis.
// Contact email: dlast@sas.upenn.edu (lastd44@gmail.com)
//
// Base class for an un-systematically shifted (i.e. CV) universe. Implement
// "Get" functions for all the quantities that you need for your analysis.
//
// This class inherits from PU::MinervaUniverse, which in turn inherits from
// PU::BaseUniverse. PU::BU defines the interface with anatuples.
// 
// Within the class, "WeightFunctions" and "MuonFunctions" are included to gain
// access to standardized weight and muon variable getters. See:
// https://cdcvs.fnal.gov/redmine/projects/minerva-sw/wiki/MinervaUniverse_Structure_
// for a full list of standardized functions you can use. In general, if a
// standard version of a function is available, you should be using it.
// =============================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

#include <iostream>

#include "PlotUtils/PhysicsVariables.h"//Included by David, unsure if needed
#include "PlotUtils/MinervaUniverse.h"
//Needed for neutron candidates business... May change at some point, but for now this is what we're working with.
#include "event/NeutCands.h"
#include "TVector3.h"

class CVUniverse : public PlotUtils::MinervaUniverse {

  public:
  #include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
  #include "PlotUtils/TruthFunctions.h" //Getq3True
  #include "PlotUtils/RecoilEnergyFunctions.h" //GetRecoilEnergy
  // ========================================================================
  // Constructor/Destructor
  // ========================================================================
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0)
      : PlotUtils::MinervaUniverse(chw, nsigma) {}

  virtual ~CVUniverse() {}

  // ========================================================================
  // Quantities defined here as constants for the sake of below. Definition
  // matched to Dan's CCQENuInclusiveME variables from:
  // `/minerva/app/users/drut1186/cmtuser/Minerva_v22r1p1_OrigCCQENuInc/Ana/CCQENu/ana_common/include/CCQENuUtils.h`
  // ========================================================================
  static constexpr double M_n = 939.56536;
  static constexpr double M_p = 938.272013;
  static constexpr double M_nucleon = (1.5*M_n+M_p)/2.5;

  static constexpr int PDG_n = 2112;
  static constexpr int PDG_p = 2212;

  // ========================================================================
  // Write a "Get" function for all quantities access by your analysis.
  // For composite quantities (e.g. Enu) use a calculator function.
  //
  // In order to properly calculate muon variables and systematics use the
  // various functions defined in MinervaUniverse.
  // E.g. GetPmu, GetEmu, etc.
  // ========================================================================

  // Quantities only needed for cuts
  // Although unlikely, in principle these quanties could be shifted by a
  // systematic. And when they are, they'll only be shifted correctly if we
  // write these accessor functions.

  //Muon kinematics
  double GetMuonP() const //GeV/c
  {
    return GetPmu()/1000.0 ;
  }

  double GetMuonPT() const //GeV/c
  {
    return GetPmu()/1000. * sin(GetThetamu());
  }

  double GetMuonPz() const //GeV/c
  {
    return GetPmu()/1000. * cos(GetThetamu());
  }

  double GetMuonPTTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * sin(GetThetalepTrue());
  }

  double GetMuonPzTrue() const //GeV/c
  {
    return GetPlepTrue()/1000. * cos(GetThetalepTrue());
  }

  double GetMuonPTrue() const //GeV/c
  {
    return GetPlepTrue()/1000.;
  }

  double GetEmuGeV() const //GeV
  {
    return GetEmu()/1000.;
  }

  double GetElepTrueGeV() const //GeV
  {
    return GetElepTrue()/1000.;
  }

  int GetInteractionType() const {
    return GetInt("mc_intType");
  }

  int GetTargetNucleon() const {
    return GetInt("mc_targetNucleon");
  }
  
  double GetBjorkenXTrue() const {
    return GetDouble("mc_Bjorkenx");
  }

  double GetBjorkenYTrue() const {
    return GetDouble("mc_Bjorkeny");
  }

  virtual bool IsMinosMatchMuon() const {
    return GetInt("has_interaction_vertex") == 1;
  }
  
  ROOT::Math::XYZTVector GetVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("vtx").data());
    return result;
  }

  ROOT::Math::XYZTVector GetTrueVertex() const
  {
    ROOT::Math::XYZTVector result;
    result.SetCoordinates(GetVec<double>("mc_vtx").data());
    return result;
  }

  virtual int GetTDead() const {
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj");
  }
  
  //TODO: If there was a spline correcting Eavail, it might not really be Eavail.
  //      Our energy correction spline, one of at least 2 I know of, corrects q0
  //      so that we get the right neutrino energy in an inclusive sample.  So,
  //      this function could be correcting for neutron energy which Eavail should
  //      not do.
  virtual double GetEavail() const {
    return GetDouble("recoilE_SplineCorrected");
  }
  
  virtual double GetQ2Reco() const{
    return GetDouble("qsquared_recoil");
  }

  //GetRecoilE is designed to match the NSF validation suite
  virtual double GetRecoilE() const {
    return GetVecElem("recoil_summed_energy", 0);
  }
  
  virtual double Getq3() const{
    double eavail = GetEavail()/pow(10,3);
    double q2 = GetQ2Reco() / pow(10,6);
    double q3mec = sqrt(eavail*eavail + q2);
    return q3mec;
  }
   
  virtual int GetCurrent() const { return GetInt("mc_current"); }

  virtual int GetTruthNuPDG() const { return GetInt("mc_incoming"); }

  virtual double GetMuonQP() const {
    return GetDouble((GetAnaToolName() + "_minos_trk_qp").c_str());
  }

  //Some functions to match CCQENuInclusive treatment of DIS weighting. Name matches same Dan area as before.
  virtual double GetTrueExperimentersQ2() const {
    double Enu = GetEnuTrue(); //MeV
    double Emu = GetElepTrue(); //MeV
    double thetaMu = GetThetalepTrue();
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double CalcTrueExperimentersQ2(double Enu, double Emu, double thetaMu) const{
    return 4.0*Enu*Emu*pow(sin(thetaMu/2.0),2.0);//MeV^2
  }

  virtual double GetTrueExperimentersW() const {
    double nuclMass = M_nucleon;
    int struckNucl = GetTargetNucleon();
    if (struckNucl == PDG_n){
      nuclMass=M_n;
    }
    else if (struckNucl == PDG_p){
      nuclMass=M_p;
    }
    double Enu = GetEnuTrue();
    double Emu = GetElepTrue();
    double thetaMu = GetThetalepTrue();
    double Q2 = CalcTrueExperimentersQ2(Enu, Emu, thetaMu);
    return TMath::Sqrt(pow(nuclMass,2) + 2.0*(Enu-Emu)*nuclMass - Q2);
  }

  // Functions added by David that have no match in the above.
  virtual int GetNTracks() const { return GetInt("multiplicity"); }

  virtual int GetNNeutBlobs() const { return GetInt((GetAnaToolName() + "_BlobIs3D_sz").c_str()); }

  virtual std::vector<double> GetEMBlobStartZVec() const { return GetVec<double>("nonvtx_iso_blobs_start_position_z_in_prong"); }

  virtual std::vector<int> GetEMBlobNHitsVec() const { return GetVec<int>("nonvtx_iso_blobs_n_hits_in_prong"); }

  virtual std::vector<double> GetEMBlobEnergyVec() const { return GetVec<double>("nonvtx_iso_blobs_energy_in_prong"); }

  virtual std::vector<double> GetEMNBlobsTotalEnergyTotalNHits(double shift = 0) const {
    std::vector<double> info;
    double nBlobs = 0;
    double totalE = shift;
    double nHits = 0;
    std::vector<double> StartZVec = GetEMBlobStartZVec();
    std::vector<double> EnergyVec = GetEMBlobEnergyVec();
    std::vector<int> NHitsVec = GetEMBlobNHitsVec();
    for (unsigned int i=0; i<StartZVec.size(); ++i){
      if (StartZVec.at(i) > 4750.0){
	nBlobs+=1.0;
	totalE+=EnergyVec.at(i);
	nHits+=(double)NHitsVec.at(i);
      }
    }
    info.push_back(nBlobs);
    info.push_back(totalE);
    info.push_back(nHits);
    return info;
  }

  virtual int GetTargetZ() const { return GetInt("mc_targetZ"); }

  virtual int GetNFSPart() const { return GetInt("mc_nFSPart"); }

  virtual std::vector<int> GetFSPartPDG() const { return GetVec<int>("mc_FSPartPDG"); }

  virtual std::vector<double> GetFSPartE() const { return GetVec<double>("mc_FSPartE"); }

  virtual std::vector<double> GetFSPartPx() const { return GetVec<double>("mc_FSPartPx"); }

  virtual std::vector<double> GetFSPartPy() const { return GetVec<double>("mc_FSPartPy"); }

  virtual std::vector<double> GetFSPartPz() const { return GetVec<double>("mc_FSPartPz"); }

  virtual double GetMaxFSNeutronKE() const {
    double max_KE = -999.;
    std::vector<int> PDGs = GetFSPartPDG();
    std::vector<double> Es = GetFSPartE();
    for (int iFS=0; iFS < PDGs.size(); ++iFS){
      double KE = Es.at(iFS)-M_n;
      if (PDGs.at(iFS) == 2112 && KE > max_KE) max_KE = KE;
    }
    return max_KE;
  }

  virtual int GetNImprovedMichel() const { return GetInt("improved_michel_vertex_type_sz"); }

  virtual int GetNuHelicity() const { return GetInt((GetAnaToolName()+"_nuHelicity").c_str()); }

  virtual double GetEnuCCQEPickledGeV() const{ //RETURNS IN MeV^2
    int charge=-1; //hard-coded since I'm focused on anti-nu
    double enu=PlotUtils::nuEnergyCCQE( GetEmu(), GetPmu(), GetThetamu(), charge)*MeVGeV;
    return enu;
  }

  virtual double GetQ2QEPickledGeV() const{ //RETURNS IN MeV^2
    int charge=-1; //hard-coded since I'm focused on anti-nu
    if (GetEnuCCQEPickledGeV()<=0.0) return 0.0;
    else{
      double q2=PlotUtils::qSquaredCCQE( GetEmu(), GetPmu(), GetThetamu(), charge)*MeVGeV*MeVGeV;
      return q2;
    }
  }

  // Functions added by David that have a confusing match to above, be careful with naming
  virtual int GetHasInteractionVertex() const { return GetInt("has_interaction_vertex"); }

  virtual std::vector<double> GetVtx() const { return GetVec<double>("vtx"); }

  virtual std::vector<double> GetTrueVtx() const { return GetVec<double>("mc_vtx"); }

  virtual double GetVtxZ() const { return GetVtx().at(2); }

  virtual double GetTrueVtxZ() const { return GetTrueVtx().at(2); }

  virtual int GetIsMinosMatchTrack() const { return GetInt("isMinosMatchTrack"); }
  
  virtual int GetIsMinosMatchTrackOLD() const { return GetInt("muon_is_minos_match_track"); }
  
  virtual int GetIsMinosMatchStub() const { return GetInt("isMinosMatchStub"); }
  
  virtual int GetIsMinosMatchStubOLD() const { return GetInt("muon_is_minos_match_stub"); }

  // Functions added by David that have a conflicting match above, be careful with naming
  double MeVGeV=0.001;

  virtual double GetCalRecoilEnergy() const{
    return GetDouble("recoil_energy_nonmuon_nonvtx100mm")+GetDouble("recoil_energy_nonmuon_nonvtx100mm_nuclTargs");
  /*
  if (GetVec<double>("recoil_summed_energy").size()==0) return -999.0;
    return (GetVec<double>("recoil_summed_energy")[0]-GetDouble("recoil_energy_nonmuon_vtx100mm"));
  */
  }

  virtual double GetNonCalRecoilEnergy() const{
    return 0.0;
  }

  virtual double GetDANRecoilEnergyGeV() const{
  //double recoilE = GetDouble("recoil_energy_nonmuon_nonvtx100mm")*MeVGeV;
    double recoilE = GetCalRecoilEnergy()*MeVGeV;
    return recoilE;
  }

  virtual double GetRecoilEnergyGeV() const{
    double recoilE = GetRecoilEnergy()*MeVGeV;
    return recoilE;
  }

  virtual int GetPTPZBin() const{
    double pT = GetMuonPT();
    double pZ = GetMuonPz();
    int bin = -1;
    if (pZ < 1.5 || pZ > 15.0 || pT < 0.0 || pT > 2.5) return bin;
    if (pZ <= 5.0){
      //if (pT <= 0.2) bin = 0;
      if (pT <= 0.25) bin = 0;
      else if (pT <= 0.4) bin = 1;
      //else if (pT <= 0.65) bin = 2;
      else if (pT <= 0.7) bin = 2;
      //else if (pT <= 0.82) bin = 3;
      else if (pT <= 0.85) bin = 3;
      else if (pT <= 1.0) bin = 4;
      else bin = 5;
    }
    else if (pZ <= 8.0){
      //if (pT <= 0.2) bin = 6;
      if (pT <= 0.25) bin = 6;
      else if (pT <= 0.4) bin = 7;
      //else if (pT <= 0.65) bin = 8;
      else if (pT <= 0.7) bin = 8;
      //else if (pT <= 0.82) bin = 9;
      else if (pT <= 0.85) bin = 9;
      else if (pT <= 1.0) bin = 10;
      else bin = 11;
    }
    else{
      //if (pT <= 0.5) bin = 12;
      if (pT <= 0.55) bin = 12;
      else bin = 13;
    }
    return bin;
  }

  virtual double GetRecoilQ2Bin() const{
    double bin = -0.5;
    double binSize = 1.0/50.0; //Hard-coded bin size for recoil energy

    double recE = GetDANRecoilEnergyGeV();
    double Q2 = GetQ2QEPickledGeV();
    if (recE < 0 || Q2 < 0 || recE > 1.0 || Q2 > 2.0) return bin;

    int binRec = (int)(recE/binSize);
    int binQ2 = -1;
    if (Q2 < 0.00625) binQ2 = 0;
    else if (Q2 < 0.0125) binQ2 = 1;
    else if (Q2 < 0.025) binQ2 = 2;
    else if (Q2 < 0.0375) binQ2 = 3;
    else if (Q2 < 0.05) binQ2 = 4;
    else if (Q2 < 0.1) binQ2 = 5;
    else if (Q2 < 0.15) binQ2 = 6;
    else if (Q2 < 0.2) binQ2 = 7;
    else if (Q2 < 0.3) binQ2 = 8;
    else if (Q2 < 0.4) binQ2 = 9;
    else if (Q2 < 0.6) binQ2 = 10;
    else if (Q2 < 0.8) binQ2 = 11;
    else if (Q2 < 1.0) binQ2 = 12;
    else if (Q2 < 1.2) binQ2 = 13;
    else binQ2 = 14;
    bin = binRec*1.0+binQ2*50.0+0.5;
    return bin;
  }

  //Neutron Candidate Business. May change at some point, but this is current to what was used before/validating against...
  virtual std::vector<double> GetNeutCandEs() const{ return GetVec<double>((GetAnaToolName()+"_BlobTotalE").c_str()); }

  //returns a NeutCands object of just the lead candidate in the interest of verifying that this doesn't affect cuts based around the event object. Further restructuring will likely just get the single neut cand into the event an allow for filling of one or the other. I could also more intelligently make the NeutCands obejct generally, but it's what it is for now. This is a short term solution to check that the code can be better optimized to not take so long by not filling more than one blob per event.
  virtual NeutronCandidates::NeutCands GetLeadNeutCandOnly(){
    std::vector<NeutronCandidates::NeutCand> cands = {};
    int nBlobs = GetNNeutBlobs();

    if (nBlobs > 0){
      std::vector<double> Es = GetNeutCandEs();
      int leadNeutEIndex = std::max_element(Es.begin(),Es.end()) - Es.begin();

      cands.push_back(GetNeutCand(leadNeutEIndex));
    }

    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  };

  virtual NeutronCandidates::NeutCand GetNeutCand(int index){
    std::vector<double> vtx = GetVtx();
    TVector3 EvtVtx;
    EvtVtx.SetXYZ(vtx.at(0),vtx.at(1),vtx.at(2));
    NeutronCandidates::intCandData intData;
    NeutronCandidates::doubleCandData doubleData;
    std::string toolName = GetAnaToolName();
    for (const auto& intMember: NeutronCandidates::GetBranchIntMap()){
      intData[intMember.first]={};
      for (const auto& branchName: intMember.second){
	intData[intMember.first].push_back(GetVecElemInt((toolName+branchName).c_str(),index));
      }
    }
    for (const auto& doubleMember: NeutronCandidates::GetBranchDoubleMap()){
      doubleData[doubleMember.first]={};
      for (const auto& branchName: doubleMember.second){
	doubleData[doubleMember.first].push_back(GetVecElem((toolName+branchName).c_str(),index));
      }
    }
    return NeutronCandidates::NeutCand(intData,doubleData,EvtVtx);
  };
  
  virtual NeutronCandidates::NeutCands GetNeutCands(){
    std::vector<NeutronCandidates::NeutCand> cands;
    int nBlobs = GetNNeutBlobs();
    for(int neutBlobIndex=0; neutBlobIndex < nBlobs; ++neutBlobIndex){
      cands.push_back(GetNeutCand(neutBlobIndex));
    }
    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  };

  //Still needed for some systematics to compile, but shouldn't be used for reweighting anymore.
  protected:
  #include "PlotUtils/WeightFunctions.h" // Get*Weight
};

#endif
