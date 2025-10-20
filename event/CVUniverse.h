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

#include "PlotUtils/ErrorHandler.h" //For ROOT::exception in case it's being used to react to non-existent files.

#include "PlotUtils/PhysicsVariables.h"//Included by David, unsure if needed
#include "PlotUtils/MinervaUniverse.h"
//Needed for neutron candidates business... May change at some point, but for now this is what we're working with.
#include "event/NeutCands.h"
#include "TVector3.h"
#include "TRandom3.h"

class CVUniverse : public PlotUtils::MinervaUniverse {
  public:
  #include "PlotUtils/MuonFunctions.h" // GetMinosEfficiencyWeight
  #include "PlotUtils/TruthFunctions.h" //Getq3True
  #include "PlotUtils/RecoilEnergyFunctions.h" //GetRecoilEnergy
  // ========================================================================
  // Constructor/Destructor
  // ========================================================================
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0)
    : PlotUtils::MinervaUniverse(chw, nsigma), m_LeadNeutIndex(-999) {
    m_Random = new TRandom3(0);
    LoadNeutronReweightHistos();
  }

  virtual ~CVUniverse() {}

  int m_LeadNeutIndex;

  std::map<int, TH2D*> m_NeutRWHists;
  
  std::map<int, TString> NeutRWCategs = {{11, "sigQE"},
					 {18, "sig2p2h"},
					 {19, "sigOther"},
					 {1, "1chargePi"},
					 {2, "1neutPi"},
					 {3, "NPi"},
					 {4, "SubThresh"},
					 {5, "TrackableProt"},
					 {-999, "Other"}};
  
  void LoadNeutronReweightHistos(){
    for (auto Categ : NeutRWCategs){
      m_NeutRWHists[Categ.first] = nullptr;
    }

    try{ 
      std::string weightFileName = "";
      //if(std::getenv("PLOTUTILSROOT")) weightFileName = std::string(std::getenv("PLOTUTILSROOT")) + "/../etc/extraWeightFiles/BURP.root";
      //if(std::getenv("PLOTUTILSROOT")) weightFileName = std::string(std::getenv("PLOTUTILSROOT")) + "/../etc/extraWeightFiles/test.root";
      if(std::getenv("PLOTUTILSROOT")) weightFileName = std::string(std::getenv("PLOTUTILSROOT")) + "/../etc/extraWeightFiles/NeutronRenorm_" + GetPlaylist() + ".root";
      std::unique_ptr<TFile> weightFile(TFile::Open(weightFileName.c_str()));
      
      for (auto Categ : NeutRWCategs){
	TH2D* hist = (TH2D*)(weightFile->Get("NeutronRenormWeight_"+Categ.second));
	if (hist) m_NeutRWHists[Categ.first] = hist;
      }
    }
    catch(const ROOT::warning& /*w*/){
      std::cout << "No Neutron Renormalization File. Doing without" << std::endl;
    }

  }
  
  virtual void OnNewEntry() override{
    m_LeadNeutIndex = -999;//Resetting to avoid any possible mishaps with the indexing of an array.
  }
  
  // ========================================================================
  // Quantities defined here as constants for the sake of below. Definition
  // matched to Dan's CCQENuInclusiveME variables from:
  // `/minerva/app/users/drut1186/cmtuser/Minerva_v22r1p1_OrigCCQENuInc/Ana/CCQENu/ana_common/include/CCQENuUtils.h`
  // ========================================================================
  static constexpr double M_n = 939.56536;
  static constexpr double M_p = 938.272013;
  static constexpr double M_pi = 139.57061;
  static constexpr double M_nucleon = (1.5*M_n+M_p)/2.5;

  static constexpr int PDG_n = 2112;
  static constexpr int PDG_p = 2212;
  
  static constexpr double MeVGeV=0.001;

  std::pair<double, double> ProjectTrackToZ(std::vector<double> point, std::vector<double> trkMom, double z){
    std::pair<double,double> ret = std::make_pair(-999,-999);
    if (point.size() < 3 || trkMom.size() !=3){
      return ret;
    }
    double dXPerZ = trkMom.at(0)/trkMom.at(2);
    double dYPerZ = trkMom.at(1)/trkMom.at(2);
    double dZ = z-point.at(2);
    double X = point.at(0) + dZ*dXPerZ;
    double Y = point.at(1) + dZ*dYPerZ;
    ret = std::make_pair(X,Y);
    return ret;
  }


  double GetTCutFromAngle(double angle, std::vector<double> trkMom, double dZ) const{
    double ret=0.0;
    if (trkMom.size() != 3) return ret;
    const double radianCorr = TMath::Pi()/180.;    
    double angleTan = TMath::Tan(radianCorr*angle);
    double dRdZ = TMath::Sqrt(trkMom.at(0)*trkMom.at(0) + trkMom.at(1)*trkMom.at(1) + trkMom.at(2)* trkMom.at(2))/(trkMom.at(2));
    double dR = dRdZ*dZ;
    ret = dR*angleTan;
    
    return ret;
  }
  
  double GetTXFromXY(double x, double y, int view) const{
    const double radianCorr = TMath::Pi()/180.;
    if (view == 1) return x;
    else if (view == 2) return x*TMath::Cos(radianCorr*60)-y*TMath::Sin(radianCorr*60);
    else if (view == 3) return x*TMath::Cos(radianCorr*60)+y*TMath::Sin(radianCorr*60);
    else return -999999999;
  }  
  
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
    int matchMuon = GetIsMinosMatchTrack();
    return (matchMuon == 1);
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
  
  virtual int GetTargetA() const { return GetInt("mc_targetA"); }

  virtual int GetNFSPart() const { return GetInt("mc_nFSPart"); }

  virtual std::vector<int> GetFSPartPDG() const { return GetVec<int>("mc_FSPartPDG"); }

  virtual std::vector<double> GetFSPartE() const { return GetVec<double>("mc_FSPartE"); }

  virtual std::vector<double> GetFSPartPx() const { return GetVec<double>("mc_FSPartPx"); }

  virtual std::vector<double> GetFSPartPy() const { return GetVec<double>("mc_FSPartPy"); }

  virtual std::vector<double> GetFSPartPz() const { return GetVec<double>("mc_FSPartPz"); }

  virtual int GetNeutronReweightCategory(double neutKE) const {
    int ret = -999;

    /*
      std::map<int, std::string> Categs = {{11, "sigQE"},
      {18, "sig2p2h"},
      {19, "sigOther"},
      {1, "1chargePi"},
      {2, "1neutPi"},
      {3, "NPi"},
      {4, "SubThresh"},
      {5, "TrackableProt"}};
    */
    
    if (GetTruthNuPDG() != -14 || GetCurrent() != 1) return ret;
    
    int genie_n_muons = 0;
    int genie_n_piPM = 0;
    int genie_n_pi0 = 0;
    int genie_n_mesons = 0;
    int genie_n_heavy_baryons = 0;
    int genie_n_photons = 0;
    int genie_n_protons_above = 0;
    int genie_n_neutrons = 0;
    int genie_n_neutrons_above = 0;
    
    std::vector<int> PDGs = GetFSPartPDG();
    std::vector<double> Es = GetFSPartE();
    
    for (unsigned int i=0; i<PDGs.size(); ++i){
      int pdg = PDGs.at(i);
      double energy = Es.at(i);
      double proton_E = M_p +120.0;//Trackable Protons at 120 MeV T_p
      double neutron_E = M_n + neutKE;
      if ( abs(pdg) == 13) genie_n_muons++;
      else if ( pdg == 22  && energy > 10) genie_n_photons++;
      else if ( abs(pdg) == 211 ) genie_n_piPM++;
      else if ( pdg == 111 ) genie_n_pi0++;
      else if (abs(pdg) == 321 || abs(pdg) == 323 || pdg == 130 || pdg == 310 || pdg == 311 || pdg == 313 ){
	genie_n_mesons++;
      }
      else if ( pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 4112 || pdg == 4122 || pdg == 4222 || pdg == 411 || pdg == 421){
	genie_n_heavy_baryons++;
      }
      else if ( pdg == 2212 && energy > proton_E) genie_n_protons_above++;
      else if ( pdg == 2112 ) {
	if (energy > M_n) genie_n_neutrons++;
	if (energy > neutron_E) genie_n_neutrons_above++;
      }
    }

    if (genie_n_muons == 1 &&
        genie_n_piPM == 0 &&
        genie_n_pi0 == 0 &&
        genie_n_mesons == 0 &&
        genie_n_heavy_baryons == 0 &&
        genie_n_photons == 0 &&
        genie_n_protons_above == 0 &&
        genie_n_neutrons_above > 0) {
      ret = 10;
      ret += GetInteractionType();
      if (ret != 11 && ret != 18) ret = 19;
    }
    else if ( (genie_n_piPM + genie_n_pi0) > 1) ret = 3;
    else if ( genie_n_piPM == 1 ) ret = 1;
    else if ( genie_n_pi0 == 1 ) ret = 2;
    else if ( genie_n_heavy_baryons == 0 && genie_n_photons == 0 && genie_n_mesons == 0 && genie_n_muons == 1){
      if ( genie_n_protons_above > 0 ) ret = 5;
      else if ( genie_n_neutrons > 0 ) ret = 4;
    }
    
    return ret;
  }
  
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
  
  virtual double GetTotalNeutronKE() const {
    double tot_KE = 0.0;
    std::vector<int> PDGs = GetFSPartPDG();
    std::vector<double> Es = GetFSPartE();
    for (int iFS=0; iFS < PDGs.size(); ++iFS){
      double KE = Es.at(iFS)-M_n;
      if (PDGs.at(iFS) == 2112) tot_KE += KE;
    }
    return tot_KE*MeVGeV;
  }

  virtual double GetAvailableEnergy() const {
    double tot_E = 0.0;
    std::vector<int> PDGs = GetFSPartPDG();
    std::vector<double> Es = GetFSPartE();
    for (int iFS=0; iFS < PDGs.size(); ++iFS){
      if (PDGs.at(iFS)==2112 || PDGs.at(iFS) > 1000000000 || fabs(PDGs.at(iFS))==13) continue;
      else if (PDGs.at(iFS)==2212) tot_E += Es.at(iFS) - M_p;
      else if (fabs(PDGs.at(iFS))==211) tot_E += Es.at(iFS) - M_pi;
      else tot_E += Es.at(iFS);
    }
    return tot_E*MeVGeV;
  }

  virtual double GetRecoilProxy() const {
    return 0.1*GetTotalNeutronKE()+GetAvailableEnergy();
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

  virtual double GetVtxX() const { return GetVtx().at(0); }

  virtual double GetTrueVtxX() const { return GetTrueVtx().at(0); }

  virtual double GetVtxY() const { return GetVtx().at(1); }

  virtual double GetTrueVtxY() const { return GetTrueVtx().at(1); }

  virtual double GetVtxZ() const { return GetVtx().at(2); }

  virtual double GetTrueVtxZ() const { return GetTrueVtx().at(2); }

  virtual int GetIsMinosMatchTrack() const { return GetInt("isMinosMatchTrack"); }
  
  virtual int GetIsMinosMatchTrackOLD() const { return GetInt("muon_is_minos_match_track"); }
  
  virtual int GetIsMinosMatchStub() const { return GetInt("isMinosMatchStub"); }
  
  virtual int GetIsMinosMatchStubOLD() const { return GetInt("muon_is_minos_match_stub"); }

  virtual double GetCalRecoilEnergy() const{
    return GetDouble("recoil_energy_nonmuon_nonvtx100mm")+GetDouble("recoil_energy_nonmuon_nonvtx100mm_nuclTargs");
    //return GetDouble("recoil_energy_nonmuon_nonvtx100mm");//TEMP FOR DAN STUFF
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

  virtual double ApplyCaloTuning(double E) const{
    return E;
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

  virtual int GetDaisyPetal(double vtxX, double vtxY) const{
    int petal = -999;
    double angleInRad = TMath::Pi()+atan2(-vtxY,-vtxX); //signs are for 0 to 2Pi
    double angleInDeg = (angleInRad*180.0)/TMath::Pi();
    petal = (angleInDeg >= 0.0 && angleInDeg < 360.0) ? (int)(angleInDeg)/30 : petal; //12 daisy petals starting from 0 at x=1, y=0 on the unit circle. hence division by 30 degrees.
    return petal;
  }

  virtual double GetTrueDaisyPetal()const{
    double vtxX = GetTrueVtxX();
    double vtxY = GetTrueVtxY();
    return GetDaisyPetal(vtxX, vtxY);
  }

  virtual double GetRecoDaisyPetal() const{
    double vtxX = GetVtxX();
    double vtxY = GetVtxY();
    return GetDaisyPetal(vtxX, vtxY);
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

      m_LeadNeutIndex = leadNeutEIndex;//Setting so other functions can use. If they are called out of order, you will see no change.

      cands.push_back(GetNeutCand(leadNeutEIndex));
    }

    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  };

  //Need to set up something where this is what gets called by choice... Might just hard-code it for now... Need to run the effect of the GENIE drop too... Maybe I need to define the function with an input probability...
  virtual NeutronCandidates::NeutCands GetLeadNeutCandOnlyWithDrop(double prob=0.25, double thresh=10.0)
  {
    if (prob < 0.0) prob = 0.0;
    if (prob > 1.0) prob = 1.0;
    
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
          if (Es.at(idx) < thresh){
              int parentID = GetVecElemInt((toolName+branchNameParent).c_str(), idx);
              int ID = GetVecElemInt((toolName+branchNamePID).c_str(), idx);
              //std::cout << "parent: " << parentID << ", self: " << ID << std::endl;                                                                                                                              
              if ((parentID==2112 || ID==2112) && m_Random->Binomial(1,prob)){
		continue;
	      }
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

  virtual NeutronCandidates::NeutCands GetLeadNeutCandOnlyWithDropNoNearVertex(double prob=0.25, double thresh=10.0)
  {
    if (prob < 0.0) prob = 0.0;
    if (prob > 1.0) prob = 1.0;
    
    std::vector<NeutronCandidates::NeutCand> cands = {};
    int nBlobs = GetNNeutBlobs();

    std::vector<double> vtx = GetVtx();

    if (nBlobs > 0){
      std::vector<double> Es = GetNeutCandEs();
      std::string toolName = GetAnaToolName();
      std::string branchNameParent = "_BlobParentMCPID";
      std::string branchNamePID = "_BlobMCPID";
      int leadNeutEIndex = -999;
      double maxE = -999;
      for (unsigned int idx = 0; idx < Es.size(); ++idx){
        if (Es.at(idx) > maxE){
          if (Es.at(idx) < thresh){
              int parentID = GetVecElemInt((toolName+branchNameParent).c_str(), idx);
              int ID = GetVecElemInt((toolName+branchNamePID).c_str(), idx);
              //std::cout << "parent: " << parentID << ", self: " << ID << std::endl;                                                                                                                              
              if ((parentID==2112 || ID==2112) && m_Random->Binomial(1,prob)){
		continue;
	      }
          }

	  std::vector<double> zpos = GetNeutZPerCluster(idx);
	  std::vector<double> tpos = GetNeutTPosPerCluster(idx);
	  std::vector<int> view = GetNeutViewPerCluster(idx);

	  bool skipCand = false;
	  
	  for (int iClus=0; iClus < zpos.size(); ++iClus){
	    if (fabs(zpos.at(iClus) - vtx.at(2)) > 50.0) continue;//Only worry about clusters within 20mm (~1 plane) of vertex position in z for now.
	    double vtxTPosInView = GetTXFromXY(vtx.at(0), vtx.at(1), view.at(iClus));
	    if (fabs(tpos.at(iClus) - vtxTPosInView) < 50.0){//Skip candidates where there is a cluster within 50mm (~ 3 strips) of the vertex position in the view of the cluster.
	      skipCand = true;
	      //std::cout << "Skipping index: " << idx << "for being too close to the vertex" << std::endl;
	      break;
	    }
	  }

	  if (skipCand){
	    continue;
	  }
	  
          maxE = Es.at(idx);
          leadNeutEIndex = idx;
        }
      }

      if (leadNeutEIndex >= 0){
        cands.push_back(GetNeutCand(leadNeutEIndex));

	//std::cout << "Index: " << leadNeutEIndex << std::endl;
	//std::cout << "E: " << maxE << std::endl;
      }

      m_LeadNeutIndex = leadNeutEIndex;
    }

    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  }

  virtual NeutronCandidates::NeutCands GetLeadNeutCandOnlyWithDropNoNearVertexOrMuon(double prob=0.25, double thresh=10.0)
  {
    if (prob < 0.0) prob = 0.0;
    if (prob > 1.0) prob = 1.0;
    
    std::vector<NeutronCandidates::NeutCand> cands = {};
    int nBlobs = GetNNeutBlobs();

    std::vector<double> vtx = GetVtx();
    std::vector<double> muonMom = {GetMuon4V().X(), GetMuon4V().Y(), GetMuon4V().Z()};

    if (nBlobs > 0){
      std::vector<double> Es = GetNeutCandEs();
      std::string toolName = GetAnaToolName();
      std::string branchNameParent = "_BlobParentMCPID";
      std::string branchNamePID = "_BlobMCPID";
      int leadNeutEIndex = -999;
      double maxE = -999;
      for (unsigned int idx = 0; idx < Es.size(); ++idx){
        if (Es.at(idx) > maxE){
          if (Es.at(idx) < thresh){
              int parentID = GetVecElemInt((toolName+branchNameParent).c_str(), idx);
              int ID = GetVecElemInt((toolName+branchNamePID).c_str(), idx);
              //std::cout << "parent: " << parentID << ", self: " << ID << std::endl;                                                                                                                              
              if ((parentID==2112 || ID==2112) && m_Random->Binomial(1,prob)){
		continue;
	      }
          }

	  std::vector<double> zpos = GetNeutZPerCluster(idx);
	  std::vector<double> tpos = GetNeutTPosPerCluster(idx);
	  std::vector<int> view = GetNeutViewPerCluster(idx);

	  bool skipCand = false;
	  
	  for (int iClus=0; iClus < zpos.size(); ++iClus){
	    if (fabs(zpos.at(iClus) - vtx.at(2)) <= 50.0){
	      double vtxTPosInView = GetTXFromXY(vtx.at(0), vtx.at(1), view.at(iClus));
	      if (fabs(tpos.at(iClus) - vtxTPosInView) < 50.0){//Skip candidates where there is a cluster within 50mm (~ 3 strips) of the vertex position in the view of the cluster.
		skipCand = true;
		break;
	      }
	    }

	    //Only worry about muon contamination mostly downstream of the vertex...
	    if (zpos.at(iClus) - vtx.at(2) >= 0.0){
	      std::pair<double, double> muonXYAtZ = ProjectTrackToZ(vtx, muonMom, zpos.at(iClus));
	      double muonTPosAtZ = GetTXFromXY(muonXYAtZ.first, muonXYAtZ.second, view.at(iClus));
	      if (fabs(tpos.at(iClus) - muonTPosAtZ) < 150.0){
		skipCand = true;
		break;
	      }
	    }
	  }

	  if (skipCand) continue;
	  
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

  virtual NeutronCandidates::NeutCands GetLeadNeutCandOnlyWithDropNoNearVertexOrMuonOrECAL(double prob=0.25, double thresh=10.0)
  {
    if (prob < 0.0) prob = 0.0;
    if (prob > 1.0) prob = 1.0;
    
    std::vector<NeutronCandidates::NeutCand> cands = {};
    int nBlobs = GetNNeutBlobs();

    std::vector<double> vtx = GetVtx();
    std::vector<double> muonMom = {GetMuon4V().X(), GetMuon4V().Y(), GetMuon4V().Z()};

    if (nBlobs > 0){
      std::vector<double> Es = GetNeutCandEs();
      std::string toolName = GetAnaToolName();
      std::string branchNameParent = "_BlobParentMCPID";
      std::string branchNamePID = "_BlobMCPID";
      int leadNeutEIndex = -999;
      double maxE = -999;
      for (unsigned int idx = 0; idx < Es.size(); ++idx){
	//Added in this requirement of the energy being above 2 MeV since I see 0 3D neutrons above that point.
        if (Es.at(idx) > maxE && Es.at(idx) >= 2.0){
          if (Es.at(idx) < thresh){
              int parentID = GetVecElemInt((toolName+branchNameParent).c_str(), idx);
              int ID = GetVecElemInt((toolName+branchNamePID).c_str(), idx);
              //std::cout << "parent: " << parentID << ", self: " << ID << std::endl;                                                                                                                              
              if ((parentID==2112 || ID==2112) && m_Random->Binomial(1,prob)){
		continue;
	      }
          }

	  std::vector<double> zpos = GetNeutZPerCluster(idx);
	  std::vector<double> tpos = GetNeutTPosPerCluster(idx);
	  std::vector<int> view = GetNeutViewPerCluster(idx);

	  bool skipCand = false;
	  
	  for (int iClus=0; iClus < zpos.size(); ++iClus){
	    //Just skip anything near/around/in the ECAL.
	    if (zpos.at(iClus) > 8422){
	      skipCand=true;
	      //std::cout << "Skipping index: " << idx << " for being in the ECAL" << std::endl;
	      break;
	    }
	    
	    if (fabs(zpos.at(iClus) - vtx.at(2)) <= 50.0){
	      double vtxTPosInView = GetTXFromXY(vtx.at(0), vtx.at(1), view.at(iClus));
	      if (fabs(tpos.at(iClus) - vtxTPosInView) < 50.0){//Skip candidates where there is a cluster within 50mm (~ 3 strips) of the vertex position in the view of the cluster.
		skipCand = true;
		//std::cout << "Skipping index: " << idx << " for being too close to the vertex" << std::endl;
		break;
	      }
	    }
	    
	    //Only worry about muon contamination mostly downstream of the vertex...
	    if (zpos.at(iClus) - vtx.at(2) >= 0.0){
	      std::pair<double, double> muonXYAtZ = ProjectTrackToZ(vtx, muonMom, zpos.at(iClus));
	      double muonTPosAtZ = GetTXFromXY(muonXYAtZ.first, muonXYAtZ.second, view.at(iClus));
	      double cutAtThisZ = GetTCutFromAngle(15, muonMom, zpos.at(iClus)-vtx.at(2));
	      if (fabs(tpos.at(iClus) - muonTPosAtZ) < cutAtThisZ){
	      //if (fabs(tpos.at(iClus) - muonTPosAtZ) < 100.0){
		skipCand = true;
		//std::cout << "Skipping index: " << idx << " for being too close to the muon" << std::endl;
		break;
	      }
	    }
	  }

	  if (skipCand) continue;

          maxE = Es.at(idx);
          leadNeutEIndex = idx;
        }
      }

      if (leadNeutEIndex >= 0){
        cands.push_back(GetNeutCand(leadNeutEIndex));

	//std::cout << "Index: " << leadNeutEIndex << std::endl;
	//std::cout << "E: " << maxE << std::endl;
      }

      m_LeadNeutIndex = leadNeutEIndex;
    }

    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  }

  //This is so that when using the CV with the drop as above, the same neutrons are dropped in all universes... see runEventLoop for implementation
  virtual NeutronCandidates::NeutCands GetLeadNeutCandOnlyFromIndex(int idx)
  {
    std::vector<NeutronCandidates::NeutCand> cands = {};
    int nBlobs = GetNNeutBlobs();

    if (nBlobs < idx){
      NeutronCandidates::NeutCands dummy(cands);
      return cands;
    }
    else{
      cands.push_back(GetNeutCand(idx));
    }

    m_LeadNeutIndex = idx;

    NeutronCandidates::NeutCands EvtCands(cands);
    return EvtCands;
  }
  
  //Returns the total energy of the candidate to fill a variable that is reco when just trying to make efficiency as function of true neutron energy :)
  virtual double GetLeadNeutCandE() const{

    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchEName = "_BlobTotalE";
      std::string branchPIDName = "_BlobMCPID";

      double factor = 1.0;
      // Change below to have a shift in the reconstructed energy for protons that make neutron candidates.
      //if (GetVecElem((toolName+branchPIDName).c_str(), m_LeadNeutIndex) == 2212) factor = 1.0-0.035;//Temp change to test shifting proton energy by 3.5% down to see if this covers the discrepancy in the shape of this variable.
      
      return factor*GetVecElem((toolName+branchEName).c_str(), m_LeadNeutIndex);
    }

    return -999;
  };
  
  //Returns the total energy of the candidate to fill a variable that is reco when just trying to make efficiency as function of true neutron energy :)
  virtual TVector3 GetLeadNeutCandPos() const{
    
    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchNameX = "_BlobBegX";
      std::string branchNameY = "_BlobBegY";
      std::string branchNameZ = "_BlobBegZ";
      double x = GetVecElem((toolName+branchNameX).c_str(), m_LeadNeutIndex);//Consistent calculation for KE from other comparisons.
      double y = GetVecElem((toolName+branchNameY).c_str(), m_LeadNeutIndex);//Consistent calculation for KE from other comparisons.
      double z = GetVecElem((toolName+branchNameZ).c_str(), m_LeadNeutIndex);//Consistent calculation for KE from other comparisons.
      
      TVector3 pos(x,y,z);
      return pos;
    }

    TVector3 dummy(-999,-999,-999);
    return dummy;
  };
  
  virtual TVector3 GetLeadNeutCandFlightPath() const{
    TVector3 pos = GetLeadNeutCandPos();
    std::vector<double> vtx = GetVtx();
    TVector3 vtx3(vtx.at(0),vtx.at(1),vtx.at(2));
    TVector3 FP = pos-vtx3;
    return FP;
  };

  virtual double GetLeadNeutCandAngleToMuon() const{
    TVector3 muonMom(GetMuon4V().X(),GetMuon4V().Y(),GetMuon4V().Z());
    TVector3 FP=GetLeadNeutCandFlightPath();
    if (FP.Mag() == 0 || muonMom.Mag() == 0) return -999;
    else return FP.Angle(muonMom);
  }; 

  virtual double GetLeadNeutCandZPos() const{
    TVector3 pos = GetLeadNeutCandPos();
    return pos.Z();
  }
  
  virtual double GetLeadNeutVtxZDist() const{
    TVector3 FP = GetLeadNeutCandFlightPath();
    return fabs(FP.Z());
  };

  virtual double GetLeadNeutVtxDist() const{
    TVector3 FP = GetLeadNeutCandFlightPath();
    return FP.Mag();
  };

  virtual std::vector<double> GetLeadNeutTPosPerCluster() const{
    std::vector<double> ret;
    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchTPos = "_BlobTPosPerCluster";
      std::string branchNClus = "_BlobNClusters";
      int idx1=0;
      
      std::vector<int> nClus = GetVec<int>((toolName+branchNClus).c_str());
      for (int iCand=0; iCand < m_LeadNeutIndex; ++iCand){
	idx1 += nClus.at(iCand);
      }
      
      std::vector<double> tPos = GetVec<double>((toolName+branchTPos).c_str());
      for (int iClus=idx1; iClus < idx1+nClus.at(m_LeadNeutIndex); ++iClus){
	ret.push_back(tPos.at(iClus));
      }
    }
    return ret;
  }

  virtual std::vector<int> GetLeadNeutViewPerCluster() const{
    std::vector<int> ret;
    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchView = "_BlobViewPerCluster";
      std::string branchNClus = "_BlobNClusters";
      int idx1=0;
      
      std::vector<int> nClus = GetVec<int>((toolName+branchNClus).c_str());
      for (int iCand=0; iCand < m_LeadNeutIndex; ++iCand){
	idx1 += nClus.at(iCand);
      }
      
      std::vector<int> view = GetVec<int>((toolName+branchView).c_str());
      for (int iClus=idx1; iClus < idx1+nClus.at(m_LeadNeutIndex); ++iClus){
	ret.push_back(view.at(iClus));
      }
    }
    return ret;
  }

  virtual std::vector<double> GetLeadNeutTimePerCluster() const{
    std::vector<double> ret;
    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchTime = "_BlobTimePerCluster";
      std::string branchNClus = "_BlobNClusters";
      int idx1=0;
      
      std::vector<int> nClus = GetVec<int>((toolName+branchNClus).c_str());
      for (int iCand=0; iCand < m_LeadNeutIndex; ++iCand){
	idx1 += nClus.at(iCand);
      }
      
      std::vector<double> time = GetVec<double>((toolName+branchTime).c_str());
      for (int iClus=idx1; iClus < idx1+nClus.at(m_LeadNeutIndex); ++iClus){
	ret.push_back(time.at(iClus));
      }
    }
    return ret;
  }
  
  virtual std::vector<double> GetLeadNeutEnergyPerCluster() const{
    std::vector<double> ret;
    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchEnergy = "_BlobEnergyPerCluster";
      std::string branchNClus = "_BlobNClusters";
      int idx1=0;
      
      std::vector<int> nClus = GetVec<int>((toolName+branchNClus).c_str());
      for (int iCand=0; iCand < m_LeadNeutIndex; ++iCand){
	idx1 += nClus.at(iCand);
      }
      
      std::vector<double> time = GetVec<double>((toolName+branchEnergy).c_str());
      for (int iClus=idx1; iClus < idx1+nClus.at(m_LeadNeutIndex); ++iClus){
	ret.push_back(time.at(iClus));
      }
    }
    return ret;
  }

  virtual std::vector<double> GetLeadNeutZPerCluster() const{
    std::vector<double> ret;
    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchZPos = "_BlobZPerCluster";
      std::string branchNClus = "_BlobNClusters";
      int idx1=0;
      
      std::vector<int> nClus = GetVec<int>((toolName+branchNClus).c_str());
      for (int iCand=0; iCand < m_LeadNeutIndex; ++iCand){
	idx1 += nClus.at(iCand);
      }
      
      std::vector<double> zPos = GetVec<double>((toolName+branchZPos).c_str());
      for (int iClus=idx1; iClus < idx1+nClus.at(m_LeadNeutIndex); ++iClus){
	ret.push_back(zPos.at(iClus));
      }
    }
    return ret;
  }

  virtual std::vector<double> GetNeutZPerCluster(int index) const{
    std::vector<double> ret;
    if (index >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchZPos = "_BlobZPerCluster";
      std::string branchNClus = "_BlobNClusters";
      int idx1=0;
      
      std::vector<int> nClus = GetVec<int>((toolName+branchNClus).c_str());
      for (int iCand=0; iCand < index; ++iCand){
	idx1 += nClus.at(iCand);
      }
      
      std::vector<double> zPos = GetVec<double>((toolName+branchZPos).c_str());
      for (int iClus=idx1; iClus < idx1+nClus.at(index); ++iClus){
	ret.push_back(zPos.at(iClus));
      }
    }
    return ret;
  }

  virtual std::vector<double> GetNeutTPosPerCluster(int index) const{
    std::vector<double> ret;
    if (index >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchTPos = "_BlobTPosPerCluster";
      std::string branchNClus = "_BlobNClusters";
      int idx1=0;
      
      std::vector<int> nClus = GetVec<int>((toolName+branchNClus).c_str());
      for (int iCand=0; iCand < index; ++iCand){
	idx1 += nClus.at(iCand);
      }
      
      std::vector<double> tPos = GetVec<double>((toolName+branchTPos).c_str());
      for (int iClus=idx1; iClus < idx1+nClus.at(index); ++iClus){
	ret.push_back(tPos.at(iClus));
      }
    }
    return ret;
  }

  virtual std::vector<int> GetNeutViewPerCluster(int index) const{
    std::vector<int> ret;
    if (index >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchView = "_BlobViewPerCluster";
      std::string branchNClus = "_BlobNClusters";
      int idx1=0;
      
      std::vector<int> nClus = GetVec<int>((toolName+branchNClus).c_str());
      for (int iCand=0; iCand < index; ++iCand){
	idx1 += nClus.at(iCand);
      }
      
      std::vector<int> view = GetVec<int>((toolName+branchView).c_str());
      for (int iClus=idx1; iClus < idx1+nClus.at(index); ++iClus){
	ret.push_back(view.at(iClus));
      }
    }
    return ret;
  }
  
  virtual double GetMATCHEDLeadNeutCandE() const{

    if (m_LeadNeutIndex >= 0){
      std::string toolName = GetAnaToolName();
      std::string branchName = "_BlobMCTopTrackE";
      
      return (GetVecElem((toolName+branchName).c_str(), m_LeadNeutIndex)-M_n);//Consistent calculation for KE from other comparisons.
    }
    return -999;
  }
  
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

  double GetNeutronNormWeight() const{
    double ret = 1.0;
    return ret;//Temporary to Remake the Renormalization Plots In The Face of Extending The Reweight Beyond The 200 MeV cutoff...
    int categ = GetNeutronReweightCategory(10.0);
    TH2D* reweightHist = m_NeutRWHists.at(categ);
    if (!reweightHist) return ret;
    int binX = reweightHist->GetXaxis()->FindBin(GetMuonPTTrue());
    int binY = reweightHist->GetYaxis()->FindBin(GetMaxFSNeutronKE());
    double val = reweightHist->GetBinContent(binX,binY);
    ret = std::max(0.0, val);
    if (ret==0.0) ret = 1.0;
    return ret;
  }
  
  private:
  TRandom3* m_Random;
  
  //Still needed for some systematics to compile, but shouldn't be used for reweighting anymore.
  protected:
  #include "PlotUtils/WeightFunctions.h" // Get*Weight
};

#endif
