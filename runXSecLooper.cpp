//TODO: ADD THE NORMALIZATION VALUE FROM THE TARGETS HERE. SELF-NORMALIZE.

#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include <PlotUtils/TargetUtils.h>
#include <PlotUtils/PlotUtilsPhysicalConstants.h>

#include "util/GetRecoTargetZ.h"

#include <cstdlib>
typedef unsigned int uint;

//Borrowed largely from NSFNukeCCInclusive Code in Anezka's working directory
double GetNormFactor(int tgtZ){
  //Borrowed instead from Aaron Bercellie. It would make sense that the number of tracker targets utilized in the GENIE event rates is the entire simulated tracker, that explains why none of the options I tried resulted in a factor that made sense. Aaron's code seemed to also 
  const double tracker_mass = PlotUtils::TargetUtils::Get().GetTrackerMass(PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true, 850.);
  const double nAtoms_C = tracker_mass * PlotUtils::TargetUtils::Get().GetTrackerElementMassFraction( 6, true ) * PlotUtils::TargetProp::AtomsPerGram::C;
  const double nAtoms_C_temp = PlotUtils::TargetUtils::Get().GetTrackerElementNAtoms( 6, PlotUtils::TargetProp::Tracker::Face, PlotUtils::TargetProp::Tracker::Back, true, 850.);
  double trackerAtomsC = PlotUtils::TargetUtils::Get().GetTrackerElementNAtoms( 6, 5980, 8422, true, 850.0);

  std::cout << "Tracker C12 Atoms 3 ways: " << nAtoms_C << ", " << nAtoms_C_temp << ", " << trackerAtomsC << std::endl;

  double nucleons = 0.0;
  if (tgtZ == 6){
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(3, tgtZ, true, 850.0);
  }
  else if (tgtZ == 26){
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(1, tgtZ, true, 850.0);
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(2, tgtZ, true, 850.0);
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(3, tgtZ, true, 850.0);
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(5, tgtZ, true, 850.0);
  }
  else if (tgtZ == 82){
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(1, tgtZ, true, 850.0);
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(2, tgtZ, true, 850.0);
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(3, tgtZ, true, 850.0);
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(4, tgtZ, true, 850.0);
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(5, tgtZ, true, 850.0);
  }
  else if (tgtZ == 8){
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(6, tgtZ, true, 850.0);
    //nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(6, 1, true, 850.0); This was wrong. The oxygen and hydrgoen are not handled separately. They are handled all at once.
  }
  else if (tgtZ == -1){
    nucleons += PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 5980, 8422. , true, 850.0 ); //Hard-coding to see effect of requiring this compared to the normalization factor calculated by xSecLooper with it's own 
  }

  //double output = (nucleons > 0) ? trackerAtomsC/nucleons : 1.0;
  double output = (nucleons > 0) ? nAtoms_C/nucleons : 1.0;
  return output;
}


std::map<TString,int> IntTypeMap = {{"QE",1},{"RES",2},{"DIS",3},{"2p2h",8},{"Other",-999}};

class MinModDepCCQEXSec : public XSec
{
public:
  MinModDepCCQEXSec(const char* name, const int tgtZ, const double neutKE = 10.0, const TString intType = "")
    :XSec(name), fTgtZ(tgtZ), fIntType(intType), fNeutKE(neutKE), fApothem(850.0)
  {
  };

  bool inApothem( ChainWrapper& chw, int entry ){
    double x=fabs(chw.GetValue("mc_vtx", entry, 0));
    double y=fabs(chw.GetValue("mc_vtx", entry, 1));

    if(x*x + y*y < fApothem*fApothem) return true;

    double lenOfSide = fApothem * ( 2 / sqrt(3) );

    if( x > fApothem )
      return false;

    if( y < lenOfSide/2.0 )
      return true;

    double slope = (lenOfSide / 2.0) / fApothem;
    if( y < lenOfSide - x*slope )
      return true;

    return false;    
  }

  //Hard-coding correct target. This uses the same function I'm using in the selection of my truth in the event loop is there something conceptually wrong aobut that...?
  bool isCorrectTgt( ChainWrapper& chw, int entry ){
    double x=chw.GetValue("mc_vtx", entry, 0);
    double y=chw.GetValue("mc_vtx", entry, 1);
    double z=chw.GetValue("mc_vtx", entry, 2);
    int trueTgtZ = chw.GetValue("mc_targetZ", entry);

    if (fTgtZ == -1){
      return (z > 5980 && z < 8422);
    }
    else if (fTgtZ == 8){
      if ((trueTgtZ != fTgtZ) && (trueTgtZ != 1)) return false;
      int trueTgtCode = util::GetTrueTgtCode(trueTgtZ, x, y, z);//Is this problematic to use the exact same definition? Feels self-fulfilling...
      return (trueTgtCode == 6666);
    }
    else {
      if (fTgtZ != trueTgtZ) return false;
      int trueTgtCode = util::GetTrueTgtCode(trueTgtZ, x, y, z);//Is this problematic to use the exact same definition? Feels self-fulfilling...
      return (trueTgtCode > 0 && trueTgtCode != 6666); //Water handled above.
    }
    return false;
  }

  bool isCCInclusiveSignal( ChainWrapper& chw, int entry )
  {
    double theta              = 0.;
    double true_muon_px   = (double)chw.GetValue("mc_primFSLepton",entry,0)/1000;
    double true_muon_py   = (double)chw.GetValue("mc_primFSLepton",entry,1)/1000;
    double true_muon_pz   = (double)chw.GetValue("mc_primFSLepton",entry,2)/1000;
    double numi_beam_angle_rad = -0.05887;
    double pyprime = -1.0*sin(numi_beam_angle_rad)*true_muon_pz + cos(numi_beam_angle_rad)*true_muon_py;
    double pzprime =  1.0*cos(numi_beam_angle_rad)*true_muon_pz + sin(numi_beam_angle_rad)*true_muon_py;
    double pSquare = pow(true_muon_px,2) + pow(pyprime,2) + pow(pzprime,2);
    double pMag = sqrt(pSquare);
    theta = acos( pzprime / pMag );
    theta *= 180./3.14159;
    
    //if(!chw.GetValue("truth_is_fiducial",entry)) return false; //Doesn't work for MasterAnaDev tuples.  What does this even mean in the targets anyway? :(
    //if( pzprime >= 1.5 && theta <= 20.0 ) return true;
    if( pMag > 1.5 && pMag < 20.0 && theta <= 17.0 ) return true;
    return false;
    
  }
  
  bool isFSSignal( ChainWrapper& chw, int entry )
  {
    int genie_n_muons = 0;
    int genie_n_mesons = 0;
    int genie_n_heavy_baryons_plus_pi0s = 0;
    int genie_n_photons = 0;
    int genie_n_protons = 0;
    int genie_n_neutrons = 0;

    int nFSPart = (int)chw.GetValue("mc_nFSPart", entry);

    int intType = (int)chw.GetValue("mc_intType",entry);

    const int corrIntType = IntTypeMap[fIntType];
    if (corrIntType == -999){
      if (intType == 1 || intType == 2 || intType == 3 || intType == 8) return false;
    }
    else if (corrIntType != 0){
      if (intType != corrIntType) return false;
    }

    for (unsigned int i=0; i<nFSPart; ++i){
      int pdg = (int)chw.GetValue("mc_FSPartPDG",entry,i);
      double energy = (double)chw.GetValue("mc_FSPartE",entry,i);
      double proton_E = 1058.272;
      double neutron_E = 939.57+fNeutKE;//Testing not hard coded neut KE
      if (abs(pdg) == 13) genie_n_muons++;
      else if ( pdg == 22  && energy > 10) genie_n_photons++;
      else if ( abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 323 || pdg == 111 || pdg == 130 || pdg == 310 || pdg == 311 || pdg == 313 ){
	genie_n_mesons++;
      }
      else if ( pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 4112 || pdg == 4122 || pdg == 4222 || pdg == 411 || pdg == 421 || pdg == 111){
	genie_n_heavy_baryons_plus_pi0s++;
      }
      else if ( pdg == 2212 && energy > proton_E) genie_n_protons++;
      else if ( pdg == 2112 && energy > neutron_E) genie_n_neutrons++;
    }

      return genie_n_muons == 1 &&
	genie_n_mesons == 0 &&
	genie_n_heavy_baryons_plus_pi0s == 0 &&
	genie_n_photons == 0 &&
	genie_n_protons == 0 &&
	genie_n_neutrons > 0; //Hard-coded to include > 0 neutrons for now.
  }

  // Override this method from the base class to decide what events to
  // include in this selection
  virtual bool passesCuts(ChainWrapper& chw, int entry)
  {
    if((int)chw.GetValue("mc_incoming", entry)!=-14) return false;
    if((int)chw.GetValue("mc_current", entry)!=1) return false;
    if(!isCCInclusiveSignal  ( chw, entry ) ) return false;
    if(!inApothem(chw,entry)) return false;
    if(!isCorrectTgt(chw,entry)) return false;
    if(!isFSSignal(chw, entry)) return false;
    
    return true;
  }

  const int fTgtZ;

  const TString fIntType;

  const double fNeutKE;

  const double fApothem;

};

int main(const int argc, const char** argv)
{

  /*
  std::cout << GetNormFactor(-1) << std::endl;
  return 1;
  */

  //Read a playlist file from the command line
  if(argc < 5)
  {
    std::cerr << "Expected at least 4 command line argument, but got " << argc - 1 << ".\n\n"
              << "USAGE: runXSecLooper <outName> <playlist> <neutKE> <MCPlaylists.txt>\n\n"
              << "MCPlaylists.txt are a list of files which shall contain one .root file per line that has a Truth tree in it.\n"
              << "This program returns 0 when it suceeds.  It produces a .root file with GENIEXSECEXTRACT in its name.\n";
    return 1;
  }

  string outName = argv[1];
  TString pList = argv[2];
  double neutKE = atof(argv[3]);

  std::vector<const char*> fileNames(argv+4,argv+argc);

  // Create the XSecLooper and tell it the input files
  // Inputs should be the merged ntuples:
  XSecLooper loop(fileNames[0]);
  for(auto whichFile = fileNames.begin()+1; whichFile != fileNames.end(); ++whichFile)
    {
      loop.addFiles(*whichFile);
    }


  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(-14);
  if (pList.Contains("5A")){
    loop.setPlaylist(PlotUtils::FluxReweighter::minervame5A);
  }
  else if (pList.Contains("6")){
    loop.setPlaylist(PlotUtils::FluxReweighter::minervame6A);
  }
  else if (pList.Contains("All")){
    loop.setPlaylist(PlotUtils::FluxReweighter::minervame6A);
  }
  else{
    cout << " Unknown playlist value." << endl;
    return 86;
  }


  // Setting the number of Universes in the GENIE error band (default 100, put 0 if you do not want to include the universes)
  loop.setNumUniv(0); 
  loop.setFiducial(5980, 8422);//

  // Add the differential cross section dsigma/ds_dpT
  //double pt_edges[] = { 0.0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 1.25, 1.5, 2.5, 4.5 };
  //int pt_nbins = 14; 

  /* REMOVING UNUSED BINNINGS SINCE THE LATEST HAS BEEN FINALIZED
  //Hard-coded to latest optimization for carbon migration matrix...
  double pt_edges_carbon[] = { 0.0, 0.09, 0.18, 0.25, 0.34, 0.425, 0.515, 0.64, 0.78, 0.94, 1.15, 1.5};
  int pt_nbins_carbon = 11;

  //Hard-coded to latest optimization for carbon migration matrix...
  double pt_edges_combined[] = { 0.0, 0.18, 0.34, 0.515, 0.64, 0.78, 0.94, 1.15, 1.5};
  int pt_nbins_combined = 8;

  //Coded for testing against the latest before the Carbon bin optimization effort
  double pt_edges_std[] = { 0.0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 1.25, 1.5, 2.5};
  int pt_nbins_std = 13;

  */
  //Hard-coded to latest optimization for carbon migration matrix...
  double pt_edges[] = {0, 0.125, 0.25, 0.38, 0.515, 0.64, 0.78, 0.94, 1.15, 1.5};
  int pt_nbins = 9;

  MinModDepCCQEXSec* ds_dpT = new MinModDepCCQEXSec("pTmu_Tracker", -1, neutKE);
  ds_dpT->setBinEdges(pt_nbins, pt_edges);
  ds_dpT->setVariable(XSec::kPTLep);
  ds_dpT->setIsFluxIntegrated(true);
  ds_dpT->setDimension(1);
  ds_dpT->setFluxIntLimits(0.0, 100.0);
  ds_dpT->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT->setNormalizationValue(GetNormFactor(-1));
  ds_dpT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT);

  MinModDepCCQEXSec* ds_dpT_QE = new MinModDepCCQEXSec("pTmu_Tracker_QE", -1, neutKE, "QE");
  ds_dpT_QE->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_QE->setVariable(XSec::kPTLep);
  ds_dpT_QE->setIsFluxIntegrated(true);
  ds_dpT_QE->setDimension(1);
  ds_dpT_QE->setFluxIntLimits(0.0, 100.0);
  ds_dpT_QE->setNormalizationType(XSec::kSelfNorm);
  ds_dpT_QE->setNormalizationValue(GetNormFactor(-1));
  ds_dpT_QE->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_QE);

  MinModDepCCQEXSec* ds_dpT_RES = new MinModDepCCQEXSec("pTmu_Tracker_RES", -1, neutKE, "RES");
  ds_dpT_RES->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_RES->setVariable(XSec::kPTLep);
  ds_dpT_RES->setIsFluxIntegrated(true);
  ds_dpT_RES->setDimension(1);
  ds_dpT_RES->setFluxIntLimits(0.0, 100.0);
  ds_dpT_RES->setNormalizationType(XSec::kSelfNorm);
  ds_dpT_RES->setNormalizationValue(GetNormFactor(-1));
  ds_dpT_RES->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_RES);

  MinModDepCCQEXSec* ds_dpT_DIS = new MinModDepCCQEXSec("pTmu_Tracker_DIS", -1, neutKE, "DIS");
  ds_dpT_DIS->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_DIS->setVariable(XSec::kPTLep);
  ds_dpT_DIS->setIsFluxIntegrated(true);
  ds_dpT_DIS->setDimension(1);
  ds_dpT_DIS->setFluxIntLimits(0.0, 100.0);
  ds_dpT_DIS->setNormalizationType(XSec::kSelfNorm);
  ds_dpT_DIS->setNormalizationValue(GetNormFactor(-1));
  ds_dpT_DIS->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_DIS);

  MinModDepCCQEXSec* ds_dpT_2p2h = new MinModDepCCQEXSec("pTmu_Tracker_2p2h", -1, neutKE, "2p2h");
  ds_dpT_2p2h->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_2p2h->setVariable(XSec::kPTLep);
  ds_dpT_2p2h->setIsFluxIntegrated(true);
  ds_dpT_2p2h->setDimension(1);
  ds_dpT_2p2h->setFluxIntLimits(0.0, 100.0);
  ds_dpT_2p2h->setNormalizationType(XSec::kSelfNorm);
  ds_dpT_2p2h->setNormalizationValue(GetNormFactor(-1));
  ds_dpT_2p2h->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_2p2h);

  MinModDepCCQEXSec* ds_dpT_Other = new MinModDepCCQEXSec("pTmu_Tracker_Other", -1, neutKE, "Other");
  ds_dpT_Other->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Other->setVariable(XSec::kPTLep);
  ds_dpT_Other->setIsFluxIntegrated(true);
  ds_dpT_Other->setDimension(1);
  ds_dpT_Other->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Other->setNormalizationType(XSec::kSelfNorm);
  ds_dpT_Other->setNormalizationValue(GetNormFactor(-1));
  ds_dpT_Other->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Other);

  MinModDepCCQEXSec* ds_dpT_C = new MinModDepCCQEXSec("pTmu_C", 6, neutKE);
  ds_dpT_C->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_C->setVariable(XSec::kPTLep);
  ds_dpT_C->setIsFluxIntegrated(true);
  ds_dpT_C->setDimension(1);
  ds_dpT_C->setFluxIntLimits(0.0, 100.0);
  ds_dpT_C->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_C->setNormalizationValue(GetNormFactor(6));
  ds_dpT_C->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_C);

  MinModDepCCQEXSec* ds_dpT_C_QE = new MinModDepCCQEXSec("pTmu_C_QE", 6, neutKE,"QE");
  ds_dpT_C_QE->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_C_QE->setVariable(XSec::kPTLep);
  ds_dpT_C_QE->setIsFluxIntegrated(true);
  ds_dpT_C_QE->setDimension(1);
  ds_dpT_C_QE->setFluxIntLimits(0.0, 100.0);
  ds_dpT_C_QE->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_C_QE->setNormalizationValue(GetNormFactor(6));
  ds_dpT_C_QE->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_C_QE);

  MinModDepCCQEXSec* ds_dpT_C_RES = new MinModDepCCQEXSec("pTmu_C_RES", 6, neutKE,"RES");
  ds_dpT_C_RES->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_C_RES->setVariable(XSec::kPTLep);
  ds_dpT_C_RES->setIsFluxIntegrated(true);
  ds_dpT_C_RES->setDimension(1);
  ds_dpT_C_RES->setFluxIntLimits(0.0, 100.0);
  ds_dpT_C_RES->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_C_RES->setNormalizationValue(GetNormFactor(6));
  ds_dpT_C_RES->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_C_RES);

  MinModDepCCQEXSec* ds_dpT_C_DIS = new MinModDepCCQEXSec("pTmu_C_DIS", 6, neutKE,"DIS");
  ds_dpT_C_DIS->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_C_DIS->setVariable(XSec::kPTLep);
  ds_dpT_C_DIS->setIsFluxIntegrated(true);
  ds_dpT_C_DIS->setDimension(1);
  ds_dpT_C_DIS->setFluxIntLimits(0.0, 100.0);
  ds_dpT_C_DIS->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_C_DIS->setNormalizationValue(GetNormFactor(6));
  ds_dpT_C_DIS->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_C_DIS);

  MinModDepCCQEXSec* ds_dpT_C_2p2h = new MinModDepCCQEXSec("pTmu_C_2p2h", 6, neutKE,"2p2h");
  ds_dpT_C_2p2h->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_C_2p2h->setVariable(XSec::kPTLep);
  ds_dpT_C_2p2h->setIsFluxIntegrated(true);
  ds_dpT_C_2p2h->setDimension(1);
  ds_dpT_C_2p2h->setFluxIntLimits(0.0, 100.0);
  ds_dpT_C_2p2h->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_C_2p2h->setNormalizationValue(GetNormFactor(6));
  ds_dpT_C_2p2h->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_C_2p2h);

  MinModDepCCQEXSec* ds_dpT_C_Other = new MinModDepCCQEXSec("pTmu_C_Other", 6, neutKE,"Other");
  ds_dpT_C_Other->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_C_Other->setVariable(XSec::kPTLep);
  ds_dpT_C_Other->setIsFluxIntegrated(true);
  ds_dpT_C_Other->setDimension(1);
  ds_dpT_C_Other->setFluxIntLimits(0.0, 100.0);
  ds_dpT_C_Other->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_C_Other->setNormalizationValue(GetNormFactor(6));
  ds_dpT_C_Other->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_C_Other);

  MinModDepCCQEXSec* ds_dpT_Fe = new MinModDepCCQEXSec("pTmu_Iron", 26, neutKE);
  ds_dpT_Fe->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Fe->setVariable(XSec::kPTLep);
  ds_dpT_Fe->setIsFluxIntegrated(true);
  ds_dpT_Fe->setDimension(1);
  ds_dpT_Fe->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Fe->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Fe->setNormalizationValue(GetNormFactor(26));
  ds_dpT_Fe->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Fe);

  MinModDepCCQEXSec* ds_dpT_Fe_QE = new MinModDepCCQEXSec("pTmu_Iron_QE", 26, neutKE, "QE");
  ds_dpT_Fe_QE->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Fe_QE->setVariable(XSec::kPTLep);
  ds_dpT_Fe_QE->setIsFluxIntegrated(true);
  ds_dpT_Fe_QE->setDimension(1);
  ds_dpT_Fe_QE->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Fe_QE->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Fe_QE->setNormalizationValue(GetNormFactor(26));
  ds_dpT_Fe_QE->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Fe_QE);

  MinModDepCCQEXSec* ds_dpT_Fe_RES = new MinModDepCCQEXSec("pTmu_Iron_RES", 26, neutKE, "RES");
  ds_dpT_Fe_RES->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Fe_RES->setVariable(XSec::kPTLep);
  ds_dpT_Fe_RES->setIsFluxIntegrated(true);
  ds_dpT_Fe_RES->setDimension(1);
  ds_dpT_Fe_RES->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Fe_RES->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Fe_RES->setNormalizationValue(GetNormFactor(26));
  ds_dpT_Fe_RES->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Fe_RES);

  MinModDepCCQEXSec* ds_dpT_Fe_DIS = new MinModDepCCQEXSec("pTmu_Iron_DIS", 26, neutKE, "DIS");
  ds_dpT_Fe_DIS->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Fe_DIS->setVariable(XSec::kPTLep);
  ds_dpT_Fe_DIS->setIsFluxIntegrated(true);
  ds_dpT_Fe_DIS->setDimension(1);
  ds_dpT_Fe_DIS->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Fe_DIS->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Fe_DIS->setNormalizationValue(GetNormFactor(26));
  ds_dpT_Fe_DIS->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Fe_DIS);

  MinModDepCCQEXSec* ds_dpT_Fe_2p2h = new MinModDepCCQEXSec("pTmu_Iron_2p2h", 26, neutKE, "2p2h");
  ds_dpT_Fe_2p2h->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Fe_2p2h->setVariable(XSec::kPTLep);
  ds_dpT_Fe_2p2h->setIsFluxIntegrated(true);
  ds_dpT_Fe_2p2h->setDimension(1);
  ds_dpT_Fe_2p2h->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Fe_2p2h->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Fe_2p2h->setNormalizationValue(GetNormFactor(26));
  ds_dpT_Fe_2p2h->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Fe_2p2h);

  MinModDepCCQEXSec* ds_dpT_Fe_Other = new MinModDepCCQEXSec("pTmu_Iron_Other", 26, neutKE, "Other");
  ds_dpT_Fe_Other->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Fe_Other->setVariable(XSec::kPTLep);
  ds_dpT_Fe_Other->setIsFluxIntegrated(true);
  ds_dpT_Fe_Other->setDimension(1);
  ds_dpT_Fe_Other->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Fe_Other->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Fe_Other->setNormalizationValue(GetNormFactor(26));
  ds_dpT_Fe_Other->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Fe_Other);

  MinModDepCCQEXSec* ds_dpT_Pb = new MinModDepCCQEXSec("pTmu_Lead", 82, neutKE);
  ds_dpT_Pb->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Pb->setVariable(XSec::kPTLep);
  ds_dpT_Pb->setIsFluxIntegrated(true);
  ds_dpT_Pb->setDimension(1);
  ds_dpT_Pb->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Pb->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Pb->setNormalizationValue(GetNormFactor(82));
  ds_dpT_Pb->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Pb);

  MinModDepCCQEXSec* ds_dpT_Pb_QE = new MinModDepCCQEXSec("pTmu_Lead_QE", 82, neutKE, "QE");
  ds_dpT_Pb_QE->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Pb_QE->setVariable(XSec::kPTLep);
  ds_dpT_Pb_QE->setIsFluxIntegrated(true);
  ds_dpT_Pb_QE->setDimension(1);
  ds_dpT_Pb_QE->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Pb_QE->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Pb_QE->setNormalizationValue(GetNormFactor(82));
  ds_dpT_Pb_QE->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Pb_QE);

  MinModDepCCQEXSec* ds_dpT_Pb_RES = new MinModDepCCQEXSec("pTmu_Lead_RES", 82, neutKE, "RES");
  ds_dpT_Pb_RES->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Pb_RES->setVariable(XSec::kPTLep);
  ds_dpT_Pb_RES->setIsFluxIntegrated(true);
  ds_dpT_Pb_RES->setDimension(1);
  ds_dpT_Pb_RES->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Pb_RES->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Pb_RES->setNormalizationValue(GetNormFactor(82));
  ds_dpT_Pb_RES->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Pb_RES);

  MinModDepCCQEXSec* ds_dpT_Pb_DIS = new MinModDepCCQEXSec("pTmu_Lead_DIS", 82, neutKE, "DIS");
  ds_dpT_Pb_DIS->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Pb_DIS->setVariable(XSec::kPTLep);
  ds_dpT_Pb_DIS->setIsFluxIntegrated(true);
  ds_dpT_Pb_DIS->setDimension(1);
  ds_dpT_Pb_DIS->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Pb_DIS->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Pb_DIS->setNormalizationValue(GetNormFactor(82));
  ds_dpT_Pb_DIS->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Pb_DIS);

  MinModDepCCQEXSec* ds_dpT_Pb_2p2h = new MinModDepCCQEXSec("pTmu_Lead_2p2h", 82, neutKE, "2p2h");
  ds_dpT_Pb_2p2h->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Pb_2p2h->setVariable(XSec::kPTLep);
  ds_dpT_Pb_2p2h->setIsFluxIntegrated(true);
  ds_dpT_Pb_2p2h->setDimension(1);
  ds_dpT_Pb_2p2h->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Pb_2p2h->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Pb_2p2h->setNormalizationValue(GetNormFactor(82));
  ds_dpT_Pb_2p2h->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Pb_2p2h);

  MinModDepCCQEXSec* ds_dpT_Pb_Other = new MinModDepCCQEXSec("pTmu_Lead_Other", 82, neutKE, "Other");
  ds_dpT_Pb_Other->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Pb_Other->setVariable(XSec::kPTLep);
  ds_dpT_Pb_Other->setIsFluxIntegrated(true);
  ds_dpT_Pb_Other->setDimension(1);
  ds_dpT_Pb_Other->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Pb_Other->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Pb_Other->setNormalizationValue(GetNormFactor(82));
  ds_dpT_Pb_Other->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Pb_Other);

  MinModDepCCQEXSec* ds_dpT_Water = new MinModDepCCQEXSec("pTmu_Water", 8, neutKE);
  ds_dpT_Water->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Water->setVariable(XSec::kPTLep);
  ds_dpT_Water->setIsFluxIntegrated(true);
  ds_dpT_Water->setDimension(1);
  ds_dpT_Water->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Water->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Water->setNormalizationValue(GetNormFactor(8));
  ds_dpT_Water->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Water);

  MinModDepCCQEXSec* ds_dpT_Water_QE = new MinModDepCCQEXSec("pTmu_Water_QE", 8, neutKE, "QE");
  ds_dpT_Water_QE->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Water_QE->setVariable(XSec::kPTLep);
  ds_dpT_Water_QE->setIsFluxIntegrated(true);
  ds_dpT_Water_QE->setDimension(1);
  ds_dpT_Water_QE->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Water_QE->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Water_QE->setNormalizationValue(GetNormFactor(8));
  ds_dpT_Water_QE->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Water_QE);

  MinModDepCCQEXSec* ds_dpT_Water_RES = new MinModDepCCQEXSec("pTmu_Water_RES", 8, neutKE, "RES");
  ds_dpT_Water_RES->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Water_RES->setVariable(XSec::kPTLep);
  ds_dpT_Water_RES->setIsFluxIntegrated(true);
  ds_dpT_Water_RES->setDimension(1);
  ds_dpT_Water_RES->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Water_RES->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Water_RES->setNormalizationValue(GetNormFactor(8));
  ds_dpT_Water_RES->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Water_RES);

  MinModDepCCQEXSec* ds_dpT_Water_DIS = new MinModDepCCQEXSec("pTmu_Water_DIS", 8, neutKE, "DIS");
  ds_dpT_Water_DIS->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Water_DIS->setVariable(XSec::kPTLep);
  ds_dpT_Water_DIS->setIsFluxIntegrated(true);
  ds_dpT_Water_DIS->setDimension(1);
  ds_dpT_Water_DIS->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Water_DIS->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Water_DIS->setNormalizationValue(GetNormFactor(8));
  ds_dpT_Water_DIS->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Water_DIS);

  MinModDepCCQEXSec* ds_dpT_Water_2p2h = new MinModDepCCQEXSec("pTmu_Water_2p2h", 8, neutKE, "2p2h");
  ds_dpT_Water_2p2h->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Water_2p2h->setVariable(XSec::kPTLep);
  ds_dpT_Water_2p2h->setIsFluxIntegrated(true);
  ds_dpT_Water_2p2h->setDimension(1);
  ds_dpT_Water_2p2h->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Water_2p2h->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Water_2p2h->setNormalizationValue(GetNormFactor(8));
  ds_dpT_Water_2p2h->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Water_2p2h);

  MinModDepCCQEXSec* ds_dpT_Water_Other = new MinModDepCCQEXSec("pTmu_Water_Other", 8, neutKE, "Other");
  ds_dpT_Water_Other->setBinEdges(pt_nbins, pt_edges);
  ds_dpT_Water_Other->setVariable(XSec::kPTLep);
  ds_dpT_Water_Other->setIsFluxIntegrated(true);
  ds_dpT_Water_Other->setDimension(1);
  ds_dpT_Water_Other->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Water_Other->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Water_Other->setNormalizationValue(GetNormFactor(8));
  ds_dpT_Water_Other->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Water_Other);

  loop.runLoop();

  // Get the output histograms and save them to file
  string geniefilename =  "GENIEXSECEXTRACT_" + outName + "_neutronKE_" + to_string(neutKE) + ".root";
  TFile fout(geniefilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i)
  {
    loop.getXSecs()[i]->getXSecHist()->Write();
    loop.getXSecs()[i]->getEvRateHist()->Write();
  }

  return 0;
}
