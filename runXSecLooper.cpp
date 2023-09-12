//TODO: ADD THE NORMALIZATION VALUE FROM THE TARGETS HERE. SELF-NORMALIZE.

#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
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
    nucleons += PlotUtils::TargetUtils::Get().GetPassiveTargetNNucleons(6, 1, true, 850.0);
  }
  else if (tgtZ == -1){
    nucleons += PlotUtils::TargetUtils::Get().GetTrackerNNucleons( 5980, 8422. , true, 850.0 ); //Hard-coding to see effect of requiring this compared to the normalization factor calculated by xSecLooper with it's own 
  }

  //double output = (nucleons > 0) ? trackerAtomsC/nucleons : 1.0;
  double output = (nucleons > 0) ? nAtoms_C/nucleons : 1.0;
  return output;
}


class MinModDepCCQEXSec : public XSec
{
public:
  MinModDepCCQEXSec(const char* name, const int tgtZ)
    :XSec(name), fTgtZ(tgtZ), fApothem(850.0)
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

    for (unsigned int i=0; i<nFSPart; ++i){
      int pdg = (int)chw.GetValue("mc_FSPartPDG",entry,i);
      double energy = (double)chw.GetValue("mc_FSPartE",entry,i);
      double proton_E = 1058.272;
      double neutron_E = 939.57+10.0;//Hard-coded neutKE for now.
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

  const double fApothem;

};

int main(const int argc, const char** argv)
{

  /*
  std::cout << GetNormFactor(-1) << std::endl;
  return 1;
  */

  //Read a playlist file from the command line
  if(argc < 4)
  {
    std::cerr << "Expected at least 3 command line argument, but got " << argc - 1 << ".\n\n"
              << "USAGE: runXSecLooper <outName> <playlist> <MCPlaylists.txt>\n\n"
              << "MCPlaylists.txt are a list of files which shall contain one .root file per line that has a Truth tree in it.\n"
              << "This program returns 0 when it suceeds.  It produces a .root file with GENIEXSECEXTRACT in its name.\n";
    return 1;
  }

  string outName = argv[1];
  TString pList = argv[2];

  std::vector<const char*> fileNames(argv+3,argv+argc);

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

  //Hard-coded to latest optimization for carbon migration matrix...
  double pt_edges_carbon[] = { 0.0, 0.09, 0.18, 0.25, 0.34, 0.425, 0.515, 0.64, 0.78, 0.94, 1.15, 1.5};
  int pt_nbins_carbon = 11;

  //Coded for testing against the latest before the Carbon bin optimization effort
  double pt_edges_std[] = { 0.0, 0.075, 0.15, 0.25, 0.325, 0.4, 0.475, 0.55, 0.7, 0.85, 1.0, 1.25, 1.5, 2.5};
  int pt_nbins_std = 13;

  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpT = new MinModDepCCQEXSec("pTmu_Tracker_std_binning", -1);
  ds_dpT->setBinEdges(pt_nbins_std, pt_edges_std);
  ds_dpT->setVariable(XSec::kPTLep);
  ds_dpT->setIsFluxIntegrated(true);
  ds_dpT->setDimension(1);
  ds_dpT->setFluxIntLimits(0.0, 100.0);
  ds_dpT->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT->setNormalizationValue(GetNormFactor(-1));
  ds_dpT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT);

  MinModDepCCQEXSec* ds_dpT_carbon_bins = new MinModDepCCQEXSec("pTmu_Tracker_carbon_binning", -1);
  ds_dpT_carbon_bins->setBinEdges(pt_nbins_carbon, pt_edges_carbon);
  ds_dpT_carbon_bins->setVariable(XSec::kPTLep);
  ds_dpT_carbon_bins->setIsFluxIntegrated(true);
  ds_dpT_carbon_bins->setDimension(1);
  ds_dpT_carbon_bins->setFluxIntLimits(0.0, 100.0);
  ds_dpT_carbon_bins->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_carbon_bins->setNormalizationValue(GetNormFactor(-1));
  ds_dpT_carbon_bins->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_carbon_bins);

  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpT_C = new MinModDepCCQEXSec("pTmu_C_std_binning", 6);
  ds_dpT_C->setBinEdges(pt_nbins_std, pt_edges_std);
  ds_dpT_C->setVariable(XSec::kPTLep);
  ds_dpT_C->setIsFluxIntegrated(true);
  ds_dpT_C->setDimension(1);
  ds_dpT_C->setFluxIntLimits(0.0, 100.0);
  ds_dpT_C->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_C->setNormalizationValue(GetNormFactor(6));
  ds_dpT_C->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_C);

  MinModDepCCQEXSec* ds_dpT_C_carbon_bins = new MinModDepCCQEXSec("pTmu_C_carbon_binning", 6);
  ds_dpT_C_carbon_bins->setBinEdges(pt_nbins_carbon, pt_edges_carbon);
  ds_dpT_C_carbon_bins->setVariable(XSec::kPTLep);
  ds_dpT_C_carbon_bins->setIsFluxIntegrated(true);
  ds_dpT_C_carbon_bins->setDimension(1);
  ds_dpT_C_carbon_bins->setFluxIntLimits(0.0, 100.0);
  ds_dpT_C_carbon_bins->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_C_carbon_bins->setNormalizationValue(GetNormFactor(6));
  ds_dpT_C_carbon_bins->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_C_carbon_bins);

  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpT_Fe = new MinModDepCCQEXSec("pTmu_Iron_std_binning", 26);
  ds_dpT_Fe->setBinEdges(pt_nbins_std, pt_edges_std);
  ds_dpT_Fe->setVariable(XSec::kPTLep);
  ds_dpT_Fe->setIsFluxIntegrated(true);
  ds_dpT_Fe->setDimension(1);
  ds_dpT_Fe->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Fe->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Fe->setNormalizationValue(GetNormFactor(26));
  ds_dpT_Fe->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Fe);

  MinModDepCCQEXSec* ds_dpT_Fe_carbon_bins = new MinModDepCCQEXSec("pTmu_Iron_carbon_binning", 26);
  ds_dpT_Fe_carbon_bins->setBinEdges(pt_nbins_carbon, pt_edges_carbon);
  ds_dpT_Fe_carbon_bins->setVariable(XSec::kPTLep);
  ds_dpT_Fe_carbon_bins->setIsFluxIntegrated(true);
  ds_dpT_Fe_carbon_bins->setDimension(1);
  ds_dpT_Fe_carbon_bins->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Fe_carbon_bins->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Fe_carbon_bins->setNormalizationValue(GetNormFactor(26));
  ds_dpT_Fe_carbon_bins->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Fe_carbon_bins);

  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpT_Pb = new MinModDepCCQEXSec("pTmu_Lead_std_binning", 82);
  ds_dpT_Pb->setBinEdges(pt_nbins_std, pt_edges_std);
  ds_dpT_Pb->setVariable(XSec::kPTLep);
  ds_dpT_Pb->setIsFluxIntegrated(true);
  ds_dpT_Pb->setDimension(1);
  ds_dpT_Pb->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Pb->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Pb->setNormalizationValue(GetNormFactor(82));
  ds_dpT_Pb->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Pb);

  MinModDepCCQEXSec* ds_dpT_Pb_carbon_bins = new MinModDepCCQEXSec("pTmu_Lead_carbon_binning", 82);
  ds_dpT_Pb_carbon_bins->setBinEdges(pt_nbins_carbon, pt_edges_carbon);
  ds_dpT_Pb_carbon_bins->setVariable(XSec::kPTLep);
  ds_dpT_Pb_carbon_bins->setIsFluxIntegrated(true);
  ds_dpT_Pb_carbon_bins->setDimension(1);
  ds_dpT_Pb_carbon_bins->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Pb_carbon_bins->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Pb_carbon_bins->setNormalizationValue(GetNormFactor(82));
  ds_dpT_Pb_carbon_bins->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Pb_carbon_bins);

  // Flux-integrated over the range 0.0 to 100.0 GeV
  MinModDepCCQEXSec* ds_dpT_Water = new MinModDepCCQEXSec("pTmu_Water_std_binning", 8);
  ds_dpT_Water->setBinEdges(pt_nbins_std, pt_edges_std);
  ds_dpT_Water->setVariable(XSec::kPTLep);
  ds_dpT_Water->setIsFluxIntegrated(true);
  ds_dpT_Water->setDimension(1);
  ds_dpT_Water->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Water->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Water->setNormalizationValue(GetNormFactor(8));
  ds_dpT_Water->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Water);

  MinModDepCCQEXSec* ds_dpT_Water_carbon_bins = new MinModDepCCQEXSec("pTmu_Water_carbon_binning", 8);
  ds_dpT_Water_carbon_bins->setBinEdges(pt_nbins_carbon, pt_edges_carbon);
  ds_dpT_Water_carbon_bins->setVariable(XSec::kPTLep);
  ds_dpT_Water_carbon_bins->setIsFluxIntegrated(true);
  ds_dpT_Water_carbon_bins->setDimension(1);
  ds_dpT_Water_carbon_bins->setFluxIntLimits(0.0, 100.0);
  ds_dpT_Water_carbon_bins->setNormalizationType(XSec::kSelfNorm);  
  ds_dpT_Water_carbon_bins->setNormalizationValue(GetNormFactor(8));
  ds_dpT_Water_carbon_bins->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_Water_carbon_bins);

  loop.runLoop();

  // Get the output histograms and save them to file
  string geniefilename =  "GENIEXSECEXTRACT_" + outName + ".root";
  TFile fout(geniefilename.c_str(), "RECREATE");
  for(uint i=0; i<loop.getXSecs().size(); ++i)
  {
    loop.getXSecs()[i]->getXSecHist()->Write();
    loop.getXSecs()[i]->getEvRateHist()->Write();
  }

  return 0;
}
