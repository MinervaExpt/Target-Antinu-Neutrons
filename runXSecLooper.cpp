//TODO: ADD THE NORMALIZATION VALUE FROM THE TARGETS HERE. SELF-NORMALIZE.

#include "GENIEXSecExtract/XSecLooper.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

#include "util/GetRecoTargetZ.h"

#include <cstdlib>
typedef unsigned int uint;

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

    if (fTgtZ != -1){
      if (fTgtZ != trueTgtZ) return false;
      int trueTgtCode = util::GetTrueTgtCode(trueTgtZ, x, y, z);//Is this problematic to use the exact same definition? Feels self-fulfilling...
      return (trueTgtCode > 0 && trueTgtCode != 6666); //Currently ignoring water still.
    }
    else {
      return (z > 5980 && z < 8422);
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
  //Read a playlist file from the command line
  if(argc < 3)
  {
    std::cerr << "Expected at least 2 command line argument, but got " << argc - 1 << ".\n\n"
              << "USAGE: runXSecLooper <outName> <MCPlaylists.txt>\n\n"
              << "MCPlaylists.txt are a list of files which shall contain one .root file per line that has a Truth tree in it.\n"
              << "This program returns 0 when it suceeds.  It produces a .root file with GENIEXSECEXTRACT in its name.\n";
    return 1;
  }

  string outName = argv[1];

  std::vector<const char*> fileNames(argv+2,argv+argc);

  // Create the XSecLooper and tell it the input files
  // Inputs should be the merged ntuples:
  XSecLooper loop(fileNames[0]);
  for(auto whichFile = fileNames.begin()+1; whichFile != fileNames.end(); ++whichFile)
    {
      loop.addFiles(*whichFile);
    }


  // Tell the XSecLooper which neutrino type we're considering (mandatory)
  loop.setNuPDG(-14);
  loop.setPlaylist(PlotUtils::FluxReweighter::minervame6A);

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
  ds_dpT->setNormalizationType(XSec::kPerNucleon);  
  ds_dpT->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT);

  MinModDepCCQEXSec* ds_dpT_carbon = new MinModDepCCQEXSec("pTmu_Tracker_carbon_binning", -1);
  ds_dpT_carbon->setBinEdges(pt_nbins_carbon, pt_edges_carbon);
  ds_dpT_carbon->setVariable(XSec::kPTLep);
  ds_dpT_carbon->setIsFluxIntegrated(true);
  ds_dpT_carbon->setDimension(1);
  ds_dpT_carbon->setFluxIntLimits(0.0, 100.0);
  ds_dpT_carbon->setNormalizationType(XSec::kPerNucleon);  
  ds_dpT_carbon->setUniverses(0); //default value, put 0 if you do not want universes to be included.
  loop.addXSec(ds_dpT_carbon);

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
