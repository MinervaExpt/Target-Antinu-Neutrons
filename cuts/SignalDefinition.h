//File: SignalDefinition
//Brief: Signal Definition for CCQE 1+ selection derived from standard anti-nu CCQE definitions.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef David_CCQESigDef_H
#define David_CCQESigDef_H

#include "PlotUtils/Cut.h"

namespace MySignal{

  template <class UNIVERSE>
  class IsAntiNu: public PlotUtils::SignalConstraint<UNIVERSE>
  {
  public:
    IsAntiNu(): PlotUtils::SignalConstraint<UNIVERSE>("Neutrino PDG = -14") {}

  private:
    bool checkConstraint(const UNIVERSE& univ) const //override
    {
      return univ.GetTruthNuPDG() == -14;
    }
  };

  template <class UNIVERSE>
  class IsCorrectFS: public PlotUtils::SignalConstraint<UNIVERSE>
  {
  public:
    IsCorrectFS(): PlotUtils::SignalConstraint<UNIVERSE>("CCQE 1+ neutron FS") {}

  private:
    bool checkConstraint(const UNIVERSE& univ) const //override
    {
      int genie_n_muons = 0;
      int genie_n_mesons = 0;
      int genie_n_heavy_baryons_plus_pi0s = 0;
      int genie_n_photons =0;
      int genie_n_protons = 0;
      int genie_n_neutrons = 0;
      
      std::vector<int> PDGs = univ.GetFSPartPDG();
      std::vector<double> Es = univ.GetFSPartE();

      for (unsigned int i=0; i<PDGs.size(); ++i){
	int pdg = PDGs.at(i);
	double energy = Es.at(i);
	double proton_E = 1058.272;
	double neutron_E = 949.57;
	if (abs(pdg) == 13) genie_n_muons++;
	else if ( pdg == 22  && energy > 10) genie_n_photons++;
	else if ( abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 323 || pdg == 111 || pdg == 130 || pdg == 310 || pdg == 311 || pdg == 313 ){
	  genie_n_mesons++;
	}
	else if ( pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 4112 || pdg == 4122 || pdg == 4222 || pdg == 411 || pdg == 421 || pdg == 111){
	  genie_n_heavy_baryons_plus_pi0s++;
	}
	else if ( pdg == 2212 && energy > proton_E ) genie_n_protons++;
	else if ( pdg == 2112 && energy > neutron_E) genie_n_neutrons++;
      }

      return genie_n_muons == 1 &&
	genie_n_mesons == 0 &&
	genie_n_heavy_baryons_plus_pi0s == 0 &&
	genie_n_photons == 0 &&
	genie_n_protons == 0 &&
	genie_n_neutrons > 0;
    }
  };

}
#endif //David_CCQESigDef_H
