//studies includes
#include "studies/Study.h"

//David's includes
#include "util/Categorized.h"
#include "event/CVUniverse.h"

class NeutronInelReweightStudy: public Study
{
  public:
    //PerMichelVarByGENIELabel fills a histogram with 1 entry per Michel with some variable calculated from that Michel.  Your function will get to see the CVUniverse, the NeutronEvent (= reconstructed Michels), and which Michel it's looping over.
  NeutronInelReweightStudy(const std::map<std::string, std::vector<CVUniverse*>>& univs, double neutKE=0.0): Study(), fNeutKE(neutKE)
    {
      // Change this to be the right set of categories breaking down signal/BKG appropriately.
      // Note that here the signal definition doesn't have
      std::map<int, std::string> Categs = {{11, "sigQE"},
					   {18, "sig2p2h"},
					   {19, "sigOther"},
					   {1, "1chargePi"},
					   {2, "1neutPi"},
					   {3, "NPi"},
					   {4, "SubThresh"},
					   {5, "TrackableProt"}};
      
      // Do whatever I need to to make this work for 2D histos
      // Change it to a fixed binning start fine and move to coarse as needed
      // My range is 0-1.5, so only worry about it there so much for now.
      // For the lead neutron range 0-1GeV should cover most, if not, everything

      const int nBinsX = 200;
      const int nBinsY = 200;
      const double hiBinX = 1.5;
      const double hiBinY = 1000;
      
      m_pT_v_LeadTn_ByCateg = new util::Categorized<HIST, int>("TrueEvRate_pT_v_LeadTn", "; pT [GeV/c]; Lead T_{n} [MeV]", Categs, nBinsX, 0.0, hiBinX, nBinsY, 0.0, hiBinY, univs);

    }

    void SaveOrDraw(TDirectory& outDir)
    {
      TDirectory* dir;
      dir = outDir.GetDirectory("NeutronInelHists");
      if (dir == NULL){
	outDir.mkdir("NeutronInelHists");
      }
      dir = outDir.GetDirectory("NeutronInelHists");
      
      //MAke sure this makes sense... Maybe need to match Variable changes
      m_pT_v_LeadTn_ByCateg->visit([dir](HIST& wrapper)
                                {
                                  wrapper.SyncCVHistos();
                                  wrapper.hist->SetDirectory(dir);
                                });

    }

    void SaveOrDrawData(TFile& outFile){
      return;
    }
  
    void SaveOrDrawMC(TFile& outFile){
      return;
    }
  
  private:
    using HIST = PlotUtils::Hist2DWrapper<CVUniverse>;

    util::Categorized<HIST, int>* m_pT_v_LeadTn_ByCateg;

    double fNeutKE;
  
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) { return; }

    //All of your plots happen here so far.
    void fillSelectedSignal(const CVUniverse& univ, const NeutronEvent& evt, const double weight) { return; }

    //Do nothing for now...  Good place for efficiency denominators in the future.
    void fillTruthSignal(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      double pT = univ.GetMuonPTTrue();
      double leadTn = univ.GetMaxFSNeutronKE();
      int label = univ.GetNeutronReweightCategory(fNeutKE);
      //int label = 1;
      (*m_pT_v_LeadTn_ByCateg)[label].FillUniverse(&univ, pT, leadTn, weight);
      return;
    }
};
