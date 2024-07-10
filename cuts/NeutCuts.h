//File: NeutCuts
//Brief: Neut Cuts for neutron selection initially derived from Tejin's cuts.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef David_NeutCuts_H
#define David_NeutCuts_H

#include "PlotUtils/Cut.h"
#include "util/GetRecoTargetZ.h"

namespace MyNeutCuts{

  template <class UNIVERSE, class EVENT>
  class LeadNeutIs3D: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor                                                          
    LeadNeutIs3D():
    PlotUtils::Cut<UNIVERSE, EVENT>("Lead Neut. Cand. is 3D")
    {
    }

  private:
  // THE cut function                                                     
    bool checkCut(const UNIVERSE& /*univ*/, EVENT& evt) const override
    {
      // Call a CVUniverse member function to make the cut                
      return evt.GetLeadingNeutCand().GetIs3D() == 1;
    }
    
  };

  //Currently a fixed angle. Could make it something
  template <class UNIVERSE, class EVENT>
  class LeadNeutIsFarFromMuon: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    LeadNeutIsFarFromMuon(): PlotUtils::Cut<UNIVERSE, EVENT>("Lead Neut. Cand. > 15deg from muon") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& evt) const override
    {
      TVector3 muonMom(univ.GetMuon4V().X(),univ.GetMuon4V().Y(),univ.GetMuon4V().Z());
      TVector3 leadFP=evt.GetLeadingNeutCand().GetFlightPath();
      if (leadFP.Mag() == 0 || muonMom.Mag() == 0) return false;
      else return leadFP.Angle(muonMom) > 0.261799388;
    }
  };

  //Fixed for now. Can make variable later
  template <class UNIVERSE, class EVENT>
  class LeadNeutZDistMin: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    LeadNeutZDistMin(): PlotUtils::Cut<UNIVERSE, EVENT>("Lead Neut. Cand. far from vtx. in z") {}

  private:
    bool checkCut(const UNIVERSE& /*univ*/, EVENT& evt) const override
    {
      return evt.GetLeadingNeutCand().GetVtxZDist() >= 100;
    }
  };

  //Fixed for now. Can make variable later
  template <class UNIVERSE, class EVENT>
  class LeadNeutDistMin: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    LeadNeutDistMin(): PlotUtils::Cut<UNIVERSE, EVENT>("Lead Neut. Cand. far from vtx.") {}

  private:
    bool checkCut(const UNIVERSE& /*univ*/, EVENT& evt) const override
    {
      return evt.GetLeadingNeutCand().GetVtxDist() >= 100;
    }
  };

  template <class UNIVERSE, class EVENT>
  class LeadNeutOutsideTgt: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    LeadNeutOutsideTgt(): PlotUtils::Cut<UNIVERSE, EVENT>("Lead Neut. Cand. Outside Tgt.") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& evt) const override
    {
      double x = evt.GetLeadingNeutCand().GetXPos();
      double y = evt.GetLeadingNeutCand().GetYPos();
      double z = evt.GetLeadingNeutCand().GetZPos();
      int TgtByZ = util::GetRecoTargetZWide(x, y, z, false);
      return TgtByZ < 1 || TgtByZ > 6;
    }

  };

  template <class UNIVERSE, class EVENT>
  class LeadNeutInTracker: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    LeadNeutInTracker(double boundary): PlotUtils::Cut<UNIVERSE, EVENT>("Lead Neut. Cand. far from vtx."), fBound(boundary) {}

  private:
    bool checkCut(const UNIVERSE& /*univ*/, EVENT& evt) const override
    {
      return evt.GetLeadingNeutCand().GetZPos() >= fBound;
    }
    
    double fBound;

  };

}
#endif //David_BlobCuts_H
