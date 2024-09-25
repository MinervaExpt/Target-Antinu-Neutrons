//File: CCQECuts
//Brief: Cuts for CCQE selection derived from Tejin's cuts.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef David_CCQECuts_H
#define David_CCQECuts_H

#include "PlotUtils/Cut.h"

namespace MyCCQECuts{

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class PMuRange: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor                                                          
    PMuRange(const std::string& name, const double pMin, const double pMax):
    PlotUtils::Cut<UNIVERSE, EVENT>(name), fMin(pMin), fMax(pMax)
    {
    }

  private:
  // THE cut function                                                     
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      // Call a CVUniverse member function to make the cut                
      return univ.GetMuonP() >= fMin && univ.GetMuonP() <= fMax;
    }
    
    const double fMin;
    const double fMax;

  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class CorrectDaisy: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor                                                          
    CorrectDaisy(const std::string& name, const int petal):
    PlotUtils::Cut<UNIVERSE, EVENT>(name), fPetal(petal)
    {
    }

  private:
  // THE cut function                                                     
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      // Call a CVUniverse member function to make the cut                
      int petal = univ.GetRecoDaisyPetal();
      return (petal == fPetal);
    }
    
    const int fPetal;

  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class IsAntiNu: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    IsAntiNu(): PlotUtils::Cut<UNIVERSE, EVENT>("Helicity == 2") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetNuHelicity() == 2;
    }
  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class IsNu: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    IsNu(): PlotUtils::Cut<UNIVERSE, EVENT>("Helicity == 1") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetNuHelicity() == 1;
    }
  };

  //May want to at some point make this a maximum so that it could
  //be tested or a sideband...
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class IsSingleTrack: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    IsSingleTrack(): PlotUtils::Cut<UNIVERSE, EVENT>("Single Track") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetNTracks() == 1;
    }
  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class NoMichels: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    NoMichels(): PlotUtils::Cut<UNIVERSE, EVENT>("No Improved Michels") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetNImprovedMichel() <= 0;
    }
  };

  // Maybe break this apart/ use the EVENT class to handle this cut. Currently implementing as before.
  template <class UNIVERSE, class EVENT>
  class AllEMBlobsCuts: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    AllEMBlobsCuts(const bool useEvent): PlotUtils::Cut<UNIVERSE, EVENT>("EM Blobs Cuts"), fUse(useEvent) {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& evt) const override
    {
      if (!fUse){
	std::vector<double> EMBlobInfo = univ.GetEMNBlobsTotalEnergyTotalNHits();
	if (EMBlobInfo.at(0) == 0) return true;
	return EMBlobInfo.at(0) < 2 && EMBlobInfo.at(1) >= 10.0*EMBlobInfo.at(2);
      }
      else{
	double nBlobs = evt.GetEMNBlobs();
	double BlobE = evt.GetEMBlobE();
	double BlobNHit = evt.GetEMBlobNHits();
	if (nBlobs == 0) return true;
	return nBlobs < 2 && BlobE >= 10.0*BlobNHit;
      }
    }

    const bool fUse;

  };

  template <class UNIVERSE, class EVENT>
  class EMNBlobsCut: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    EMNBlobsCut(const bool useEvent): PlotUtils::Cut<UNIVERSE, EVENT>("EM N Blobs Cut"), fUse(useEvent) {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& evt) const override
    {
      if (!fUse){
	std::vector<double> EMBlobInfo = univ.GetEMNBlobsTotalEnergyTotalNHits();
	return EMBlobInfo.at(0) < 2;
      }
      else{
	double nBlobs = evt.GetEMNBlobs();
	return nBlobs < 2;
      }
    }

    const bool fUse;

  };

  template <class UNIVERSE, class EVENT>
  class EMBlobENHitRatioCut: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    EMBlobENHitRatioCut(const bool useEvent): PlotUtils::Cut<UNIVERSE, EVENT>("EM Blob E/NHit Cut"), fUse(useEvent) {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& evt) const override
    {
      if (!fUse){
	std::vector<double> EMBlobInfo = univ.GetEMNBlobsTotalEnergyTotalNHits();
	if (EMBlobInfo.at(0) == 0) return true;
	return EMBlobInfo.at(1) >= 10.0*EMBlobInfo.at(2);
      }
      else{
	double nBlobs = evt.GetEMNBlobs();
	double BlobE = evt.GetEMBlobE();
	double BlobNHit = evt.GetEMBlobNHits();
	if (nBlobs == 0) return true;
	return BlobE >= 10.0*BlobNHit;
      }
    }

    const bool fUse;

  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class RecoilCut: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    RecoilCut(): PlotUtils::Cut<UNIVERSE, EVENT>("Recoil Cut") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      double Q2GeV = univ.GetQ2QEPickledGeV();
      double recoilEGeV = univ.GetDANRecoilEnergyGeV();
      if (Q2GeV < 0.0 || recoilEGeV < 0.0) return false;
      else if (Q2GeV < 0.175) return (recoilEGeV < 0.13);
      //else if (Q2GeV < 0.3) return (recoilEGeV < (0.04+0.43*Q2GeV)); Tejin modification. Removed since it's hydrogen specific. Need to revisit this optimization at some point though.
      else if (Q2GeV < 1.4) return (recoilEGeV < (0.08+0.3*Q2GeV));
      else return (recoilEGeV < 0.5);
    }

  };

  //Looser recoil cut useful for the definition of a recoil sideband.
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class LooseRecoilCut: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
  LooseRecoilCut(double offset): PlotUtils::Cut<UNIVERSE, EVENT>("Loose (SB) Recoil Cut"), fOffset(offset) {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      double Q2GeV = univ.GetQ2QEPickledGeV();
      double recoilEGeV = univ.GetDANRecoilEnergyGeV();
      if (Q2GeV < 0.0 || recoilEGeV < 0.0) return false;
      else if (Q2GeV < 0.175) return (recoilEGeV < 0.13+fOffset);
      //else if (Q2GeV < 0.3) return (recoilEGeV < (0.04+0.43*Q2GeV)); Tejin modification. Removed since it's hydrogen specific. Need to revisit this optimization at some point though.
      else if (Q2GeV < 1.4) return (recoilEGeV < (0.08+fOffset+0.3*Q2GeV));
      else return (recoilEGeV < 0.5+fOffset);
    }

    const double fOffset;

  };

  //Looser recoil cut useful for the definition of a recoil sideband.
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class RemoveRecoilBand: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
  RemoveRecoilBand(double offset): PlotUtils::Cut<UNIVERSE, EVENT>("Remove Band Recoil Cut"), fOffset(offset) {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      double Q2GeV = univ.GetQ2QEPickledGeV();
      double recoilEGeV = univ.GetDANRecoilEnergyGeV();
      if (Q2GeV < 0.0 || recoilEGeV < 0.0) return false;
      else if (Q2GeV < 0.175) return (recoilEGeV < 0.13 || recoilEGeV >= 0.13+fOffset);
      //else if (Q2GeV < 0.3) return (recoilEGeV < (0.04+0.43*Q2GeV)); Tejin modification. Removed since it's hydrogen specific. Need to revisit this optimization at some point though.
      else if (Q2GeV < 1.4) return (recoilEGeV < (0.08+0.3*Q2GeV) || recoilEGeV >= (0.08+fOffset+0.3*Q2GeV));
      else return (recoilEGeV < 0.5 || recoilEGeV >= 0.5+fOffset);
    }

    const double fOffset;

  };

}
#endif //David_CCQECuts_H
