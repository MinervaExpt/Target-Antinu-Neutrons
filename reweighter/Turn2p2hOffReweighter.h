//File: Turn2p2hOffReweighter.h
//Brief: Turn off all 2p2h events in the regular sample so that I can add in a signal-only sample to replace them.
//Author: Andrew Olivier aolivier@ur.rochester.edu, David Last

#ifndef TURN2P2HOFFREWEIGHTER_CPP
#define TURN2P2HOFFREWEIGHTER_CPP

//PlotUtils includes
#include "PlotUtils/NSFDefaults.h"
#include "PlotUtils/MnvTuneSystematics.h"

//Reweighter includes
#include "PlotUtils/Reweighter.h"

template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class Turn2p2hOffReweighter: public PlotUtils::Reweighter<UNIVERSE, EVENT>
{
  public:
    Turn2p2hOffReweighter(): PlotUtils::Reweighter<UNIVERSE, EVENT>()
    {
    }

    virtual ~Turn2p2hOffReweighter() = default;

    double GetWeight(const UNIVERSE& univ, const EVENT& /*event*/) const override
    {
      return (univ.GetInt("mc_intType")==8)?0:1;
    }

    std::string GetName() const override { return "Turn2p2hOff"; }
    bool DependsReco() const override { return false; }
};

#endif //TURN2P2HOFFREWEIGHTER_CPP
