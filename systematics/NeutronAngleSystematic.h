//This is meant to hold what it needs to run a universe which weights to the ratio of data/MC in bins of pT and neutron angle. Need to think about how exactly to get that ratio in a useful fashion.
//A Histogram to sample should be fine... just need to make some decisions about the definition of the angular weight.

#include "event/CVUniverse.h"

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;

class NeutronAngle: public CVUniverse{
public:
NeutronAngle(PlotUtils::ChainWrapper* chw): CVUniverse(chw)
  {
  }
  
  virtual ~NeutronAngle() = default;
  
  std::string ShortName() const override
  {
    return "NeutronAngle";
  }
  
  std::string LatexName() const override
  {
    return "Neutron Angle";
  }
};

UniverseMap GetNeutronDroppingUnivs(PlotUtils::ChainWrapper* chw)
{
  UniverseMap error_bands;
  error_bands["NeutronAngle"].push_back(new NeutronAngle(chw));

  return error_bands;
};
