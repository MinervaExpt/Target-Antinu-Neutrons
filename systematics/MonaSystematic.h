#ifndef MONASystematics_h
#define MONASystematics_h

//==============================================================================
// Get Several standard MINERvA systematics
//==============================================================================

#include "event/CVUniverse.h"
#include "PlotUtils/GenericVerticalSystematic.h"
#include "PlotUtils/NeutronInelasticReweighter.h"

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;

std::map<std::string,std::vector<int>> MonaMapDefault = {{"nGamma",{1000060120, 2112}},
							 {"threeAlpha",{1000020040, 100002040, 100002040, 2112}},
							 {"Bnp",{1000050110, 2112, 2212}}};

UniverseMap GetMonaSystematicMap(PlotUtils::ChainWrapper* chain)
{
  // return map
  UniverseMap error_bands;

  error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault)), 1.0));
  error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault)), -1.0));

  return error_bands;
}

#endif  // MONASystematics_h
