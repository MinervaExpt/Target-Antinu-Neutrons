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
							 {"threeAlpha",{1000020040, 1000020040, 1000020040, 2112}},
							 {"Bnp",{1000050110, 2112, 2212}}};


std::map<std::string, double> MonaMapCVNSigma = {{"nGamma",0.0},
						 {"threeAlpha",0.0},
						 {"Bnp",0.0}};

std::map<std::string, double> MonaMapBnpUp = {{"nGamma",0.0},
					      {"threeAlpha",0.0},
					      {"Bnp",1.0}};

std::map<std::string, double> MonaMapBnpDown = {{"nGamma",0.0},
						{"threeAlpha",0.0},
						{"Bnp",-1.0}};

std::map<std::string, double> MonaMapNGammaUp = {{"nGamma",1.0},
						 {"threeAlpha",0.0},
						 {"Bnp",0.0}};

std::map<std::string, double> MonaMapNGammaDown = {{"nGamma",-1.0},
						   {"threeAlpha",0.0},
						   {"Bnp",0.0}};

std::map<std::string, double> MonaMap3AlphaUp = {{"nGamma",0.0},
						 {"threeAlpha",1.0},
						 {"Bnp",0.0}};

std::map<std::string, double> MonaMap3AlphaDown = {{"nGamma",0.0},
						   {"threeAlpha",-1.0},
						   {"Bnp",0.0}};

UniverseMap GetMonaSystematicMap(PlotUtils::ChainWrapper* chain)
{
  // return map
  UniverseMap error_bands;

  //error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault)), 1.0));
  //error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault)), -1.0));

  error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,MonaMapCVNSigma,1)), 1.0));

  /*
  error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,MonaMapBnpUp,1)), 1.0));
  error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,MonaMapBnpDown,1)), 1.0));
  error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,MonaMapNGammaUp,1)), 1.0));
  error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,MonaMapNGammaDown,1)), 1.0));
  error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,MonaMap3AlphaUp,1)), 1.0));
  error_bands["NeutronInelasticsReweight"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,MonaMap3AlphaDown,1)), 1.0));
  */  
  //Checking that these bands give the same result as the single universe band above as they should... I can perform this check even with the modified systematics since the mode in the reweighter is currently set to be the same.
  //error_bands["NeutronInelasticsReweight_1"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,1)), 1.0, "1"));
  //error_bands["NeutronInelasticsReweight_1"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,1)), -1.0, "1"));

  /*
  error_bands["NeutronInelasticsReweight_2"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,2)), 1.0, "2"));
  error_bands["NeutronInelasticsReweight_2"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,2)), -1.0, "2"));

  error_bands["NeutronInelasticsReweight_3"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,3)), 1.0, "3"));
  error_bands["NeutronInelasticsReweight_3"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,3)), -1.0, "3"));

  error_bands["NeutronInelasticsReweight_4"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,4)), 1.0, "4"));
  error_bands["NeutronInelasticsReweight_4"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,4)), -1.0, "4"));

  error_bands["NeutronInelasticsReweight_5"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,5)), 1.0, "5"));
  error_bands["NeutronInelasticsReweight_5"].push_back(new PlotUtils::GenericVerticalUniverse<CVUniverse, PlotUtils::detail::empty>(chain, std::unique_ptr<PlotUtils::Reweighter<CVUniverse, PlotUtils::detail::empty>>(new NeutronInelasticReweighter<CVUniverse>(MonaMapDefault,5)), -1.0, "5"));
  */
  return error_bands;
}

#endif  // MONASystematics_h
