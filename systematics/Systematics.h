#ifndef Systematics_h
#define Systematics_h

//==============================================================================
// Get Several standard MINERvA systematics
//==============================================================================

#include "event/CVUniverse.h"
#include "systematics/MonaSystematic.h"
#include "PlotUtils/FluxSystematics.h"
#include "PlotUtils/GenieSystematics.h"
#include "PlotUtils/MinosEfficiencySystematics.h"
#include "PlotUtils/MuonSystematics.h"
#include "PlotUtils/MnvTuneSystematics.h"
#include "PlotUtils/MuonResolutionSystematics.h"
#include "PlotUtils/AngleSystematics.h"
#include "PlotUtils/ResponseSystematics.h"

typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;

UniverseMap removeElasticBand(UniverseMap unprunedBands){
  UniverseMap prunedBands;
  for (const auto& band: unprunedBands){
    size_t pos=0;
    std::cout << "Systematic band from GENIE: " << band.first << std::endl;
    std::cout << "Not Found? " << ((pos=(band.first).find("FrElas_N"))==std::string::npos) << std::endl;
    if ((pos=(band.first).find("FrElas_N")) == std::string::npos) prunedBands[band.first] = band.second;
    else std::cout << "Skipping band broken by FSI bug fix" << std::endl;
  }
  return prunedBands;
}

UniverseMap GetStandardSystematics(PlotUtils::ChainWrapper* chain, std::string response_name_tag = "", bool doMona = true, bool pruneElasN = false)
{
  // return map
  UniverseMap error_bands;

  // CV
  error_bands[std::string("cv")].push_back(new CVUniverse(chain));

  //========================================================================
  // FLUX
  //========================================================================
  UniverseMap bands_flux =
      PlotUtils::GetFluxSystematicsMap<CVUniverse>(chain, CVUniverse::GetNFluxUniverses());
  error_bands.insert(bands_flux.begin(), bands_flux.end());

  //========================================================================
  // GENIE
  //========================================================================
  // Standard
  UniverseMap bands_genie =
      PlotUtils::GetGenieSystematicsMap<CVUniverse>(chain); //PlotUtils::GetStandardGenieSystematicsMap<CVUniverse>(chain);
  if (pruneElasN){
    std::cout << "Entering pruning" << std::endl;
    bands_genie = removeElasticBand(bands_genie);
    std::cout << "After pruning" << std::endl;
  }
  error_bands.insert(bands_genie.begin(), bands_genie.end());

  //========================================================================
  // MnvTunes
  //========================================================================
  // RPA
  UniverseMap bands_rpa = PlotUtils::GetRPASystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_rpa.begin(), bands_rpa.end());

  // 2P2H
  UniverseMap bands_2p2h = PlotUtils::Get2p2hSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_2p2h.begin(), bands_2p2h.end());

  //========================================================================
  // Muons
  //========================================================================
  // Muon reco in MINERvA -- Catchall systematic for pmu reco in minerva.
  // Lateral-only. Shifts pmu.
  UniverseMap bands_muon_minerva =
      PlotUtils::GetMinervaMuonSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_minerva.begin(), bands_muon_minerva.end());

  // Muons in MINOS -- Catchall systematic for wiggle solution -- correlates
  // flux universes and minos muon momentum reco.
  // Lateral AND Vertical systematic. Shifts Pmu and GetFluxAndCVUniverse.
  //
  // Expect a non-zero systematic even when no pmu involved.
  UniverseMap bands_muon_minos =
     PlotUtils::GetMinosMuonSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_minos.begin(), bands_muon_minos.end());

  // Vertical only
  UniverseMap bands_minoseff =
      PlotUtils::GetMinosEfficiencySystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_minoseff.begin(), bands_minoseff.end());

  UniverseMap bands_muon_resolution = PlotUtils::GetMuonResolutionSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_muon_resolution.begin(), bands_muon_resolution.end());

  UniverseMap bands_geant = PlotUtils::GetGeantHadronSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_geant.begin(), bands_geant.end());

  // Beam angle
  UniverseMap bands_angle = PlotUtils::GetAngleSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_angle.begin(), bands_angle.end());

  //Andrew's MONA Neutron Inelastic Reweight
  if (doMona){
    UniverseMap bands_mona = GetMonaSystematicMap(chain);
    error_bands.insert(bands_mona.begin(), bands_mona.end());
  }

  // Hadron inelastics cross sections
  //TODO: There's some special recoil function I need to write for the response systematics to work correctly
  UniverseMap bands_response = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain, response_name_tag);
  //UniverseMap bands_response = PlotUtils::GetResponseSystematicsMap<CVUniverse>(chain);
  error_bands.insert(bands_response.begin(), bands_response.end());

  return error_bands;
}

#endif  // Systematics_h
