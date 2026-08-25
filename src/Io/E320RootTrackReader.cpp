#include "TrackingPipeline/Io/E320RootTrackReader.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "TFile.h"
#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"
#include "TrackingPipeline/MagneticField/ConstantMagField.hpp"
#include "TrackingPipeline/MagneticField/IdealQuadrupoleMagField.hpp"

using namespace Acts::UnitLiterals;

E320::E320RootTrackReader::E320RootTrackReader(const Config& config,
                                               Acts::Logging::Level level)
    : IReader(),
      m_cfg(config),
      m_logger(Acts::getDefaultLogger(name(), level)) {
  if (m_cfg.filePaths.empty()) {
    throw std::invalid_argument("Missing filename");
  }
  if (m_cfg.treeName.empty()) {
    throw std::invalid_argument("Missing tree name");
  }

  if (m_cfg.filePaths.size() == 1) {
    m_file = new TFile(m_cfg.filePaths.at(0).c_str());
    m_tree = m_file->Get<TTree>(m_cfg.treeName.c_str());
  } else {
    m_chainOwner = new TChain(m_cfg.treeName.c_str());
    // Add the files to the chain
    for (const auto& path : m_cfg.filePaths) {
      m_chainOwner->Add(path.c_str());
    }
    m_tree = dynamic_cast<TTree*>(m_chainOwner);
  }

  // Set event Id branch
  m_tree->SetBranchAddress("eventId", &m_eventId);
  if (m_tree->GetBranch("eventId") == nullptr) {
    throw std::invalid_argument("Missing eventId branch");
  }
  auto nEntries = static_cast<std::size_t>(m_tree->GetEntries());

  if (!m_cfg.mergeIntoOneEvent) {
    // Add the first entry
    m_tree->GetEntry(0);
    m_eventMap.emplace_back(m_eventId, 0, 0);

    for (std::size_t i = 0; i < nEntries; ++i) {
      m_tree->GetEntry(i);
      if (m_eventId != std::get<0>(m_eventMap.back())) {
        std::get<2>(m_eventMap.back()) = i;
        m_eventMap.emplace_back(m_eventId, i, i);
      }
      if (i == nEntries - 1) {
        std::get<2>(m_eventMap.back()) = nEntries;
      }
    }
  } else {
    m_eventMap = {{0, 0, nEntries}};
  }

  // Sort by event id
  std::ranges::sort(m_eventMap, [](const auto& a, const auto& b) {
    return std::get<0>(a) < std::get<0>(b);
  });

  ACTS_DEBUG("Event range: " << availableEvents().first << " - "
                             << availableEvents().second);

  //------------------------------------------------------------------
  // Set the rest of the branches

  // Event meta data
  m_tree->SetBranchAddress("eudaqTrgN", &m_eudaqTrgN);
  m_tree->SetBranchAddress("eudaqDAQNumber", &m_eudaqDAQNumber);
  m_tree->SetBranchAddress("eudaqRunStartTs", &m_eudaqRunStartTs);
  m_tree->SetBranchAddress("eudaqRunEndTs", &m_eudaqRunEndTs);
  m_tree->SetBranchAddress("epicsParity", &m_epicsParity);
  m_tree->SetBranchAddress("epicsPulseId", &m_epicsPulseId);
  m_tree->SetBranchAddress("epicsDAQNumber", &m_epicsDAQNumber);

  // Magnet configuration as seen by track
  m_tree->SetBranchAddress("quad1Grad", &m_quad1Grad);
  m_tree->SetBranchAddress("quad2Grad", &m_quad2Grad);
  m_tree->SetBranchAddress("quad3Grad", &m_quad3Grad);
  m_tree->SetBranchAddress("xCorrectorStrength", &m_xCorrectorStrength);
  m_tree->SetBranchAddress("dipoleStrength", &m_dipoleStrength);

  // Measurement hits in the surface frame
  m_tree->SetBranchAddress("trackHitsLocal", &m_trackHitsLocal);

  // Measurement hits in the global frame
  m_tree->SetBranchAddress("trackHitsGlobal", &m_trackHitsGlobal);

  // Direction measurements in the track frame
  m_tree->SetBranchAddress("onSurfaceTrackDirection",
                           &m_onSurfaceTrackDirection);

  // Covariances of the track hits
  m_tree->SetBranchAddress("trackHitCovs", &m_trackHitCovs);

  // Covariances of the track agnles
  m_tree->SetBranchAddress("trackAngleCovs", &m_trackAngleCovs);

  // Geometry ids of the track hits
  m_tree->SetBranchAddress("geometryIds", &m_geometryIds);

  // KF predicted track hits in the surface frame
  m_tree->SetBranchAddress("predictedTrackHitsLocal",
                           &m_predictedTrackHitsLocal);
  m_tree->SetBranchAddress("filteredTrackHitsLocal", &m_filteredTrackHitsLocal);
  m_tree->SetBranchAddress("smoothedTrackHitsLocal", &m_smoothedTrackHitsLocal);

  // KF predicted track hits in the global frame
  m_tree->SetBranchAddress("predictedTrackHitsGlobal",
                           &m_predictedTrackHitsGlobal);
  m_tree->SetBranchAddress("filteredTrackHitsGlobal",
                           &m_filteredTrackHitsGlobal);
  m_tree->SetBranchAddress("smoothedTrackHitsGlobal",
                           &m_smoothedTrackHitsGlobal);

  // KF predicted on surface momenta in the track frame
  m_tree->SetBranchAddress("predictedOnSurfaceMomentum",
                           &m_predictedOnSurfaceMomentum);
  m_tree->SetBranchAddress("filteredOnSurfaceMomentum",
                           &m_filteredOnSurfaceMomentum);
  m_tree->SetBranchAddress("smoothedOnSurfaceMomentum",
                           &m_smoothedOnSurfaceMomentum);

  // KF residuals with respect to the measurements
  m_tree->SetBranchAddress("predictedHitResiduals", &m_predictedHitResiduals);
  m_tree->SetBranchAddress("filteredHitResiduals", &m_filteredHitResiduals);
  m_tree->SetBranchAddress("smoothedHitResiduals", &m_smoothedHitResiduals);

  m_tree->SetBranchAddress("predictedAngleResiduals",
                           &m_predictedAngleResiduals);
  m_tree->SetBranchAddress("filteredAngleResiduals", &m_filteredAngleResiduals);
  m_tree->SetBranchAddress("smoothedAngleResiduals", &m_smoothedAngleResiduals);

  // KF pulls with respect to the measurements
  m_tree->SetBranchAddress("predictedHitPulls", &m_predictedHitPulls);
  m_tree->SetBranchAddress("filteredHitPulls", &m_filteredHitPulls);
  m_tree->SetBranchAddress("smoothedHitPulls", &m_smoothedHitPulls);

  m_tree->SetBranchAddress("predictedAnglePulls", &m_predictedAnglePulls);
  m_tree->SetBranchAddress("filteredAnglePulls", &m_filteredAnglePulls);
  m_tree->SetBranchAddress("smoothedAnglePulls", &m_smoothedAnglePulls);

  // Guessed bound track parameters
  m_tree->SetBranchAddress("boundTrackParametersGuess",
                           &m_boundTrackParametersGuess);
  m_tree->SetBranchAddress("boundTrackCovGuess", &m_boundTrackCovGuess);

  // KF predicted bound track parameters
  m_tree->SetBranchAddress("boundTrackParametersEst",
                           &m_boundTrackParametersEst);
  m_tree->SetBranchAddress("boundTrackCovEst", &m_boundTrackCovEst);

  // Initial guess of the momentum at the IP
  m_tree->SetBranchAddress("originMomentumGuess", &m_originMomentumGuess);

  // Initial guess of the vertex at the IP
  m_tree->SetBranchAddress("vertexGuess", &m_vertexGuess);

  // KF predicted momentum at the IP
  m_tree->SetBranchAddress("originMomentumEst", &m_originMomentumEst);

  // KF predicted vertex at the IP
  m_tree->SetBranchAddress("vertexEst", &m_vertexEst);

  // Chi2 of the track with respect ot the measurement
  m_tree->SetBranchAddress("chi2Predicted", &m_chi2Predicted);
  m_tree->SetBranchAddress("chi2Filtered", &m_chi2Filtered);
  m_tree->SetBranchAddress("chi2Smoothed", &m_chi2Smoothed);

  // Number of degrees of freedom of the track
  m_tree->SetBranchAddress("ndf", &m_ndf);

  // Track ID
  m_tree->SetBranchAddress("trackId", &m_trackId);

  // PDG ID
  m_tree->SetBranchAddress("pdgId", &m_pdgId);

  // Charge
  m_tree->SetBranchAddress("charge", &m_charge);

  //------------------------------------------------------------------

  // Initialize constraint surfaces geometry ids
  const auto& goInst = *GeometryOptions::instance();
  m_constraintSurfacesGeoIds.insert(goInst.bpm0Parameters.geoId);
  m_constraintSurfacesGeoIds.insert(goInst.bpm1Parameters.geoId);
  m_constraintSurfacesGeoIds.insert(goInst.bpm2Parameters.geoId);
  m_constraintSurfacesGeoIds.insert(goInst.bpm3Parameters.geoId);
  m_constraintSurfacesGeoIds.insert(goInst.ipSurfaceParameters.geoId);

  // Initialize the data handles
  m_outputSourceLinks.initialize(m_cfg.outputSourceLinks);

  m_outputSeedsGuess.initialize(m_cfg.outputSeedsGuess);
  m_outputTrackParametersGuess.initialize(m_cfg.outputTrackParametersGuess);

  m_outputSeedsEst.initialize(m_cfg.outputSeedsEst);
  m_outputTrackParametersEst.initialize(m_cfg.outputTrackParametersEst);

  m_outputMagneticFieldParameters.initialize(
      m_cfg.outputMagneticFieldParameters);

  m_outputEventMetaData.initialize(m_cfg.outputEventMetaData);
}

std::pair<std::size_t, std::size_t> E320::E320RootTrackReader::availableEvents()
    const {
  return {std::get<0>(m_eventMap.front()), std::get<0>(m_eventMap.back()) + 1};
}

ProcessCode E320::E320RootTrackReader::read(const AlgorithmContext& ctx) {
  auto it = std::ranges::find_if(m_eventMap, [&](const auto& a) {
    return std::get<0>(a) == ctx.eventNumber;
  });

  if (it == m_eventMap.end()) {
    // Explicitly warn if it happens for the first or last event as that might
    // indicate a human error
    if ((ctx.eventNumber == availableEvents().first) &&
        (ctx.eventNumber == availableEvents().second - 1)) {
      ACTS_WARNING("Reading empty event: " << ctx.eventNumber);
    } else {
      ACTS_DEBUG("Reading empty event: " << ctx.eventNumber);
    }

    m_outputSourceLinks(ctx, {});

    m_outputSeedsGuess(ctx, {});
    m_outputTrackParametersGuess(ctx, {});

    m_outputSeedsEst(ctx, {});
    m_outputTrackParametersEst(ctx, {});

    // Return success flag
    return ProcessCode::SUCCESS;
  }

  // lock the mutex
  std::lock_guard<std::mutex> lock(m_read_mutex);

  ACTS_DEBUG("Reading event: " << std::get<0>(*it)
                               << " stored in entries: " << std::get<1>(*it)
                               << " - " << std::get<2>(*it));

  // Create the measurements
  std::vector<Acts::SourceLink> sourceLinks{};

  IndexSeeds seedsGuess{};
  std::vector<Acts::CurvilinearTrackParameters> trackParametersGuess{};

  IndexSeeds seedsEst{};
  std::vector<Acts::CurvilinearTrackParameters> trackParametersEst{};

  std::vector<std::shared_ptr<MagneticFieldStore>> magFieldStores;

  E320RootDataReader::EventMetaData eventMetaData{};

  const auto& goInst = *E320::GeometryOptions::instance();
  std::size_t longIdx = goInst.longIdx;
  std::size_t shortIdx = goInst.shortIdx;

  std::size_t quad1Id = goInst.quad1Id;
  std::size_t quad2Id = goInst.quad2Id;
  std::size_t quad3Id = goInst.quad3Id;
  std::size_t xCorrectorId = goInst.xCorrectorId;
  std::size_t dipoleId = goInst.dipoleId;

  std::size_t eventId = std::get<0>(*it);
  const Constraints& constraints = m_cfg.constraints;
  for (auto entry = std::get<1>(*it); entry < std::get<2>(*it); entry++) {
    m_tree->GetEntry(entry);

    if (constraints.requireEpicsParity &&
        m_epicsParity != constraints.requiredEpicsParity) {
      continue;
    }
    if (m_chi2Predicted < constraints.minPredictedChi2 ||
        m_chi2Predicted > constraints.maxPredictedChi2) {
      continue;
    }
    if (m_chi2Filtered < constraints.minFilteredChi2 ||
        m_chi2Filtered > constraints.maxFilteredChi2) {
      continue;
    }
    if (m_chi2Smoothed < constraints.minSmoothedChi2 ||
        m_chi2Smoothed > constraints.maxSmoothedChi2) {
      continue;
    }
    if (m_vertexEst->Y() < constraints.minVertexEstLong ||
        m_vertexEst->Y() > constraints.maxVertexEstLong) {
      continue;
    }
    if (m_vertexEst->Z() < constraints.minVertexEstShort ||
        m_vertexEst->Z() > constraints.maxVertexEstShort) {
      continue;
    }
    if (m_originMomentumEst->P() < constraints.minAbsMomentumEst ||
        m_originMomentumEst->P() > constraints.maxAbsMomentumEst) {
      continue;
    }

    bool passed = true;
    for (std::size_t j = 0; j < m_smoothedHitResiduals->size(); j++) {
      std::size_t geoId = m_geometryIds->at(j);
      const auto& smoothedRes = m_smoothedHitResiduals->at(j);
      const auto& [minResX, maxResX, minResY, maxResY] =
          constraints.smoothedResidualsRanges.at(geoId);
      if (smoothedRes.X() < minResX || smoothedRes.X() > maxResX ||
          smoothedRes.Y() < minResY || smoothedRes.Y() > maxResY) {
        passed = false;
        break;
      }
    }
    if (!passed) {
      continue;
    }

    eventMetaData =
        E320RootDataReader::EventMetaData{.eudaqTrgN = m_eudaqTrgN,
                                          .eudaqDAQNumber = m_eudaqDAQNumber,
                                          .eudaqRunStartTs = m_eudaqRunStartTs,
                                          .eudaqRunEndTs = m_eudaqRunEndTs,
                                          .epicsParity = m_epicsParity,
                                          .epicsPulseId = m_epicsPulseId,
                                          .epicsDAQNumber = m_epicsDAQNumber};

    // Magnetic fields
    auto mStore = std::make_shared<MagneticFieldStore>();
    mStore->store.reserve(5);

    // Quad gradients T/m
    mStore->store.insert(
        {goInst.quad1Id,
         Acts::MagneticFieldProvider::Cache(
             std::in_place_type<IdealQuadrupoleMagField::Cache>,
             m_quad1Grad * Acts::UnitConstants::T / Acts::UnitConstants::m)});
    mStore->store.insert(
        {goInst.quad2Id,
         Acts::MagneticFieldProvider::Cache(
             std::in_place_type<IdealQuadrupoleMagField::Cache>,
             m_quad2Grad * Acts::UnitConstants::T / Acts::UnitConstants::m)});
    mStore->store.insert(
        {goInst.quad3Id,
         Acts::MagneticFieldProvider::Cache(
             std::in_place_type<IdealQuadrupoleMagField::Cache>,
             m_quad3Grad * Acts::UnitConstants::T / Acts::UnitConstants::m)});

    // XCOR strength T
    mStore->store.insert(
        {goInst.xCorrectorId,
         Acts::MagneticFieldProvider::Cache(
             std::in_place_type<ConstantMagField::Cache>,
             m_xCorrectorStrength * Acts::UnitConstants::T, longIdx)});

    // Dipole strength T
    mStore->store.insert(
        {goInst.dipoleId,
         Acts::MagneticFieldProvider::Cache(
             std::in_place_type<ConstantMagField::Cache>,
             m_dipoleStrength * Acts::UnitConstants::T, shortIdx)});

    // Store the fields
    magFieldStores.push_back(mStore);

    // Track origin covariances
    Acts::BoundMatrix originCovGuess;
    Acts::BoundMatrix originCovEst;
    for (std::size_t i = 0; i < Acts::eBoundSize; i++) {
      for (std::size_t j = 0; j < Acts::eBoundSize; j++) {
        originCovGuess(i, j) = (*m_boundTrackCovGuess)(i, j);
        originCovEst(i, j) = (*m_boundTrackCovEst)(i, j);
      }
    }

    Acts::ParticleHypothesis hypothesis =
        Acts::ParticleHypothesis(Acts::PdgParticle(m_pdgId));

    std::vector<std::size_t> trackSourceLinkIndices{};
    for (std::size_t i = 0; i < m_trackHitsGlobal->size(); i++) {
      if (m_constraintSurfacesGeoIds.contains(m_geometryIds->at(i))) {
        Acts::ActsVector<ExtendedLocalSize> measLocal(
            m_trackHitsLocal->at(i).X(), m_trackHitsLocal->at(i).Y(),
            m_onSurfaceTrackDirection->at(i).Phi(),
            m_onSurfaceTrackDirection->at(i).Theta());

        Acts::ActsVector<ExtendedGlobalSize> measGlobal(
            m_trackHitsGlobal->at(i).X(), m_trackHitsGlobal->at(i).Y(),
            m_trackHitsGlobal->at(i).Z(), m_onSurfaceTrackDirection->at(i).X(),
            m_onSurfaceTrackDirection->at(i).Y(),
            m_onSurfaceTrackDirection->at(i).Z());

        Acts::SquareMatrix2 trackHitCov;
        trackHitCov << m_trackHitCovs->at(i)(0, 0), m_trackHitCovs->at(i)(0, 1),
            m_trackHitCovs->at(i)(1, 0), m_trackHitCovs->at(i)(1, 1);
        Acts::SquareMatrix2 trackAngleCov;
        trackAngleCov << m_trackAngleCovs->at(i)(0, 0),
            m_trackAngleCovs->at(i)(0, 1), m_trackAngleCovs->at(i)(1, 0),
            m_trackAngleCovs->at(i)(1, 1);

        Acts::ActsSquareMatrix<ExtendedLocalSize> measCov =
            Acts::ActsSquareMatrix<ExtendedLocalSize>::Zero();
        measCov.topLeftCorner(2, 2) = trackHitCov;
        measCov.bottomRightCorner(2, 2) = trackAngleCov;

        Acts::GeometryIdentifier geoId;
        geoId.setSensitive(m_geometryIds->at(i));
        ExtendedSourceLink esl(measLocal, measGlobal, measCov, geoId, eventId,
                               sourceLinks.size(), m_cfg.backwards);
        trackSourceLinkIndices.push_back(sourceLinks.size());
        Acts::SourceLink sl(esl);
        sourceLinks.push_back(sl);
      } else {
        Acts::Vector2 trackHitLocal(m_trackHitsLocal->at(i).X(),
                                    m_trackHitsLocal->at(i).Y());
        Acts::Vector3 trackHitGlobal(m_trackHitsGlobal->at(i).X(),
                                     m_trackHitsGlobal->at(i).Y(),
                                     m_trackHitsGlobal->at(i).Z());
        Acts::SquareMatrix2 trackHitCov;
        trackHitCov << m_trackHitCovs->at(i)(0, 0), m_trackHitCovs->at(i)(0, 1),
            m_trackHitCovs->at(i)(1, 0), m_trackHitCovs->at(i)(1, 1);

        // trackHitCov = Acts::Vector2(5_um, 5_um).cwiseAbs2().asDiagonal();

        Acts::GeometryIdentifier geoId;
        geoId.setSensitive(m_geometryIds->at(i));
        SimpleSourceLink ssl(trackHitLocal, trackHitGlobal, trackHitCov, geoId,
                             eventId, sourceLinks.size());
        trackSourceLinkIndices.push_back(sourceLinks.size());
        Acts::SourceLink sl(ssl);
        sourceLinks.push_back(sl);
      }
    }

    // Track parameters initial guess
    Acts::Vector4 vertexGuess(m_vertexGuess->X(), m_vertexGuess->Y(),
                              m_vertexGuess->Z(), 0);
    Acts::Vector3 ipDirectionGuess(m_originMomentumGuess->X(),
                                   m_originMomentumGuess->Y(),
                                   m_originMomentumGuess->Z());
    ipDirectionGuess.normalize();

    seedsGuess.emplace_back(trackSourceLinkIndices, trackParametersGuess.size(),
                            static_cast<int>(seedsGuess.size()));
    trackParametersGuess.emplace_back(vertexGuess, ipDirectionGuess,
                                      m_charge / m_originMomentumGuess->P(),
                                      originCovGuess, hypothesis);

    // Estimated track parameters
    Acts::Vector4 vertexEst(m_vertexEst->X(), m_vertexEst->Y(),
                            m_vertexEst->Z(), 0);
    Acts::Vector3 ipDirectionEst(m_originMomentumEst->X(),
                                 m_originMomentumEst->Y(),
                                 m_originMomentumEst->Z());
    ipDirectionEst.normalize();

    seedsEst.emplace_back(trackSourceLinkIndices, trackParametersEst.size(),
                          static_cast<int>(seedsEst.size()));
    trackParametersEst.emplace_back(vertexEst, ipDirectionEst,
                                    m_charge / m_originMomentumEst->P(),
                                    originCovEst, hypothesis);
  }

  ACTS_DEBUG("Read " << sourceLinks.size() << " source links");

  ACTS_DEBUG("Read " << seedsGuess.size() << " guess seeds");
  ACTS_DEBUG("Read " << trackParametersGuess.size()
                     << " guess track parameters");

  ACTS_DEBUG("Read " << seedsEst.size() << " estimated seeds");
  ACTS_DEBUG("Read " << trackParametersEst.size()
                     << " estimated track parameters");

  m_outputSourceLinks(ctx, std::move(sourceLinks));

  m_outputSeedsGuess(ctx, std::move(seedsGuess));
  m_outputTrackParametersGuess(ctx, std::move(trackParametersGuess));

  m_outputSeedsEst(ctx, std::move(seedsEst));
  m_outputTrackParametersEst(ctx, std::move(trackParametersEst));

  m_outputMagneticFieldParameters(ctx, std::move(magFieldStores));

  m_outputEventMetaData(ctx, std::move(eventMetaData));

  // Return success flag
  return ProcessCode::SUCCESS;
}
