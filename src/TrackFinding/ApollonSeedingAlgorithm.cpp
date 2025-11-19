#include "TrackingPipeline/TrackFinding/ApollonSeedingAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <cstddef>
#include <functional>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include <omp.h>

#include "TrackingPipeline/EventData/DataContainers.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"
#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

using go = ApollonGeometry::GeometryOptions;

using namespace Acts::UnitLiterals;

ApollonSeedingAlgorithm::ApollonSeedingAlgorithm(const Config& config,
                                                 Acts::Logging::Level level)
    : IAlgorithm("ApollonSeedingAlgorithm", level), m_cfg(config) {
  m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
  m_outputSeeds.initialize(m_cfg.outputSeeds);

  Acts::BoundVector ipStdDev;
  ipStdDev[Acts::eBoundLoc0] = 100_um;
  ipStdDev[Acts::eBoundLoc1] = 100_um;
  ipStdDev[Acts::eBoundTime] = 25_ns;
  ipStdDev[Acts::eBoundPhi] = 2_degree;
  ipStdDev[Acts::eBoundTheta] = 2_degree;
  ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
  m_ipCov = ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

  const auto& goInst = *go::instance();

  m_det1FirstLayerPoint[goInst.primaryIdx] = goInst.ipTc1Distance - 0.3_mm;
  m_det1FirstLayerPoint[goInst.longIdx] = goInst.tc1CenterLong;
  m_det1FirstLayerPoint[goInst.shortIdx] = goInst.tc1CenterShort;
  m_det1FirstLayerNormal = goInst.primaryDir;

  m_det2FirstLayerPoint[goInst.primaryIdx] = goInst.ipTc2Distance;
  m_det2FirstLayerPoint[goInst.longIdx] = goInst.tc2CenterLong;
  m_det2FirstLayerPoint[goInst.shortIdx] = goInst.tc2CenterShort;
  m_det2FirstLayerNormal = goInst.primaryDir;

  m_det2LastLayerPoint[goInst.primaryIdx] =
      goInst.ipTc2Distance + 4 * goInst.interChipDistance;
  m_det2LastLayerPoint[goInst.longIdx] = goInst.tc2CenterLong;
  m_det2LastLayerPoint[goInst.shortIdx] = goInst.tc2CenterShort;
  m_det2LastLayerNormal = goInst.primaryDir;

  m_dipoleEntrancePoint[goInst.primaryIdx] =
      goInst.dipoleCenterPrimary - goInst.dipoleHalfPrimary;
  m_dipoleEntrancePoint[goInst.longIdx] = goInst.dipoleCenterLong;
  m_dipoleEntrancePoint[goInst.shortIdx] = goInst.dipoleCenterShort;
  m_dipoleEntranceNormal = goInst.primaryDir;

  m_dipoleExitPoint[goInst.primaryIdx] =
      goInst.dipoleCenterPrimary + goInst.dipoleHalfPrimary;
  m_dipoleExitPoint[goInst.longIdx] = goInst.dipoleCenterLong;
  m_dipoleExitPoint[goInst.shortIdx] = goInst.dipoleCenterShort;
  m_dipoleExitNormal = goInst.primaryDir;

  m_dAngleMax = M_PI;
  m_dOrthoMax = 2 * goInst.tcHalfShort;
}

ProcessCode ApollonSeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  const auto& inputSourceLinks = m_inputSourceLinks(ctx);

  ACTS_DEBUG("Received " << inputSourceLinks.size() << " source links");

  if (inputSourceLinks.empty()) {
    ACTS_DEBUG("Input is empty. Skipping");
    m_outputSeeds(ctx, Seeds());
    return ProcessCode::SUCCESS;
  }

  Seeds outSeeds;
  std::vector<HoughTransformSeeder::HTSeed> htSeedsDet1;
  std::vector<HoughTransformSeeder::HTSeed> htSeedsDet2;
  const auto& goInst = *go::instance();
  if (m_cfg.scope == SeedingScope::detector1 ||
      m_cfg.scope == SeedingScope::fullDetector) {
    HoughTransformSeeder::Options htSeederOpt1;
    htSeederOpt1.boundBoxCenterX = goInst.tc1CenterPrimary;
    htSeederOpt1.boundBoxCenterY = goInst.tc1CenterLong;
    htSeederOpt1.boundBoxCenterZ = goInst.tc1CenterShort;

    htSeederOpt1.firstLayerId = goInst.tc1Parameters.front().geoId;
    htSeederOpt1.lastLayerId = goInst.tc1Parameters.back().geoId;
    htSeederOpt1.nLayers = goInst.tc1Parameters.size();

    htSeederOpt1.minCount = m_cfg.minXCountDet1;
    htSeederOpt1.maxChi2 = m_cfg.maxLineChi2Det1;

    std::vector<std::reference_wrapper<const Acts::SourceLink>>
        sourceLinkRefsDet1;
    sourceLinkRefsDet1.reserve(inputSourceLinks.size());
    for (const auto& sl : inputSourceLinks) {
      if (sl.get<SimpleSourceLink>().geometryId().sensitive() <=
          htSeederOpt1.lastLayerId) {
        sourceLinkRefsDet1.push_back(std::cref(sl));
      }
    }
    sourceLinkRefsDet1.shrink_to_fit();
    htSeedsDet1 = m_cfg.htSeeder->findSeeds(sourceLinkRefsDet1, htSeederOpt1);

    ACTS_DEBUG("Found " << htSeedsDet1.size()
                        << " seed connection candidates in the first detector");
  }
  if (m_cfg.scope == SeedingScope::detector2 ||
      m_cfg.scope == SeedingScope::fullDetector) {
    HoughTransformSeeder::Options htSeederOpt2;
    htSeederOpt2.boundBoxCenterX = goInst.tc2CenterPrimary;
    htSeederOpt2.boundBoxCenterY = goInst.tc2CenterLong;
    htSeederOpt2.boundBoxCenterZ = goInst.tc2CenterShort;

    htSeederOpt2.firstLayerId = goInst.tc2Parameters.front().geoId;
    htSeederOpt2.lastLayerId = goInst.tc2Parameters.back().geoId;
    htSeederOpt2.nLayers = goInst.tc2Parameters.size();

    htSeederOpt2.minCount = m_cfg.minXCountDet2;
    htSeederOpt2.maxChi2 = m_cfg.maxLineChi2Det2;

    std::vector<std::reference_wrapper<const Acts::SourceLink>>
        sourceLinkRefsDet2;
    sourceLinkRefsDet2.reserve(inputSourceLinks.size());
    for (const auto& sl : inputSourceLinks) {
      if (sl.get<SimpleSourceLink>().geometryId().sensitive() >=
          htSeederOpt2.firstLayerId) {
        sourceLinkRefsDet2.push_back(std::cref(sl));
      }
    }
    sourceLinkRefsDet2.shrink_to_fit();
    htSeedsDet2 = m_cfg.htSeeder->findSeeds(sourceLinkRefsDet2, htSeederOpt2);

    ACTS_DEBUG(
        "Found " << htSeedsDet2.size()
                 << " seed connection candidates in the second detector");
  }
  if (m_cfg.scope == SeedingScope::detector1 ||
      m_cfg.scope == SeedingScope::detector2) {
    std::vector<HoughTransformSeeder::HTSeed> htSeeds =
        (m_cfg.scope == SeedingScope::detector1) ? std::move(htSeedsDet1)
                                                 : std::move(htSeedsDet2);
    Acts::Vector3 detFirstLayerPoint = (m_cfg.scope == SeedingScope::detector1)
                                           ? m_det1FirstLayerPoint
                                           : m_det2FirstLayerPoint;
    Acts::Vector3 detFirstLayerNormal = (m_cfg.scope == SeedingScope::detector1)
                                            ? m_det1FirstLayerNormal
                                            : m_det2FirstLayerNormal;

    outSeeds.reserve(htSeeds.size());
    #pragma omp parallel for num_threads(32)
    for (std::size_t i = 0; i < htSeeds.size(); i++) {
      const auto& [point, dir, sl] = htSeeds.at(i);
      double dVertex =
          (detFirstLayerPoint - Acts::Vector3(0.3_mm, 0, 0) - point)
              .dot(detFirstLayerNormal) /
          dir.dot(detFirstLayerNormal);
      Acts::Vector3 vertex3 = point + dir * dVertex;
      Acts::Vector4 vertex(vertex3.x(), vertex3.y(), vertex3.z(), 0);

      outSeeds.emplace_back(std::move(sl),
                            Acts::CurvilinearTrackParameters(
                                vertex, dir, -1_e / 1_GeV, m_ipCov,
                                Acts::ParticleHypothesis::electron()),
                            i);
    }
  } else {
    outSeeds = scanEnergy(htSeedsDet1, htSeedsDet2);
  }
  ACTS_DEBUG("Found " << outSeeds.size() << " seeds");
  ACTS_DEBUG("Sending " << outSeeds.size() << " seeds");
  m_outputSeeds(ctx, std::move(outSeeds));

  return ProcessCode::SUCCESS;
}

Seeds ApollonSeedingAlgorithm::scanEnergy(
    const std::vector<HoughTransformSeeder::HTSeed>& det1Seeds,
    const std::vector<HoughTransformSeeder::HTSeed>& det2Seeds) const {
  Seeds outSeeds;
  outSeeds.reserve(det1Seeds.size());

  std::vector<bool> activeIdxs;
  activeIdxs.reserve(det2Seeds.size());
  for (std::size_t i = 0; i < det2Seeds.size(); i++) {
    activeIdxs.push_back(true);
  }

  const auto& goInst = *go::instance();

  std::size_t primaryIdx = goInst.primaryIdx;
  std::size_t longIdx = goInst.longIdx;
  std::size_t shortIdx = goInst.shortIdx;

  double B = goInst.dipoleFieldShort;
  double sign = (B > 0) ? -1 : 1;

  Acts::Vector3 dipoleExitPoint(0, 0, 0);
  dipoleExitPoint[primaryIdx] =
      goInst.dipoleCenterPrimary + goInst.dipoleHalfPrimary;

  double xPrime = goInst.dipoleCenterPrimary + goInst.dipoleHalfPrimary;

  Acts::Vector3 dipoleExitDir(0, 0, 0);
#pragma omp parallel for num_threads(32)
  for (std::size_t i = 0; i < det1Seeds.size(); i++) {
    const auto& [point1, dir1, sl1] = det1Seeds.at(i);
    double dDet1FirstLayer =
        (m_det1FirstLayerPoint - point1).dot(m_det1FirstLayerNormal) /
        dir1.dot(m_det1FirstLayerNormal);
    Acts::Vector3 det1FirstLayerPoint = point1 + dir1 * dDet1FirstLayer;

    double dDipoleEntrance =
        (m_dipoleEntrancePoint - point1).dot(m_dipoleEntranceNormal) /
        dir1.dot(m_dipoleEntranceNormal);
    Acts::Vector3 dipoleEntrancePoint = point1 + dir1 * dDipoleEntrance;

    double x = dipoleEntrancePoint.x();
    double y = dipoleEntrancePoint.y();
    double z = dipoleEntrancePoint.z();

    double dirX = dir1.x();
    double dirY = dir1.y();
    double dirZ = dir1.z();

    double dirTrans = std::hypot(dirX, dirY);
    double asinDir = std::asin(dirY / dirTrans);

    std::size_t idxMin = 0;
    double dTotMin = std::numeric_limits<double>().max();
    double minP = 0;
    if (B == 0) {
      minP = m_cfg.minScanEnergy;

      double dDet2First =
          (m_det2FirstLayerPoint - point1).dot(m_det2FirstLayerNormal) /
          dir1.dot(m_det2FirstLayerNormal);
      Acts::Vector3 det2FirstLayerPoint1 = point1 + dir1 * dDet2First;

      double dDet2Last =
          (m_det2LastLayerPoint - point1).dot(m_det2LastLayerNormal) /
          dir1.dot(m_det2LastLayerNormal);
      Acts::Vector3 det2LastLayerPoint1 = point1 + dir1 * dDet2Last;

      for (std::size_t idx = 0; idx < activeIdxs.size(); idx++) {
        if (!activeIdxs.at(idx)) {
          continue;
        }
        const auto& [point2, dir2, sl2] = det2Seeds.at(idx);
        double dDet2First2 =
            (m_det2FirstLayerPoint - point2).dot(m_det2FirstLayerNormal) /
            dir2.dot(m_det2FirstLayerNormal);
        Acts::Vector3 det2FirstLayerPoint2 = point2 + dir2 * dDet2First2;

        double dDet2Last2 =
            (m_det2LastLayerPoint - point2).dot(m_det2LastLayerNormal) /
            dir2.dot(m_det2LastLayerNormal);
        Acts::Vector3 det2LastLayerPoint2 = point2 + dir2 * dDet2Last2;

        double dTot = (det2FirstLayerPoint1 - det2FirstLayerPoint2).norm() +
                      (det2LastLayerPoint1 - det2LastLayerPoint2).norm();
        if (dTot < dTotMin) {
          dTotMin = dTot;
          idxMin = idx;
        }
      }
    } else {
      for (double P = m_cfg.minScanEnergy; P <= m_cfg.maxScanEnergy;
           P += m_cfg.energyScanStep) {
        double x0 = x - P * dirY / B;
        double y0 = y + P * dirX / B;
        double rho2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);

        double deltaX = xPrime - x0;
        double deltaY = std::sqrt(rho2 - deltaX * deltaX);

        double yPrime = y0 + sign * deltaY;

        double arg = -deltaX / (yPrime - y0);
        double denom = std::sqrt(1 + arg * arg);
        dipoleExitDir[primaryIdx] = 1.0 / denom;
        dipoleExitDir[longIdx] = arg / denom;
        dipoleExitDir[shortIdx] = dirZ;
        dipoleExitDir /= dipoleExitDir.norm();

        double zPrime =
            z + P * dirZ / B *
                    (std::asin(dipoleExitDir[longIdx] / dirTrans) - asinDir);

        dipoleExitPoint[longIdx] = yPrime;
        dipoleExitPoint[shortIdx] = zPrime;

        double dDet2First = (m_det2FirstLayerPoint - dipoleExitPoint)
                                .dot(m_det2FirstLayerNormal) /
                            dipoleExitDir.dot(m_det2FirstLayerNormal);
        Acts::Vector3 det2FirstLayerPoint =
            dipoleExitPoint + dipoleExitDir * dDet2First;

        double dDet2Last = (m_det2LastLayerPoint - dipoleExitPoint)
                               .dot(m_det2LastLayerNormal) /
                           dipoleExitDir.dot(m_det2LastLayerNormal);
        Acts::Vector3 det2LastLayerPoint =
            dipoleExitPoint + dipoleExitDir * dDet2Last;

        for (std::size_t idx = 0; idx < activeIdxs.size(); idx++) {
          if (!activeIdxs.at(idx)) {
            continue;
          }
          const auto& [point2, dir2, sl2] = det2Seeds.at(idx);
          double dDet2First2 =
              (m_det2FirstLayerPoint - point2).dot(m_det2FirstLayerNormal) /
              dir2.dot(m_det2FirstLayerNormal);
          Acts::Vector3 det2FirstLayerPoint2 = point2 + dir2 * dDet2First2;

          double dDet2Last2 =
              (m_det2LastLayerPoint - point2).dot(m_det2LastLayerNormal) /
              dir2.dot(m_det2LastLayerNormal);
          Acts::Vector3 det2LastLayerPoint2 = point2 + dir2 * dDet2Last2;

          double dTot = (det2FirstLayerPoint - det2FirstLayerPoint2).norm() +
                        (det2LastLayerPoint - det2LastLayerPoint2).norm();
          if (dTot < dTotMin) {
            dTotMin = dTot;
            idxMin = idx;
            minP = P;
          }
        }
      }
    }
    if (dTotMin < m_cfg.maxConnectionDistance) {
      std::size_t totSize =
          sl1.size() + det2Seeds.at(idxMin).sourceLinks.size();
      if (totSize < m_cfg.minSeedSize || totSize > m_cfg.maxSeedSize) {
        continue;
      }
      std::set<int> geoIdSet;
      for (const auto& sl : sl1) {
        geoIdSet.insert(sl.get<SimpleSourceLink>().geometryId().sensitive());
      }
      for (const auto& sl : det2Seeds.at(idxMin).sourceLinks) {
        geoIdSet.insert(sl.get<SimpleSourceLink>().geometryId().sensitive());
      }
      if (geoIdSet.size() < m_cfg.minLayers ||
          geoIdSet.size() > m_cfg.maxLayers) {
        continue;
      }

      std::vector<Acts::SourceLink> seedSourceLinks;
      seedSourceLinks.reserve(totSize);
      seedSourceLinks.insert(seedSourceLinks.end(), sl1.begin(), sl1.end());
      seedSourceLinks.insert(seedSourceLinks.end(),
                             det2Seeds.at(idxMin).sourceLinks.begin(),
                             det2Seeds.at(idxMin).sourceLinks.end());

      Acts::CurvilinearTrackParameters startPars(
          Acts::Vector4(det1FirstLayerPoint.x(), det1FirstLayerPoint.y(),
                        det1FirstLayerPoint.z(), 0),
          dir1, -1_e / minP, m_ipCov, Acts::ParticleHypothesis::electron());

      outSeeds.emplace_back(std::move(seedSourceLinks), std::move(startPars),
                            outSeeds.size());
    }
  }
  return outSeeds;
}
