#include "TrackingPipeline/TrackFinding/ApollonSeedingAlgorithm.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include <Acts/Definitions/Algebra.hpp>

#include <cstddef>
#include <functional>
#include <limits>
#include <span>

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

  m_det1FirstLayerPoint[goInst.primaryIdx] = goInst.ipTc1Distance;
  m_det1FirstLayerPoint[goInst.longIdx] = goInst.tc1CenterLong;
  m_det1FirstLayerPoint[goInst.shortIdx] = goInst.tc1CenterShort;
  m_det1FirstLayerNormal = goInst.primaryDir;

  m_det2FirstLayerPoint[goInst.primaryIdx] = goInst.ipTc2Distance;
  m_det2FirstLayerPoint[goInst.longIdx] = goInst.tc2CenterLong;
  m_det2FirstLayerPoint[goInst.shortIdx] = goInst.tc2CenterShort;
  m_det2FirstLayerNormal = goInst.primaryDir;

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

  const auto& goInst = *go::instance();
  HoughTransformSeeder::Options htSeederOpt;
  htSeederOpt.boundBoxCenterX = goInst.tc1CenterPrimary;
  htSeederOpt.boundBoxCenterY = goInst.tc1CenterLong;
  htSeederOpt.boundBoxCenterZ = goInst.tc1CenterShort;

  std::vector<std::reference_wrapper<const Acts::SourceLink>>
      sourceLinkRefsDet1;
  sourceLinkRefsDet1.reserve(inputSourceLinks.size());
  // std::cout << "\n\n\n-------------------------------------\n";
  // std::cout << "DET 1 SOURCE LINKS\n";
  for (const auto& sl : inputSourceLinks) {
    if (sl.get<SimpleSourceLink>().geometryId().sensitive() < 20) {
      sourceLinkRefsDet1.push_back(std::cref(sl));

      const auto& ssl = sl.get<SimpleSourceLink>();
      // std::cout << "GEOID " << ssl.geometryId() << "\n";
      // std::cout << "LOCAL HIT " << ssl.parametersLoc().transpose() << "\n";
      // std::cout << "GLOBAL HIT " << ssl.parametersGlob().transpose() << "\n";
    }
  }
  sourceLinkRefsDet1.shrink_to_fit();
  auto htSeedsDet1 = m_cfg.htSeeder->findSeeds(sourceLinkRefsDet1, htSeederOpt);
  // std::cout << "\n\n\n";

  ACTS_DEBUG("Found " << htSeedsDet1.size()
                      << " seed connection candidates in the first detector");

  // Seeds outSeeds;
  // outSeeds.reserve(htSeedsDet1.size());
  // for (int i = 0; i < htSeedsDet1.size(); i++) {
  //   const auto& seed = htSeedsDet1.at(i);
  //   outSeeds.emplace_back(
  //       seed.sourceLinks,
  //       Acts::CurvilinearTrackParameters(
  //           Acts::Vector4(), Acts::Vector3::UnitX(), -1_e / 1_GeV, m_ipCov,
  //           Acts::ParticleHypothesis::electron()),
  //       i);
  // }

  std::vector<std::reference_wrapper<const Acts::SourceLink>>
      sourceLinkRefsDet2;
  sourceLinkRefsDet2.reserve(inputSourceLinks.size());
  // std::cout << "\n\n\n-------------------------------------\n";
  // std::cout << "DET 2 SOURCE LINKS\n";
  for (const auto& sl : inputSourceLinks) {
    if (sl.get<SimpleSourceLink>().geometryId().sensitive() >= 20) {
      sourceLinkRefsDet2.push_back(std::cref(sl));

      const auto& ssl = sl.get<SimpleSourceLink>();
      // std::cout << "GEOID " << ssl.geometryId() << "\n";
      // std::cout << "LOCAL HIT " << ssl.parametersLoc().transpose() << "\n";
      // std::cout << "GLOBAL HIT " << ssl.parametersGlob().transpose() << "\n";
    }
  }
  sourceLinkRefsDet2.shrink_to_fit();
  auto htSeedsDet2 = m_cfg.htSeeder->findSeeds(sourceLinkRefsDet2, htSeederOpt);
  // std::cout << "\n\n\n";

  ACTS_DEBUG("Found " << htSeedsDet2.size()
                      << " seed connection candidates in the second detector");

  Seeds outSeeds = scanEnergy(htSeedsDet1, htSeedsDet2);
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
  for (const auto& [point1, dir1, sl1] : det1Seeds) {
    // std::cout << "POINT1 " << point1.transpose() << "\n";
    // std::cout << "DIR1 " << dir1.transpose() << "\n";

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

    // std::cout << "DIPOLE ENTRANCE POINT " << dipoleEntrancePoint.transpose()
    //           << "\n";

    double dirX = dir1.x();
    double dirY = dir1.y();
    double dirZ = dir1.z();

    // std::cout << "DIPOLE ENTRANCE DIR " << dir1.transpose() << "\n";

    double dirTrans = std::hypot(dirX, dirY);
    double asinDir = std::asin(dirY / dirTrans);

    for (double P = m_cfg.minScanEnergy; P <= m_cfg.maxScanEnergy;
         P += m_cfg.energyScanStep) {
      std::size_t idxMin = 0;
      double dTotMin = std::numeric_limits<double>().max();
      double minP = 0;

      // std::cout << "--------------------------------\n";
      // std::cout << "CHECKING MOMENTUM " << P << "\n";
      double x0 = x - P * dirY / B;
      double y0 = y + P * dirX / B;

      // std::cout << "x0 " << x0 << "\n";
      // std::cout << "y0 " << y0 << "\n";

      double rho2 = (x - x0) * (x - x0) + (y - y0) * (y - y0);

      // std::cout << "rho2 " << rho2 << "\n";

      double deltaX = xPrime - x0;
      // std::cout << "deltaX " << deltaX << "\n";

      double deltaY = std::sqrt(rho2 - deltaX * deltaX);
      // std::cout << "deltaY " << deltaY << "\n";

      double yPrime = y0 + sign * deltaY;
      // std::cout << "yPrime " << yPrime << "\n";

      double arg = -deltaX / (yPrime - y0);
      // std::cout << "arg " << arg << "\n";
      double denom = std::sqrt(1 + arg * arg);
      // std::cout << "denom " << denom << "\n";

      dipoleExitDir[primaryIdx] = 1.0 / denom;
      dipoleExitDir[longIdx] = arg / denom;
      dipoleExitDir[shortIdx] = dirZ;

      double zPrime =
          z + P * dirZ / B *
                  (std::asin(dipoleExitDir[longIdx] / dirTrans) - asinDir);

      dipoleExitPoint[longIdx] = yPrime;
      dipoleExitPoint[shortIdx] = zPrime;
      // std::cout << "DIPOLE EXIT POINT " << dipoleExitPoint.transpose() << "\n";
      // std::cout << "DIPOLE EXIT DIR " << dipoleExitDir.transpose() << "\n";
      // std::cout << "DET 2 FL POINT " << m_det2FirstLayerPoint.transpose()
      //           << "\n";

      double dDet2 = (m_det2FirstLayerPoint - dipoleExitPoint)
                         .dot(m_det2FirstLayerNormal) /
                     dipoleExitDir.dot(m_det2FirstLayerNormal);

      // std::cout << "DET 2 DISTANCE " << dDet2 << "\n";
      Acts::Vector3 det2FirstLayerPoint =
          dipoleExitPoint + dipoleExitDir * dDet2;

      // std::cout << "DET2 FIRST LAYER POINT " <<
      // det2FirstLayerPoint.transpose()
      //           << "\n";
      for (std::size_t idx = 0; idx < activeIdxs.size(); idx++) {
        if (!activeIdxs.at(idx)) {
          continue;
        }
        const auto& [point2, dir2, sl2] = det2Seeds.at(idx);
        // std::cout << "DET 2 LINE POINT " << point2.transpose() << "\n";
        // std::cout << "DET 2 LINE DIR " << dir2.transpose() << "\n";
        Acts::Vector3 dirCrossProd = dipoleExitDir.cross(dir2);
        double dAngle = std::asin(dirCrossProd.norm());
        double dOrtho =
            std::abs((point2 - det2FirstLayerPoint).dot(dirCrossProd));
        // std::cout << "DANGLE " << dAngle << "\n";
        // std::cout << "DORTHO " << dOrtho << "\n";
        double dTot = dAngle / m_dAngleMax + dOrtho / m_dOrthoMax;
        // std::cout << "DTOT " << dTot << "\n";
        if (dTot < dTotMin) {
          dTotMin = dTot;
          idxMin = idx;
          minP = P;
        }
      }
      if (dTotMin < m_cfg.maxConnectionDistance) {
        // activeIdxs.at(idxMin) = false;

        std::size_t totSize =
            sl1.size() + det2Seeds.at(idxMin).sourceLinks.size();
        // std::cout << "TOT SIZE " << totSize << "\n";
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
        // std::cout << "TOT LAYERS " << geoIdSet.size() << "\n";
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
            Acts::Vector4(det1FirstLayerPoint.x() - 0.1_mm, det1FirstLayerPoint.y(),
                          det1FirstLayerPoint.z(), 0),
            dir1, -1_e / minP, m_ipCov, Acts::ParticleHypothesis::electron());

        outSeeds.emplace_back(std::move(seedSourceLinks), std::move(startPars),
                              outSeeds.size());
      }
    }
  }
  return outSeeds;
}
