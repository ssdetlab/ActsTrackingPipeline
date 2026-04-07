#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <Eigen/src/Core/Matrix.h>
#include <TFile.h>
#include <TTree.h>
#include <omp.h>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/detail/BinningValueUtils.hpp"

std::vector<std::size_t> HoughTransformSeeder::findLineSourceLinks(
    const std::vector<Acts::SourceLink> &sourceLinks,
    const std::vector<std::size_t> &sourceLinksIndices,
    const Acts::Vector3 &pointBL, const Acts::Vector3 &dirBL,
    const Acts::Vector3 &pointTL, const Acts::Vector3 &dirTL,
    const Acts::Vector3 &pointBR, const Acts::Vector3 &dirBR,
    const Acts::Vector3 &pointTR, const Acts::Vector3 &dirTR,
    const Acts::Vector3 &shift) {
  std::size_t primaryIdx = m_cfg.primaryIdx;
  std::size_t longIdx = m_cfg.longIdx;
  std::size_t shortIdx = m_cfg.shortIdx;

  Acts::Vector3 normalPrimary = detail::indexToDirection(primaryIdx);

  std::unordered_map<Acts::GeometryIdentifier, std::array<Acts::Vector3, 4>>
      geoRefPointMap;
  geoRefPointMap.reserve(m_geoPosMap.size());
  for (const auto &[id, x] : m_geoPosMap) {
    Acts::Vector3 meas(x, 0, 0);
    double ddBL =
        (meas - pointBL).dot(normalPrimary) / dirBL.dot(normalPrimary);
    double ddTL =
        (meas - pointTL).dot(normalPrimary) / dirTL.dot(normalPrimary);

    double ddBR =
        (meas - pointBR).dot(normalPrimary) / dirBR.dot(normalPrimary);
    double ddTR =
        (meas - pointTR).dot(normalPrimary) / dirTR.dot(normalPrimary);

    Acts::Vector3 refPointBL = pointBL + dirBL * ddBL;
    Acts::Vector3 refPointTL = pointTL + dirTL * ddTL;

    Acts::Vector3 refPointBR = pointBR + dirBR * ddBR;
    Acts::Vector3 refPointTR = pointTR + dirTR * ddTR;

    geoRefPointMap[id] = {refPointBL, refPointTL, refPointTR, refPointBR};
  }

  std::vector<std::size_t> seedSourceLinkIdxs;
  for (std::size_t i = 0; i < sourceLinksIndices.size(); i++) {
    std::size_t idx = sourceLinksIndices.at(i);
    const auto &ssl = sourceLinks.at(idx).get<SimpleSourceLink>();
    Acts::Vector3 meas = ssl.parametersGlob() - shift;
    const auto &id = ssl.geometryId();

    const Acts::Vector3 &point0 = geoRefPointMap.at(id).at(0);
    const Acts::Vector3 &point1 = geoRefPointMap.at(id).at(1);
    const Acts::Vector3 &point2 = geoRefPointMap.at(id).at(2);
    const Acts::Vector3 &point3 = geoRefPointMap.at(id).at(3);

    double maxLong = std::max(
        {point0(longIdx), point1(longIdx), point2(longIdx), point3(longIdx)});
    double minLong = std::min(
        {point0(longIdx), point1(longIdx), point2(longIdx), point3(longIdx)});

    double maxShort = std::max({point0(shortIdx), point1(shortIdx),
                                point2(shortIdx), point3(shortIdx)});
    double minShort = std::min({point0(shortIdx), point1(shortIdx),
                                point2(shortIdx), point3(shortIdx)});

    if (meas(longIdx) < maxLong && meas(longIdx) > minLong &&
        meas(shortIdx) < maxShort && meas(shortIdx) > minShort) {
      seedSourceLinkIdxs.push_back(idx);
    }
  }

  return seedSourceLinkIdxs;
}

HoughTransformSeeder::HoughTransformSeeder(const Config &cfg) : m_cfg(cfg) {}

std::vector<HoughTransformSeeder::HTSeed> HoughTransformSeeder::findSeeds(
    const Acts::GeometryContext &gctx,
    const std::vector<Acts::SourceLink> &sourceLinks,
    const std::vector<std::size_t> &sourceLinksIndices, const Options &opt) {
  std::vector<HoughTransformSeeder::HTSeed> seeds;

  std::size_t primaryIdx = m_cfg.primaryIdx;
  std::size_t longIdx = m_cfg.longIdx;
  std::size_t shortIdx = m_cfg.shortIdx;

  double maxThetaLong =
      std::atan(opt.boundBoxHalfLong / opt.boundBoxHalfPrimary);
  double maxThetaShort =
      std::atan(opt.boundBoxHalfShort / opt.boundBoxHalfPrimary);

  double maxRhoLong = 2 * (opt.boundBoxHalfLong * std::sin(maxThetaLong) +
                           opt.boundBoxHalfPrimary * std::cos(maxThetaLong));
  double maxRhoShort = 2 * (opt.boundBoxHalfShort * std::sin(maxThetaShort) +
                            opt.boundBoxHalfPrimary * std::cos(maxThetaShort));

  double deltaThetaLong = M_PI / opt.nCellsThetaLong;
  double deltaRhoLong = 2 * maxRhoLong / opt.nCellsRhoLong;

  double deltaThetaShort = M_PI / opt.nCellsThetaShort;
  double deltaRhoShort = 2 * maxRhoShort / opt.nCellsRhoShort;

  Acts::Vector3 shift;
  shift(primaryIdx) = opt.boundBoxCenterPrimary;
  shift(longIdx) = opt.boundBoxCenterLong;
  shift(shortIdx) = opt.boundBoxCenterShort;

  for (const auto &[id, surf] : opt.surfaceMap) {
    m_geoPosMap.insert({id, (surf->center(gctx) - shift)(primaryIdx)});
  }

  VotingMap votingMap;
  votingMap.reserve(sourceLinksIndices.size());
  fillVotingMap(votingMap, sourceLinks, sourceLinksIndices, opt, shift,
                deltaThetaShort, deltaRhoShort, deltaThetaLong, deltaRhoLong,
                maxRhoLong, maxRhoShort);

  for (auto bIt = votingMap.begin(); bIt != votingMap.end(); ++bIt) {
    const auto &[cell, count] = *bIt;
    std::uint16_t lBound = 0;
    std::uint16_t rBound = 0;
    if (count < opt.minXCount) {
      continue;
    } else {
      lBound = 2;
      rBound = 2;
    }
    auto [nThetaLong, nRhoLong, nThetaShort, nRhoShort] = cell;

    // --------------------------------------------------
    double thetaLongBL = deltaThetaLong * nThetaLong - M_PI_2;
    double rhoLongBL = deltaRhoLong * nRhoLong - maxRhoLong;

    double thetaShortBL = deltaThetaShort * nThetaShort - M_PI_2;
    double rhoShortBL = deltaRhoShort * nRhoShort - maxRhoShort;

    double aLongBL = -1.0 / std::tan(thetaLongBL);
    double bLongBL = rhoLongBL / std::sin(thetaLongBL);

    double aShortBL = -1.0 / std::tan(thetaShortBL);
    double bShortBL = rhoShortBL / std::sin(thetaShortBL);

    Acts::Vector3 dirBL;
    dirBL(primaryIdx) = 1;
    dirBL(longIdx) = aLongBL;
    dirBL(shortIdx) = aShortBL;
    dirBL.normalize();

    Acts::Vector3 pointBL;
    pointBL(primaryIdx) = 0;
    pointBL(longIdx) = bLongBL;
    pointBL(shortIdx) = bShortBL;

    // --------------------------------------------------
    double thetaLongTL = deltaThetaLong * nThetaLong - M_PI_2;
    double rhoLongTL = deltaRhoLong * (nRhoLong + 1) - maxRhoLong;

    double thetaShortTL = deltaThetaShort * nThetaShort - M_PI_2;
    double rhoShortTL = deltaRhoShort * (nRhoShort + 1) - maxRhoShort;

    double aLongTL = -1.0 / std::tan(thetaLongTL);
    double bLongTL = rhoLongTL / std::sin(thetaLongTL);

    double aShortTL = -1.0 / std::tan(thetaShortTL);
    double bShortTL = rhoShortTL / std::sin(thetaShortTL);

    Acts::Vector3 dirTL;
    dirTL(primaryIdx) = 1;
    dirTL(longIdx) = aLongTL;
    dirTL(shortIdx) = aShortTL;
    dirTL.normalize();

    Acts::Vector3 pointTL;
    pointTL(primaryIdx) = 0;
    pointTL(longIdx) = bLongTL;
    pointTL(shortIdx) = bShortTL;

    // --------------------------------------------------
    double thetaLongBR = deltaThetaLong * (nThetaLong + 1) - M_PI_2;
    double rhoLongBR = deltaRhoLong * nRhoLong - maxRhoLong;

    double thetaShortBR = deltaThetaShort * (nThetaShort + 1) - M_PI_2;
    double rhoShortBR = deltaRhoShort * nRhoShort - maxRhoShort;

    double aLongBR = -1.0 / std::tan(thetaLongBR);
    double bLongBR = rhoLongBR / std::sin(thetaLongBR);

    double aShortBR = -1.0 / std::tan(thetaShortBR);
    double bShortBR = rhoShortBR / std::sin(thetaShortBR);

    Acts::Vector3 dirBR;
    dirBR(primaryIdx) = 1;
    dirBR(longIdx) = aLongBR;
    dirBR(shortIdx) = aShortBR;
    dirBR.normalize();

    Acts::Vector3 pointBR;
    pointBR(primaryIdx) = 0;
    pointBR(longIdx) = bLongBR;
    pointBR(shortIdx) = bShortBR;

    // --------------------------------------------------
    double thetaLongTR = deltaThetaLong * (nThetaLong + 1) - M_PI_2;
    double rhoLongTR = deltaRhoLong * (nRhoLong + 1) - maxRhoLong;

    double thetaShortTR = deltaThetaShort * (nThetaShort + 1) - M_PI_2;
    double rhoShortTR = deltaRhoShort * (nRhoShort + 1) - maxRhoShort;

    double aLongTR = -1.0 / std::tan(thetaLongTR);
    double bLongTR = rhoLongTR / std::sin(thetaLongTR);

    double aShortTR = -1.0 / std::tan(thetaShortTR);
    double bShortTR = rhoShortTR / std::sin(thetaShortTR);

    Acts::Vector3 dirTR;
    dirTR(primaryIdx) = 1;
    dirTR(longIdx) = aLongTR;
    dirTR(shortIdx) = aShortTR;
    dirTR.normalize();

    Acts::Vector3 pointTR;
    pointTR(primaryIdx) = 0;
    pointTR(longIdx) = bLongTR;
    pointTR(shortIdx) = bShortTR;

    std::vector<std::size_t> seedSourceLinksIndices = findLineSourceLinks(
        sourceLinks, sourceLinksIndices, pointBL, dirBL, pointTL, dirTL,
        pointBR, dirBR, pointTR, dirTR, shift);
    if (seedSourceLinksIndices.empty()) {
      continue;
    }

    double thetaLong = deltaThetaLong * (nThetaLong + 0.5) - M_PI_2;
    double rhoLong = deltaRhoLong * (nRhoLong + 0.5) - maxRhoLong;

    double thetaShort = deltaThetaShort * (nThetaShort + 0.5) - M_PI_2;
    double rhoShort = deltaRhoShort * (nRhoShort + 0.5) - maxRhoShort;

    double aLong = -1.0 / std::tan(thetaLong);
    double bLong = rhoLong / std::sin(thetaLong);

    double aShort = -1.0 / std::tan(thetaShort);
    double bShort = rhoShort / std::sin(thetaShort);

    Acts::Vector3 dir;
    dir(primaryIdx) = 1;
    dir(longIdx) = aLong;
    dir(shortIdx) = aShort;
    dir.normalize();

    Acts::Vector3 point;
    point(primaryIdx) = 0;
    point(longIdx) = bLong;
    point(primaryIdx) = bShort;
    point += shift;

    seeds.emplace_back(std::move(point), std::move(dir),
                       std::move(seedSourceLinksIndices));
  }
  return seeds;
}

void HoughTransformSeeder::fillVotingMap(
    VotingMap &votingMap, const std::vector<Acts::SourceLink> &sourceLinks,
    const std::vector<std::size_t> &sourceLinksIndices, const Options &opt,
    const Acts::Vector3 &shift, double deltaThetaShort, double deltaRhoShort,
    double deltaThetaLong, double deltaRhoLong, double maxRhoLong,
    double maxRhoShort) {
  std::unordered_map<int, std::vector<std::size_t>> clusters;
  clusters.reserve(sourceLinksIndices.size());
  for (std::size_t i = 0; i < sourceLinksIndices.size(); i++) {
    std::size_t idx = sourceLinksIndices.at(i);
    int geoId =
        sourceLinks.at(idx).get<SimpleSourceLink>().geometryId().sensitive();
    clusters[geoId].push_back(idx);
  }
  std::vector<int> geoIds;
  geoIds.reserve(clusters.size());
  for (const auto &[geoId, sls] : clusters) {
    geoIds.push_back(geoId);
  }

  std::size_t primaryIdx = m_cfg.primaryIdx;
  std::size_t longIdx = m_cfg.longIdx;
  std::size_t shortIdx = m_cfg.shortIdx;
  for (std::size_t i = 0; i < geoIds.size() - 1; i++) {
    int fId = geoIds.at(i);
    const auto &fClusters = clusters.at(fId);
    for (std::size_t j = i + 1; j < geoIds.size(); j++) {
      int sId = geoIds.at(j);
      const auto &sClusters = clusters.at(sId);

      for (std::size_t n = 0; n < fClusters.size(); n++) {
        Acts::Vector3 flPoint = sourceLinks.at(fClusters.at(n))
                                    .get<SimpleSourceLink>()
                                    .parametersGlob() -
                                shift;
        double flPrimary = flPoint(primaryIdx);
        double flLong = flPoint(longIdx);
        double flShort = flPoint(shortIdx);
        for (std::size_t m = 0; m < sClusters.size(); m++) {
          Acts::Vector3 slPoint = sourceLinks.at(sClusters.at(m))
                                      .get<SimpleSourceLink>()
                                      .parametersGlob() -
                                  shift;
          double diffPrimary = slPoint(primaryIdx) - flPrimary;

          double thetaLong =
              std::atan(diffPrimary / (flLong - slPoint(longIdx)));
          double rhoLong =
              flLong * std::sin(thetaLong) + flPrimary * std::cos(thetaLong);

          double thetaShort =
              std::atan(diffPrimary / (flShort - slPoint(shortIdx)));
          double rhoShort =
              flShort * std::sin(thetaShort) + flPrimary * std::cos(thetaShort);

          std::uint16_t cellThetaLong =
              std::ceil((thetaLong + M_PI_2) / deltaThetaLong) - 1;
          std::uint16_t cellRhoLong =
              std::ceil((rhoLong + maxRhoLong) / deltaRhoLong) - 1;
          std::uint16_t cellThetaShort =
              std::ceil((thetaShort + M_PI_2) / deltaThetaShort) - 1;
          std::uint16_t cellRhoShort =
              std::ceil((rhoShort + maxRhoShort) / deltaRhoShort) - 1;

          votingMap[{cellThetaLong, cellRhoLong, cellThetaShort,
                     cellRhoShort}]++;
        }
      }
    }
  }
}
