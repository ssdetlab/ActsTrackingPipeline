#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

namespace {

/// @brief check the Hough space axis range and binning
///
/// @param min lower edge of the axis
/// @param max upper edge of the axis
/// @param nCells number of cells along the axis
/// @param name axis name for the error message
void checkHoughAxis(double min, double max, std::size_t nCells,
                    const std::string &name) {
  if (!(min < max)) {
    throw std::invalid_argument("Invalid Hough space " + name + " range");
  }
  if (nCells == 0 || nCells > std::numeric_limits<std::uint16_t>::max()) {
    throw std::invalid_argument("Invalid Hough space " + name + " binning");
  }
}

/// @brief get the Hough space cell index along an axis
///
/// @param value value to bin
/// @param min lower edge of the axis
/// @param delta cell size along the axis
/// @param nCells number of cells along the axis
///
/// @return cell index, -1 if the value is outside of [min, min + nCells * delta)
std::int64_t getCellIdx(double value, double min, double delta,
                        std::size_t nCells) {
  double idx = std::floor((value - min) / delta);
  if (!(idx >= 0) || idx >= static_cast<double>(nCells)) {
    return -1;
  }
  return static_cast<std::int64_t>(idx);
}

}  // namespace

std::vector<std::size_t> HoughTransformSeeder::findLineSourceLinks(
    const std::vector<Acts::SourceLink> &sourceLinks,
    const std::vector<std::size_t> &sourceLinksIndices,
    const Acts::Vector3 &shift, double minSlopeLong, double maxSlopeLong,
    double minInterceptLong, double maxInterceptLong, double minSlopeShort,
    double maxSlopeShort, double minInterceptShort, double maxInterceptShort) {
  std::size_t primaryIdx = m_cfg.primaryIdx;
  std::size_t longIdx = m_cfg.longIdx;
  std::size_t shortIdx = m_cfg.shortIdx;

  std::vector<std::size_t> seedSourceLinkIdxs;
  for (std::size_t i = 0; i < sourceLinksIndices.size(); i++) {
    std::size_t idx = sourceLinksIndices.at(i);
    const auto &ssl = sourceLinks.at(idx).get<SimpleSourceLink>();
    Acts::Vector3 meas = ssl.parametersGlob() - shift;
    double primary = meas(primaryIdx);

    // The transverse coordinate of a line, slope * primary + intercept,
    // is linear in both the slope and the intercept, so its extrema over
    // the cell are reached at the cell corners
    double minLong = std::min(minSlopeLong * primary, maxSlopeLong * primary) +
                     minInterceptLong;
    double maxLong = std::max(minSlopeLong * primary, maxSlopeLong * primary) +
                     maxInterceptLong;

    double minShort =
        std::min(minSlopeShort * primary, maxSlopeShort * primary) +
        minInterceptShort;
    double maxShort =
        std::max(minSlopeShort * primary, maxSlopeShort * primary) +
        maxInterceptShort;

    if (meas(longIdx) >= minLong && meas(longIdx) <= maxLong &&
        meas(shortIdx) >= minShort && meas(shortIdx) <= maxShort) {
      seedSourceLinkIdxs.push_back(idx);
    }
  }

  return seedSourceLinkIdxs;
}

HoughTransformSeeder::HoughTransformSeeder(const Config &cfg) : m_cfg(cfg) {}

std::vector<HoughTransformSeeder::HTSeed> HoughTransformSeeder::findSeeds(
    const Acts::GeometryContext & /*gctx*/,
    const std::vector<Acts::SourceLink> &sourceLinks,
    const std::vector<std::size_t> &sourceLinksIndices, const Options &opt) {
  checkHoughAxis(opt.minSlopeLong, opt.maxSlopeLong, opt.nCellsSlopeLong,
                 "long slope");
  checkHoughAxis(opt.minInterceptLong, opt.maxInterceptLong,
                 opt.nCellsInterceptLong, "long intercept");
  checkHoughAxis(opt.minSlopeShort, opt.maxSlopeShort, opt.nCellsSlopeShort,
                 "short slope");
  checkHoughAxis(opt.minInterceptShort, opt.maxInterceptShort,
                 opt.nCellsInterceptShort, "short intercept");

  std::vector<HoughTransformSeeder::HTSeed> seeds;

  std::size_t primaryIdx = m_cfg.primaryIdx;
  std::size_t longIdx = m_cfg.longIdx;
  std::size_t shortIdx = m_cfg.shortIdx;

  double deltaSlopeLong =
      (opt.maxSlopeLong - opt.minSlopeLong) / opt.nCellsSlopeLong;
  double deltaInterceptLong =
      (opt.maxInterceptLong - opt.minInterceptLong) / opt.nCellsInterceptLong;

  double deltaSlopeShort =
      (opt.maxSlopeShort - opt.minSlopeShort) / opt.nCellsSlopeShort;
  double deltaInterceptShort = (opt.maxInterceptShort - opt.minInterceptShort) /
                               opt.nCellsInterceptShort;

  Acts::Vector3 shift;
  shift(primaryIdx) = opt.boundBoxCenterPrimary;
  shift(longIdx) = opt.boundBoxCenterLong;
  shift(shortIdx) = opt.boundBoxCenterShort;

  VotingMap votingMap;
  votingMap.reserve(sourceLinksIndices.size());
  fillVotingMap(votingMap, sourceLinks, sourceLinksIndices, opt, shift,
                deltaSlopeLong, deltaInterceptLong, deltaSlopeShort,
                deltaInterceptShort);

  for (auto bIt = votingMap.begin(); bIt != votingMap.end(); ++bIt) {
    const auto &[cell, count] = *bIt;
    if (count < opt.minXCount) {
      continue;
    }
    auto [nSlopeLong, nInterceptLong, nSlopeShort, nInterceptShort] = cell;

    // Cell edges
    double minSlopeLong = opt.minSlopeLong + deltaSlopeLong * nSlopeLong;
    double maxSlopeLong = minSlopeLong + deltaSlopeLong;
    double minInterceptLong =
        opt.minInterceptLong + deltaInterceptLong * nInterceptLong;
    double maxInterceptLong = minInterceptLong + deltaInterceptLong;

    double minSlopeShort = opt.minSlopeShort + deltaSlopeShort * nSlopeShort;
    double maxSlopeShort = minSlopeShort + deltaSlopeShort;
    double minInterceptShort =
        opt.minInterceptShort + deltaInterceptShort * nInterceptShort;
    double maxInterceptShort = minInterceptShort + deltaInterceptShort;

    std::vector<std::size_t> seedSourceLinksIndices = findLineSourceLinks(
        sourceLinks, sourceLinksIndices, shift, minSlopeLong, maxSlopeLong,
        minInterceptLong, maxInterceptLong, minSlopeShort, maxSlopeShort,
        minInterceptShort, maxInterceptShort);
    if (seedSourceLinksIndices.empty()) {
      continue;
    }

    // Cell center
    Acts::Vector3 dir;
    dir(primaryIdx) = 1;
    dir(longIdx) = minSlopeLong + 0.5 * deltaSlopeLong;
    dir(shortIdx) = minSlopeShort + 0.5 * deltaSlopeShort;
    dir.normalize();

    Acts::Vector3 point;
    point(primaryIdx) = 0;
    point(longIdx) = minInterceptLong + 0.5 * deltaInterceptLong;
    point(shortIdx) = minInterceptShort + 0.5 * deltaInterceptShort;
    point += shift;

    seeds.emplace_back(std::move(point), std::move(dir),
                       std::move(seedSourceLinksIndices), count);
  }
  return seeds;
}

void HoughTransformSeeder::fillVotingMap(
    VotingMap &votingMap, const std::vector<Acts::SourceLink> &sourceLinks,
    const std::vector<std::size_t> &sourceLinksIndices, const Options &opt,
    const Acts::Vector3 &shift, double deltaSlopeLong,
    double deltaInterceptLong, double deltaSlopeShort,
    double deltaInterceptShort) {
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
  for (std::size_t i = 0; i + 1 < geoIds.size(); i++) {
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

          // The Hough space lines of the measurements at the same
          // primary position are parallel
          if (diffPrimary == 0) {
            continue;
          }

          // Intersection of the Hough space lines
          // intercept = transverse - slope * primary
          double slopeLong = (slPoint(longIdx) - flLong) / diffPrimary;
          double interceptLong = flLong - slopeLong * flPrimary;

          double slopeShort = (slPoint(shortIdx) - flShort) / diffPrimary;
          double interceptShort = flShort - slopeShort * flPrimary;

          std::int64_t cellSlopeLong = getCellIdx(
              slopeLong, opt.minSlopeLong, deltaSlopeLong, opt.nCellsSlopeLong);
          std::int64_t cellInterceptLong =
              getCellIdx(interceptLong, opt.minInterceptLong,
                         deltaInterceptLong, opt.nCellsInterceptLong);
          std::int64_t cellSlopeShort =
              getCellIdx(slopeShort, opt.minSlopeShort, deltaSlopeShort,
                         opt.nCellsSlopeShort);
          std::int64_t cellInterceptShort =
              getCellIdx(interceptShort, opt.minInterceptShort,
                         deltaInterceptShort, opt.nCellsInterceptShort);
          if (cellSlopeLong < 0 || cellInterceptLong < 0 ||
              cellSlopeShort < 0 || cellInterceptShort < 0) {
            continue;
          }

          votingMap[{static_cast<std::uint16_t>(cellSlopeLong),
                     static_cast<std::uint16_t>(cellInterceptLong),
                     static_cast<std::uint16_t>(cellSlopeShort),
                     static_cast<std::uint16_t>(cellInterceptShort)}]++;
        }
      }
    }
  }
}
