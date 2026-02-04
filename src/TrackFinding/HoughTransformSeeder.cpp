#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <tuple>
#include <vector>

#include <Eigen/src/Core/Matrix.h>
#include <omp.h>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/Geometry/detail/BinningValueUtils.hpp"

void constructTracks(const std::shared_ptr<IdxTree::Node> &root,
                     std::vector<int> &track,
                     std::vector<std::vector<int>> &tracks) {
  track.push_back(root->m_idx);
  if (root->children.size() == 0) {
    tracks.push_back(track);
    track.pop_back();
    return;
  }
  for (auto &child : root->children) {
    constructTracks(child, track, tracks);
  }
  track.pop_back();
}

Eigen::MatrixXd HoughTransformSeeder::constructCov(
    const std::vector<SourceLinkRef> &sourceLinks, const Acts::Vector3 &dir,
    const Options &opt) const {
  Eigen::MatrixXd D = Eigen::MatrixXd::Constant(2 * sourceLinks.size(),
                                                2 * sourceLinks.size(), 0);
  double mcsVarFactor = std::pow(
      opt.primaryInterchipDistance * dir(m_cfg.primaryIdx) * opt.thetaRms, 2);
  for (std::size_t i = 0; i < sourceLinks.size(); i++) {
    const auto &ssl = sourceLinks.at(i).get().get<SimpleSourceLink>();
    double mcpVar = mcsVarFactor * i * (1 + 3 * i + 2 * i * i) / 6.0;
    Acts::SquareMatrix2 varMat =
        ssl.covariance() + Acts::SquareMatrix2::Identity() * mcpVar;
    D.block(i * 2, i * 2, 2, 2) = varMat;
    for (std::size_t j = 0; j < i; j++) {
      double covPiPj = mcsVarFactor * j * (j + 1) * (1 - j + 3 * i) / 6.0;
      Acts::SquareMatrix2 covMat = Acts::SquareMatrix2::Identity() * covPiPj;
      D.block(i * 2, j * 2, 2, 2) = covMat;
      D.block(j * 2, i * 2, 2, 2) = covMat;
    }
  }
  return std::move(D);
}

double HoughTransformSeeder::globalChi2Fit(
    const std::vector<SourceLinkRef> &sourceLinks, Acts::Vector3 &pos,
    Acts::Vector3 &dir, Acts::ActsSquareMatrix<6> &cov,
    const Options &opt) const {
  Eigen::MatrixXd G = Eigen::MatrixXd::Constant(2 * sourceLinks.size(), 4, 0);

  std::size_t primaryIdx = m_cfg.primaryIdx;
  std::size_t longIdx = m_cfg.longIdx;
  std::size_t shortIdx = m_cfg.shortIdx;

  double mcsVarFactor = std::pow(
      opt.primaryInterchipDistance * dir(primaryIdx) * opt.thetaRms, 2);
  Eigen::VectorXd X(2 * sourceLinks.size());
  for (std::size_t i = 0; i < sourceLinks.size(); i++) {
    const auto &ssl = sourceLinks.at(i).get().get<SimpleSourceLink>();
    const Acts::Vector3 &parameters = ssl.parametersGlob();
    X(2 * i) = parameters(longIdx);
    X(2 * i + 1) = parameters(shortIdx);

    G(2 * i, 0) = 1;
    G(2 * i, 2) = opt.primaryInterchipDistance * i;
    G(2 * i + 1, 1) = 1;
    G(2 * i + 1, 3) = opt.primaryInterchipDistance * i;
  }
  Eigen::LDLT<Eigen::MatrixXd> ldltD(constructCov(sourceLinks, dir, opt));
  Acts::SquareMatrix4 B = G.transpose() * ldltD.solve(G);
  Acts::Vector4 estimates = B.ldlt().solve(G.transpose() * ldltD.solve(X));
  double tLong = estimates(2);
  double tShort = estimates(3);

  pos(primaryIdx) = opt.boundBoxCenterPrimary - m_cfg.boundBoxHalfPrimary;
  pos(longIdx) = estimates(0);
  pos(shortIdx) = estimates(1);

  dir(primaryIdx) = 1.0;
  dir(longIdx) = tLong;
  dir(shortIdx) = tShort;
  dir.normalize();

  double denom = std::pow(1 + tLong * tLong + tShort * tShort, 1.5);

  Acts::ActsMatrix<6, 4> jacToGlob;
  jacToGlob.setZero();
  jacToGlob(1, 0) = 1;
  jacToGlob(2, 1) = 1;

  jacToGlob(3, 2) = -tLong / denom;
  jacToGlob(3, 3) = -tShort / denom;

  jacToGlob(4, 2) = (1 + tShort * tShort) / denom;
  jacToGlob(4, 3) = -tLong * tShort / denom;

  jacToGlob(5, 2) = -tLong * tShort / denom;
  jacToGlob(5, 3) = (1 + tLong * tLong) / denom;

  cov = jacToGlob * B.inverse() * jacToGlob.transpose();
  return (X - G * estimates).transpose() * ldltD.solve(X - G * estimates);
}

std::vector<std::pair<int, int>> HoughTransformSeeder::findLineSourceLinks(
    const std::span<SourceLinkRef> &sourceLinks, const Acts::Vector3 &pointBL,
    const Acts::Vector3 &dirBL, const Acts::Vector3 &pointTR,
    const Acts::Vector3 &dirTR, const Acts::Vector3 &shift) const {
  std::size_t primaryIdx = m_cfg.primaryIdx;
  std::size_t longIdx = m_cfg.longIdx;
  std::size_t shortIdx = m_cfg.shortIdx;

  Acts::Vector3 normalPrimary = detail::indexToDirection(primaryIdx);

  std::vector<std::pair<int, int>> seedSlIdxs;
  for (std::size_t idx = 0; idx < sourceLinks.size(); idx++) {
    const auto &ssl = sourceLinks[idx].get().get<SimpleSourceLink>();
    Acts::Vector3 meas = ssl.parametersGlob() - shift;

    double ddBL =
        (meas - pointBL).dot(normalPrimary) / dirBL.dot(normalPrimary);
    double ddTR =
        (meas - pointTR).dot(normalPrimary) / dirTR.dot(normalPrimary);

    Acts::Vector3 refPointBL = pointBL + dirBL * ddBL;
    Acts::Vector3 refPointTR = pointTR + dirTR * ddTR;

    if (meas(longIdx) < std::max(refPointBL(longIdx), refPointTR(longIdx)) &&
        meas(longIdx) > std::min(refPointBL(longIdx), refPointTR(longIdx)) &&
        meas(shortIdx) < std::max(refPointBL(shortIdx), refPointTR(shortIdx)) &&
        meas(shortIdx) > std::min(refPointBL(shortIdx), refPointTR(shortIdx))) {
      seedSlIdxs.push_back({idx, ssl.geometryId().sensitive()});
    }
  }

  return seedSlIdxs;
}

HoughTransformSeeder::HoughTransformSeeder(const Config &cfg) : m_cfg(cfg) {
  double maxThetaLong =
      std::atan(m_cfg.boundBoxHalfLong / m_cfg.boundBoxHalfPrimary);
  double maxThetaShort =
      std::atan(m_cfg.boundBoxHalfShort / m_cfg.boundBoxHalfPrimary);

  m_maxRhoLong = 2 * (m_cfg.boundBoxHalfLong * std::sin(maxThetaLong) +
                      m_cfg.boundBoxHalfPrimary * std::cos(maxThetaLong));
  m_maxRhoShort = 2 * (m_cfg.boundBoxHalfShort * std::sin(maxThetaShort) +
                       m_cfg.boundBoxHalfPrimary * std::cos(maxThetaShort));

  m_deltaThetaLong = M_PI / m_cfg.nCellsThetaLong;
  m_deltaRhoLong = m_maxRhoLong / m_cfg.nCellsRhoLong;

  m_deltaThetaShort = M_PI / m_cfg.nCellsThetaShort;
  m_deltaRhoShort = m_maxRhoShort / m_cfg.nCellsRhoShort;
}

std::vector<HoughTransformSeeder::HTSeed> HoughTransformSeeder::findSeeds(
    std::span<SourceLinkRef> sourceLinks, const Options &opt) const {
  std::vector<HoughTransformSeeder::HTSeed> seeds;

  std::size_t primaryIdx = m_cfg.primaryIdx;
  std::size_t longIdx = m_cfg.longIdx;
  std::size_t shortIdx = m_cfg.shortIdx;

  Acts::Vector3 shift;
  shift(primaryIdx) = opt.boundBoxCenterPrimary;
  shift(longIdx) = opt.boundBoxCenterLong;
  shift(shortIdx) = opt.boundBoxCenterShort;

  VotingMap votingMap;
  votingMap.reserve(sourceLinks.size());
  fillVotingMap(votingMap, sourceLinks, opt, shift);

  for (std::size_t bIdx = 0; bIdx < votingMap.bucket_count(); ++bIdx) {
    for (auto bIt = votingMap.begin(bIdx); bIt != votingMap.end(bIdx); ++bIt) {
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
      double thetaLongBL =
          m_deltaThetaLong * (nThetaLong - lBound + 1) - M_PI_2;
      double rhoLongBL =
          m_deltaRhoLong * (nRhoLong - lBound + 1) - m_maxRhoLong / 2.0;

      double thetaShortBL =
          m_deltaThetaShort * (nThetaShort - lBound + 1) - M_PI_2;
      double rhoShortBL =
          m_deltaRhoShort * (nRhoShort - lBound + 1) - m_maxRhoShort / 2.0;

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
      double thetaLongTR = m_deltaThetaLong * (nThetaLong + rBound) - M_PI_2;
      double rhoLongTR =
          m_deltaRhoLong * (nRhoLong + rBound) - m_maxRhoLong / 2.0;

      double thetaShortTR = m_deltaThetaShort * (nThetaShort + rBound) - M_PI_2;
      double rhoShortTR =
          m_deltaRhoShort * (nRhoShort + rBound) - m_maxRhoShort / 2.0;

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

      std::vector<std::pair<int, int>> seedSlIdxs = findLineSourceLinks(
          sourceLinks, pointBL, dirBL, pointTR, dirTR, shift);
      if (seedSlIdxs.empty()) {
        continue;
      }
      std::size_t slIdxsSize = seedSlIdxs.size();

      double thetaLong = m_deltaThetaLong * (nThetaLong + 0.5) - M_PI_2;
      double rhoLong = m_deltaRhoLong * (nRhoLong + 0.5) - m_maxRhoLong / 2.0;

      double thetaShort = m_deltaThetaShort * (nThetaShort + 0.5) - M_PI_2;
      double rhoShort =
          m_deltaRhoShort * (nRhoShort + 0.5) - m_maxRhoShort / 2.0;

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

      std::sort(
          seedSlIdxs.begin(), seedSlIdxs.end(),
          [](const auto &a, const auto &b) { return a.second < b.second; });
      auto rootEndIt = std::find_if(
          seedSlIdxs.begin(), seedSlIdxs.end(),
          [&opt](const auto &a) { return (a.second != opt.firstLayerId); });

      // double rc = 0;
      double chi2 = 0;
      for (auto it = seedSlIdxs.begin(); it != rootEndIt; it++) {
        std::vector<int> trackContainer;
        trackContainer.reserve(opt.minSeedSize);

        std::vector<std::vector<int>> splitSeedSlIdxs;
        IdxTree idxTree(seedSlIdxs, it, rootEndIt);
        constructTracks(idxTree.m_root, trackContainer, splitSeedSlIdxs);
        for (const auto &seed : splitSeedSlIdxs) {
          if (seed.size() < opt.minSeedSize || seed.size() > opt.maxSeedSize) {
            continue;
          }
          std::vector<Acts::SourceLink> seedSourceLinks;
          std::vector<SourceLinkRef> seedSourceLinksRefs;
          seedSourceLinks.reserve(slIdxsSize);
          seedSourceLinksRefs.reserve(slIdxsSize);

          Acts::Vector3 newDir = dir;
          Acts::Vector3 newPoint = point;
          Acts::ActsSquareMatrix<6> newCov;

          for (auto idx : seed) {
            seedSourceLinksRefs.push_back(std::cref(sourceLinks[idx]));
            seedSourceLinks.push_back(sourceLinks[idx]);
          }
          for (std::size_t l = 0; l < m_cfg.nGLSIterations; l++) {
            chi2 = globalChi2Fit(seedSourceLinksRefs, newPoint, newDir, newCov,
                                 opt);
          }
          if (chi2 > opt.maxChi2) {
            continue;
          }
          seeds.emplace_back(std::move(newPoint), std::move(newDir),
                             std::move(newCov), std::move(seedSourceLinks));
        }
      }
    }
  }
  return seeds;
}

void HoughTransformSeeder::fillVotingMap(VotingMap &votingMap,
                                         std::span<SourceLinkRef> sourceLinks,
                                         const Options &opt,
                                         const Acts::Vector3 &shift) const {
  std::unordered_map<int, std::vector<SourceLinkRef>> clusters;
  std::vector<int> geoIds;
  clusters.reserve(sourceLinks.size());
  for (auto sl : sourceLinks) {
    int geoId = sl.get().get<SimpleSourceLink>().geometryId().sensitive();
    clusters[geoId].push_back(sl);
  }
  if (clusters.size() < opt.minSeedSize) {
    return;
  }
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
        Acts::Vector3 flPoint =
            fClusters[n].get().get<SimpleSourceLink>().parametersGlob() - shift;
        double flPrimary = flPoint(primaryIdx);
        double flLong = flPoint(longIdx);
        double flShort = flPoint(shortIdx);
        for (std::size_t m = 0; m < sClusters.size(); m++) {
          Acts::Vector3 slPoint =
              sClusters[m].get().get<SimpleSourceLink>().parametersGlob() -
              shift;
          double diffPrimary = slPoint.x() - flPrimary;

          double thetaLong =
              std::atan(diffPrimary / (flLong - slPoint(longIdx)));
          double rhoLong =
              flLong * std::sin(thetaLong) + flPrimary * std::cos(thetaLong);

          double thetaShort =
              std::atan(diffPrimary / (flShort - slPoint(shortIdx)));
          double rhoShort =
              flShort * std::sin(thetaShort) + flPrimary * std::cos(thetaShort);

          std::uint16_t cellThetaLong =
              std::ceil((thetaLong + M_PI_2) / m_deltaThetaLong) - 1;
          std::uint16_t cellRhoLong =
              std::ceil((rhoLong + m_maxRhoLong / 2.0) / m_deltaRhoLong) - 1;
          std::uint16_t cellThetaShort =
              std::ceil((thetaShort + M_PI_2) / m_deltaThetaShort) - 1;
          std::uint16_t cellRhoShort =
              std::ceil((rhoShort + m_maxRhoShort / 2.0) / m_deltaRhoShort) - 1;

          votingMap[{cellThetaLong, cellRhoLong, cellThetaShort,
                     cellRhoShort}]++;
        }
      }
    }
  }
}
