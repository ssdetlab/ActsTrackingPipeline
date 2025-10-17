#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <tuple>
#include <vector>

#include <omp.h>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

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

double computeDistance(const Acts::Vector3 &point, const Acts::Vector3 &dir,
                       const Acts::Vector3 &meas) {
  return (point + dir.dot(meas - point) * dir - meas).norm();
}

double HoughTransformSeeder::orthogonalLeastSquares(
    const std::vector<SourceLinkRef> &sourceLinks, Acts::Vector3 &a,
    Acts::Vector3 &b, const Acts::Vector3 &shift) const {
  Acts::Vector3 meanVector(0, 0, 0);

  Eigen::MatrixXf points = Eigen::MatrixXf::Constant(sourceLinks.size(), 3, 0);
  for (int i = 0; i < points.rows(); i++) {
    Acts::Vector3 parameters =
        sourceLinks.at(i).get().get<SimpleSourceLink>().parametersGlob();
    points(i, 0) = parameters.x();
    points(i, 1) = parameters.y();
    points(i, 2) = parameters.z();

    meanVector += parameters;
  }
  meanVector /= sourceLinks.size();
  a = meanVector - shift;

  Eigen::MatrixXf centered = points.rowwise() - points.colwise().mean();
  Eigen::MatrixXf scatter = (centered.adjoint() * centered);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig(scatter);
  Eigen::MatrixXf eigvecs = eig.eigenvectors();

  b[0] = eigvecs(0, 2);
  b[1] = eigvecs(1, 2);
  b[2] = eigvecs(2, 2);
  return eig.eigenvalues()(2);
}

std::vector<std::pair<int, int>> HoughTransformSeeder::findLineSourceLinks(
    const std::span<SourceLinkRef> &sourceLinks, const Acts::Vector3 &pointBL,
    const Acts::Vector3 &dirBL, const Acts::Vector3 &pointTR,
    const Acts::Vector3 &dirTR, const Acts::Vector3 &shift) const {
  std::vector<std::pair<int, int>> seedSlIdxs;
  for (std::size_t idx = 0; idx < sourceLinks.size(); idx++) {
    const auto &ssl = sourceLinks[idx].get().get<SimpleSourceLink>();
    Acts::Vector3 meas = ssl.parametersGlob() - shift;

    double ddBL = (meas - pointBL).dot(Acts::Vector3::UnitX()) /
                  dirBL.dot(Acts::Vector3::UnitX());
    double ddTR = (meas - pointTR).dot(Acts::Vector3::UnitX()) /
                  dirTR.dot(Acts::Vector3::UnitX());

    Acts::Vector3 refPointBL = pointBL + dirBL * ddBL;
    Acts::Vector3 refPointTR = pointTR + dirTR * ddTR;

    if (meas.y() < std::max(refPointBL.y(), refPointTR.y()) &&
        meas.y() > std::min(refPointBL.y(), refPointTR.y()) &&
        meas.z() < std::max(refPointBL.z(), refPointTR.z()) &&
        meas.z() > std::min(refPointBL.z(), refPointTR.z())) {
      seedSlIdxs.push_back({idx, ssl.geometryId().sensitive()});
    }
  }

  return seedSlIdxs;
}

HoughTransformSeeder::HoughTransformSeeder(const Config &cfg) : m_cfg(cfg) {
  double maxThetaY = std::atan(m_cfg.boundBoxHalfY / m_cfg.boundBoxHalfX);
  double maxThetaZ = std::atan(m_cfg.boundBoxHalfZ / m_cfg.boundBoxHalfX);

  m_maxRhoY = 2 * (m_cfg.boundBoxHalfY * std::sin(maxThetaY) +
                   m_cfg.boundBoxHalfX * std::cos(maxThetaY));
  m_maxRhoZ = 2 * (m_cfg.boundBoxHalfZ * std::sin(maxThetaZ) +
                   m_cfg.boundBoxHalfX * std::cos(maxThetaZ));

  m_deltaThetaZ = M_PI / m_cfg.nCellsThetaX;
  m_deltaRhoZ = m_maxRhoZ / m_cfg.nCellsRhoX;

  m_deltaThetaY = M_PI / m_cfg.nCellsThetaY;
  m_deltaRhoY = m_maxRhoY / m_cfg.nCellsRhoY;
}

std::vector<HoughTransformSeeder::HTSeed> HoughTransformSeeder::findSeeds(
    std::span<SourceLinkRef> sourceLinks, const Options &opt) const {
  std::vector<HoughTransformSeeder::HTSeed> seeds;

  Acts::Vector3 shift(opt.boundBoxCenterX + m_cfg.boundBoxHalfX,
                      opt.boundBoxCenterY, opt.boundBoxCenterZ);
  VotingMap votingMap;
  votingMap.reserve(sourceLinks.size());
  fillVotingMap(votingMap, sourceLinks, 32, opt, shift);

#pragma omp parallel for num_threads(32)
  for (std::size_t bIdx = 0; bIdx < votingMap.bucket_count(); ++bIdx) {
    for (auto bIt = votingMap.begin(bIdx); bIt != votingMap.end(bIdx); ++bIt) {
      const auto &[cell, count] = *bIt;
      std::uint16_t lBound = 0;
      std::uint16_t rBound = 0;
      if (count < opt.minCount) {
        continue;
      } else {
        lBound = 2;
        rBound = 2;
      }
      auto [nThetaY, nRhoY, nThetaZ, nRhoZ] = cell;

      // --------------------------------------------------
      double thetaYBL = m_deltaThetaY * (nThetaY - lBound + 1) - M_PI_2;
      double rhoYBL = m_deltaRhoY * (nRhoY - lBound + 1) - m_maxRhoY / 2.0;

      double thetaZBL = m_deltaThetaZ * (nThetaZ - lBound + 1) - M_PI_2;
      double rhoZBL = m_deltaRhoZ * (nRhoZ - lBound + 1) - m_maxRhoZ / 2.0;

      double ayBL = -1.0 / std::tan(thetaYBL);
      double byBL = rhoYBL / std::sin(thetaYBL);

      double azBL = -1.0 / std::tan(thetaZBL);
      double bzBL = rhoZBL / std::sin(thetaZBL);

      Acts::Vector3 dirBL = Acts::Vector3(1, ayBL, azBL).normalized();
      Acts::Vector3 pointBL = Acts::Vector3(0, byBL, bzBL);

      // --------------------------------------------------
      double thetaYTR = m_deltaThetaY * (nThetaY + rBound) - M_PI_2;
      double rhoYTR = m_deltaRhoY * (nRhoY + rBound) - m_maxRhoY / 2.0;

      double thetaZTR = m_deltaThetaZ * (nThetaZ + rBound) - M_PI_2;
      double rhoZTR = m_deltaRhoZ * (nRhoZ + rBound) - m_maxRhoZ / 2.0;

      double ayTR = -1.0 / std::tan(thetaYTR);
      double byTR = rhoYTR / std::sin(thetaYTR);

      double azTR = -1.0 / std::tan(thetaZTR);
      double bzTR = rhoZTR / std::sin(thetaZTR);

      Acts::Vector3 dirTR = Acts::Vector3(1, ayTR, azTR).normalized();
      Acts::Vector3 pointTR = Acts::Vector3(0, byTR, bzTR);

      std::vector<std::pair<int, int>> seedSlIdxs = findLineSourceLinks(
          sourceLinks, pointBL, dirBL, pointTR, dirTR, shift);
      if (seedSlIdxs.empty()) {
        continue;
      }
      std::size_t slIdxsSize = seedSlIdxs.size();

      double thetaY = m_deltaThetaY * (nThetaY + 0.5) - M_PI_2;
      double rhoY = m_deltaRhoY * (nRhoY + 0.5) - m_maxRhoY / 2.0;

      double thetaZ = m_deltaThetaZ * (nThetaZ + 0.5) - M_PI_2;
      double rhoZ = m_deltaRhoZ * (nRhoZ + 0.5) - m_maxRhoZ / 2.0;

      double ay = -1.0 / std::tan(thetaY);
      double by = rhoY / std::sin(thetaY);

      double az = -1.0 / std::tan(thetaZ);
      double bz = rhoZ / std::sin(thetaZ);

      Acts::Vector3 dir = Acts::Vector3(1, ay, az).normalized();
      Acts::Vector3 point = Acts::Vector3(0, by, bz);

      std::sort(
          seedSlIdxs.begin(), seedSlIdxs.end(),
          [](const auto &a, const auto &b) { return a.second < b.second; });
      auto rootEndIt = std::find_if(
          seedSlIdxs.begin(), seedSlIdxs.end(),
          [&opt](const auto &a) { return (a.second != opt.firstLayerId); });

      double rc = 0;
      for (auto it = seedSlIdxs.begin(); it != rootEndIt; it++) {
        std::vector<int> trackContainer;
        trackContainer.reserve(m_cfg.minSeedSize);

        std::vector<std::vector<int>> splitSeedSlIdxs;
        IdxTree idxTree(seedSlIdxs, it, rootEndIt);
        constructTracks(idxTree.m_root, trackContainer, splitSeedSlIdxs);
        for (const auto &seed : splitSeedSlIdxs) {
          if (seed.size() < m_cfg.minSeedSize ||
              seed.size() > m_cfg.maxSeedSize) {
            continue;
          }
          std::vector<Acts::SourceLink> seedSourceLinks;
          std::vector<SourceLinkRef> seedSourceLinksRefs;
          seedSourceLinks.reserve(slIdxsSize);
          seedSourceLinksRefs.reserve(slIdxsSize);

          Acts::Vector3 newDir = dir;
          Acts::Vector3 newPoint = point;

          for (auto idx : seed) {
            seedSourceLinksRefs.push_back(std::cref(sourceLinks[idx]));
            seedSourceLinks.push_back(sourceLinks[idx]);
          }
          for (int l = 0; l < 2; l++) {
            rc = orthogonalLeastSquares(seedSourceLinksRefs, newPoint, newDir,
                                        shift);
          }
          if (rc == 0) {
            continue;
          }
          double chi2 = 0;
          for (auto sl : seedSourceLinksRefs) {
            double dist = computeDistance(
                newPoint, newDir,
                sl.get().get<SimpleSourceLink>().parametersGlob() - shift);
            chi2 += dist * dist;
          }
          if (chi2 > 1e-2) {
            continue;
          }
#pragma omp critical
          {
            seeds.emplace_back(std::move(newPoint + shift), std::move(newDir),
                               std::move(seedSourceLinks));
          }
        }
      }
    }
  }
  return seeds;
}

void HoughTransformSeeder::fillVotingMap(VotingMap &votingMap,
                                         std::span<SourceLinkRef> sourceLinks,
                                         int nThreads, const Options &opt,
                                         const Acts::Vector3 &shift) const {
  std::unordered_map<int, std::vector<SourceLinkRef>> clusters;
  std::vector<int> geoIds;
  clusters.reserve(sourceLinks.size());
  for (auto sl : sourceLinks) {
    int geoId = sl.get().get<SimpleSourceLink>().geometryId().sensitive();
    clusters[geoId].push_back(sl);
  }
  if (clusters.size() < m_cfg.minSeedSize) {
    return;
  }
  geoIds.reserve(clusters.size());
  for (const auto &[geoId, sls] : clusters) {
    geoIds.push_back(geoId);
  }

#pragma omp parallel for num_threads(nThreads)
  for (std::size_t i = 0; i < geoIds.size() - 1; i++) {
    int fId = geoIds.at(i);
    const auto &fClusters = clusters.at(fId);
    for (std::size_t j = i + 1; j < geoIds.size(); j++) {
      int sId = geoIds.at(j);
      const auto &sClusters = clusters.at(sId);

      for (std::size_t n = 0; n < fClusters.size(); n++) {
        Acts::Vector3 flPoint =
            fClusters[n].get().get<SimpleSourceLink>().parametersGlob() - shift;
        double flX = flPoint.x();
        double flY = flPoint.y();
        double flZ = flPoint.z();
        for (std::size_t m = 0; m < sClusters.size(); m++) {
          Acts::Vector3 slPoint =
              sClusters[m].get().get<SimpleSourceLink>().parametersGlob() -
              shift;
          double diffX = slPoint.x() - flX;

          double thetaY = std::atan(diffX / (flY - slPoint.y()));
          double rhoY = flY * std::sin(thetaY) + flX * std::cos(thetaY);

          double thetaZ = std::atan(diffX / (flZ - slPoint.z()));
          double rhoZ = flZ * std::sin(thetaZ) + flX * std::cos(thetaZ);

          std::uint16_t cellThetaY =
              std::ceil((thetaY + M_PI_2) / m_deltaThetaY) - 1;
          std::uint16_t cellRhoY =
              std::ceil((rhoY + m_maxRhoY / 2.0) / m_deltaRhoY) - 1;
          std::uint16_t cellThetaZ =
              std::ceil((thetaZ + M_PI_2) / m_deltaThetaZ) - 1;
          std::uint16_t cellRhoZ =
              std::ceil((rhoZ + m_maxRhoZ / 2.0) / m_deltaRhoZ) - 1;

#pragma omp critical
          {
            votingMap[{cellThetaY, cellRhoY, cellThetaZ, cellRhoZ}]++;
          }
        }
      }
    }
  }
}
