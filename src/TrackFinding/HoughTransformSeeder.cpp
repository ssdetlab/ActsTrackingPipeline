#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

#include <Acts/Definitions/Algebra.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <tuple>
#include <vector>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

double HoughTransformSeeder::orthogonalLeastSquares(
    const std::vector<SourceLinkRef>& sourceLinks, Acts::Vector3& a,
    Acts::Vector3& b) {
  double rc = 0.0;
  Acts::Vector3 meanVector(0, 0, 0);

  Eigen::MatrixXf points = Eigen::MatrixXf::Constant(sourceLinks.size(), 3, 0);
  for (int i = 0; i < points.rows(); i++) {
    Acts::Vector3 parameters =
        sourceLinks.at(i).get().get<SimpleSourceLink>().parametersGlob();
    points(i, 0) = parameters.x();
    points(i, 1) = parameters.y();
    points(i, 2) = parameters.z();

    meanVector += parameters - m_shift;
  }
  meanVector /= sourceLinks.size();
  a = meanVector;

  Eigen::MatrixXf centered = points.rowwise() - points.colwise().mean();
  Eigen::MatrixXf scatter = (centered.adjoint() * centered);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eig(scatter);
  Eigen::MatrixXf eigvecs = eig.eigenvectors();

  b[0] = eigvecs(0, 2);
  b[1] = eigvecs(1, 2);
  b[2] = eigvecs(2, 2);
  rc = eig.eigenvalues()(2);

  return rc;
}

double HoughTransformSeeder::computeDistance(const Acts::Vector3& point,
                                             const Acts::Vector3& dir,
                                             const Acts::Vector3& meas) {
  return (point + dir.dot(meas - point) * dir - meas).norm();
}

std::pair<Acts::Vector3, Acts::Vector3> HoughTransformSeeder::findMaxVotedLine(
    const HoughTransformSeeder::VotingMap& votingMap) {
  int maxDirIdx = 0;
  int maxCellXIdx = 0;
  int maxCellYIdx = 0;
  int maxCount = 0;
  for (auto [dirIdx, cells] : votingMap) {
    for (auto [cell, count] : cells) {
      if (count > maxCount) {
        maxCount = count;
        maxDirIdx = dirIdx;
        std::tie(maxCellXIdx, maxCellYIdx) = cell;
      }
    }
  }

  Acts::Vector3 dir = m_dirs.at(maxDirIdx);
  double denom = (1 + dir.z());
  Acts::Vector3 point =
      m_deltaX * maxCellXIdx *
          Acts::Vector3(1 - dir.x() * dir.x() / denom,
                        -dir.x() * dir.y() / denom, -dir.x()) +
      m_deltaY * maxCellYIdx *
          Acts::Vector3(-dir.x() * dir.y() / denom,
                        1 - dir.y() * dir.y() / denom, -dir.y());

  return {point, dir};
}

std::vector<int> HoughTransformSeeder::findLineSourceLinks(
    std::span<SourceLinkRef> sourceLinks, const std::vector<bool>& activeIdxs,
    const Acts::Vector3& point, const Acts::Vector3& dir, double dist) {
  std::vector<int> seedSlIdxs;
  for (int idx = 0; idx < activeIdxs.size(); idx++) {
    if (!activeIdxs.at(idx)) {
      continue;
    }

    sourceLinks[idx].get();
    if (computeDistance(
            point, dir,
            sourceLinks[idx].get().get<SimpleSourceLink>().parametersGlob() -
                m_shift) < dist) {
      seedSlIdxs.push_back(idx);
    }
  }
  return seedSlIdxs;
}

HoughTransformSeeder::HoughTransformSeeder(const Config& cfg) : m_cfg(cfg) {
  double r = (1 + std::sqrt(5)) / 2;
  m_dirs.reserve(m_cfg.nDirections);
  for (int i = 0; i < m_cfg.nDirections; i++) {
    double phi = std::fmod(2 * M_PI * i / r, 2 * M_PI) - M_PI;
    double theta = std::acos(1 - 2 * (i + 0.5) / m_cfg.nDirections);

    if (phi < m_cfg.phiMin || phi > m_cfg.phiMax || theta < m_cfg.thetaMin ||
        theta > m_cfg.thetaMax) {
      continue;
    }
    m_dirs.emplace_back(std::cos(phi) * std::sin(theta),
                        std::sin(phi) * std::sin(theta), std::cos(theta));
  }
  m_dirs.shrink_to_fit();

  double diag = std::sqrt(m_cfg.boundBoxHalfX * m_cfg.boundBoxHalfX +
                          m_cfg.boundBoxHalfY * m_cfg.boundBoxHalfY +
                          m_cfg.boundBoxHalfZ * m_cfg.boundBoxHalfZ);
  m_deltaX = diag / m_cfg.nCellsX;
  m_deltaY = diag / m_cfg.nCellsY;
  m_minDelta = std::min(m_deltaX, m_deltaY);
}

std::vector<HoughTransformSeeder::HTSeed> HoughTransformSeeder::findSeeds(
    std::span<SourceLinkRef> sourceLinks, const Options& opt) {
  std::vector<HoughTransformSeeder::HTSeed> seeds;

  m_shift = Acts::Vector3(opt.boundBoxCenterX, opt.boundBoxCenterY,
                          opt.boundBoxCenterZ);

  // std::cout << "SHIFT " << m_shift.transpose() << "\n";

  std::vector<bool> activeIdxs;
  activeIdxs.reserve(sourceLinks.size());
  int nActive = 0;
  for (std::size_t i = 0; i < sourceLinks.size(); i++) {
    activeIdxs.push_back(true);
  }
  nActive = activeIdxs.size();
  // std::cout << "N ACTIVE " << nActive << "\n";

  VotingMap votingMap;
  fillVotingMap(votingMap, sourceLinks, true);
  // std::cout << "FILLED VOTING MAP SIZE " << votingMap.size() << "\n";

  while (nActive >= m_cfg.minSeedSize) {
    auto [point, dir] = findMaxVotedLine(votingMap);
    // std::cout << "MAX VOTED POINT " << point.transpose() << "\n";
    // std::cout << "MAX VOTED DIR " << dir.transpose() << "\n";

    std::vector<int> seedSlIdxs;
    std::vector<SourceLinkRef> seedSourceLinksRefs;
    for (std::size_t it = 0; it < m_cfg.nLSIterations; it++) {
      seedSlIdxs = findLineSourceLinks(sourceLinks, activeIdxs, point, dir,
                                       3 * (m_cfg.nLSIterations - it) * m_minDelta);
      // std::cout << "LINE SOURCE LINKS " << seedSlIdxs.size() << "\n";

      seedSourceLinksRefs.reserve(seedSlIdxs.size());
      for (auto idx : seedSlIdxs) {
        seedSourceLinksRefs.push_back(std::cref(sourceLinks[idx]));
      }
      double rc = orthogonalLeastSquares(seedSourceLinksRefs, point, dir);
      if (rc == 0) {
        break;
      }
    }
    for (auto idx : seedSlIdxs) {
      activeIdxs.at(idx) = false;
      nActive--;
    }
    if (seedSlIdxs.size() < m_cfg.minSeedSize ||
        seedSlIdxs.size() > m_cfg.maxSeedSize) {
      break;
    }
    // std::cout << "ACTIVE LEFT " << nActive << "\n";

    std::vector<Acts::SourceLink> seedSourceLinks;
    seedSourceLinks.reserve(seedSlIdxs.size());
    for (auto idx : seedSlIdxs) {
      seedSourceLinks.push_back(sourceLinks[idx]);
    }

    fillVotingMap(votingMap, seedSourceLinksRefs, false);

    // std::cout << "SEED SIZE " << seedSourceLinks.size() << "\n";
    seeds.emplace_back(std::move(point + m_shift), std::move(dir),
                       std::move(seedSourceLinks));
  }

  return seeds;
}

void HoughTransformSeeder::fillVotingMap(VotingMap& votingMap,
                                         std::span<SourceLinkRef> sourceLinks,
                                         bool add) {
  for (auto sl : sourceLinks) {
    for (int i = 0; i < m_dirs.size(); i++) {
      const Acts::Vector3& dir = m_dirs.at(i);
      Acts::Vector3 point =
          sl.get().get<SimpleSourceLink>().parametersGlob() - m_shift;
      double denom = (1 + dir.z());
      double xPrime = (1 - dir.x() * dir.x() / denom) * point.x() -
                      (dir.x() * dir.y() / denom) * point.y() -
                      dir.x() * point.z();
      double yPrime = -(dir.x() * dir.y() / denom) * point.x() +
                      (1 - dir.y() * dir.y() / denom) * point.y() -
                      dir.y() * point.z();

      int cellX = std::ceil(xPrime / m_deltaX);
      int cellY = std::ceil(yPrime / m_deltaY);

      if (add) {
        votingMap[i][{cellX, cellY}]++;
      } else {
        votingMap[i][{cellX, cellY}]--;
      }
    }
  }
}
