#include "TrackingPipeline/TrackFinding/HoughTransformSeeder.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include <Acts/Utilities/VectorHelpers.hpp>

#include <algorithm>
#include <cstddef>
#include <execution>
#include <functional>
#include <stdexcept>
#include <tuple>
#include <vector>

#include <omp.h>

#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"

void mapMerge(HoughTransformSeeder::VotingMap& out,
              HoughTransformSeeder::VotingMap& in) {
  out.merge(in);
}

#pragma omp declare reduction(                                                \
        mapAdd : HoughTransformSeeder::VotingMap : mapMerge(omp_out, omp_in)) \
    initializer(omp_priv = omp_orig)

double HoughTransformSeeder::orthogonalLeastSquares(
    const std::vector<SourceLinkRef>& sourceLinks, Acts::Vector3& a,
    Acts::Vector3& b) const {
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
                                             const Acts::Vector3& meas) const {
  return (point + dir.dot(meas - point) * dir - meas).norm();
}

std::pair<Acts::Vector3, Acts::Vector3> HoughTransformSeeder::findMaxVotedLine(
    const HoughTransformSeeder::VotingMap& votingMap,
    const std::vector<Acts::Vector3>& dirs) const {
  auto [maxCount, maxDirIdx, maxCellXIdx, maxCellYIdx] = std::transform_reduce(
      std::execution::par, votingMap.begin(), votingMap.end(),
      HoughTransformSeeder::Index{0, 0, 0, 0},
      [](const HoughTransformSeeder::Index& a,
         const HoughTransformSeeder::Index& b) {
        if (std::get<0>(a) != std::get<0>(b)) {
          return (std::get<0>(a) < std::get<0>(b)) ? b : a;
        }
        if (std::get<1>(a) != std::get<1>(b)) {
          return (std::get<1>(a) < std::get<1>(b)) ? b : a;
        }
        if (std::get<2>(a) != std::get<2>(b)) {
          return (std::get<2>(a) < std::get<2>(b)) ? b : a;
        }
        if (std::get<3>(a) != std::get<3>(b)) {
          return (std::get<3>(a) < std::get<3>(b)) ? b : a;
        }
        return a;
      },
      [](const auto& kv) {
        int localVMax = 0;
        int localFMax = 0;
        int localSMax = 0;
        for (const auto& [subKey, val] : kv.second) {
          if (val > localVMax) {
            localVMax = val;
            localFMax = subKey.first;
            localSMax = subKey.second;
          }
        }
        return std::make_tuple(localVMax, kv.first, localFMax, localSMax);
      });

  Acts::Vector3 dir = dirs.at(maxDirIdx);
  double denom = 1 + dir.z();
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
    const Acts::Vector3& point, const Acts::Vector3& dir, double dist) const {
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
  std::cout << "FIND SEEDS CALL\n";
  std::cout << "SOURCE LINKS " << sourceLinks.size() << "\n";
  std::cout << "LAYERS " << opt.nLayers << "\n";

  std::vector<SourceLinkRef> flPoints;
  std::vector<SourceLinkRef> llPoints;
  flPoints.reserve(sourceLinks.size() / opt.nLayers);
  llPoints.reserve(sourceLinks.size() / opt.nLayers);
  std::cout << "ITERATING SOURCE LINKS\n";
  for (auto sl : sourceLinks) {
    int geoId = sl.get().get<SimpleSourceLink>().geometryId().sensitive();
    if (geoId == opt.firstLayerId) {
      flPoints.push_back(sl);
    } else if (geoId == opt.lastLayerId) {
      llPoints.push_back(sl);
    }
  }
  if (flPoints.empty() || llPoints.empty()) {
    return {};
  }
  flPoints.shrink_to_fit();
  llPoints.shrink_to_fit();
  std::cout << "POINTS 1 " << flPoints.size() << "\n";
  std::cout << "POINTS 2 " << llPoints.size() << "\n";

  std::vector<Acts::Vector3> dirs;
  dirs.reserve(flPoints.size() * llPoints.size());
  for (auto fsl : flPoints) {
    for (auto lsl : llPoints) {
      Acts::Vector3 dir = (lsl.get().get<SimpleSourceLink>().parametersGlob() -
                           fsl.get().get<SimpleSourceLink>().parametersGlob())
                              .normalized();
      double theta = Acts::VectorHelpers::theta(dir);
      double phi = Acts::VectorHelpers::phi(dir);
      if (theta < opt.minTheta || theta > opt.maxTheta) {
        continue;
      }
      if (phi < opt.minPhi || phi > opt.maxPhi) {
        continue;
      }
      dirs.push_back(dir);
    }
  }
  dirs.shrink_to_fit();
  if (dirs.empty()) {
    return {};
  }
  std::cout << "DIRS " << dirs.size() << "\n";
  // throw std::runtime_error("err");

  m_shift = Acts::Vector3(opt.boundBoxCenterX, opt.boundBoxCenterY,
                          opt.boundBoxCenterZ);

  std::vector<bool> activeIdxs;
  activeIdxs.reserve(sourceLinks.size());
  int nActive = 0;
  for (std::size_t i = 0; i < sourceLinks.size(); i++) {
    activeIdxs.push_back(true);
  }
  nActive = activeIdxs.size();

  VotingMap votingMap;
  fillVotingMap(votingMap, dirs, sourceLinks, 1, 16);

  while (nActive >= m_cfg.minSeedSize) {
    auto [point, dir] = findMaxVotedLine(votingMap, dirs);

    std::vector<int> seedSlIdxs;
    std::vector<SourceLinkRef> seedSourceLinksRefs;
    double rc = 0;
    for (std::size_t it = 0; it < m_cfg.nLSIterations; it++) {
      seedSourceLinksRefs.clear();
      seedSlIdxs =
          findLineSourceLinks(sourceLinks, activeIdxs, point, dir, 5e-2);

      seedSourceLinksRefs.reserve(seedSlIdxs.size());
      for (auto idx : seedSlIdxs) {
        seedSourceLinksRefs.push_back(std::cref(sourceLinks[idx]));
      }
      rc = orthogonalLeastSquares(seedSourceLinksRefs, point, dir);
      if (rc == 0) {
        break;
      }
    }
    std::size_t slIdxsSize = seedSlIdxs.size();
    for (auto idx : seedSlIdxs) {
      activeIdxs.at(idx) = false;
      nActive--;
    }
    if (slIdxsSize < m_cfg.minSeedSize || slIdxsSize > m_cfg.maxSeedSize) {
      break;
    }

    std::vector<Acts::SourceLink> seedSourceLinks;
    seedSourceLinks.reserve(slIdxsSize);
    for (auto idx : seedSlIdxs) {
      seedSourceLinks.push_back(sourceLinks[idx]);
    }

    fillVotingMap(votingMap, dirs, seedSourceLinksRefs, -1);

    seeds.emplace_back(std::move(point + m_shift), std::move(dir),
                       std::move(seedSourceLinks));
  }

  return seeds;
}

void HoughTransformSeeder::fillVotingMap(VotingMap& votingMap,
                                         const std::vector<Acts::Vector3>& dirs,
                                         std::span<SourceLinkRef> sourceLinks,
                                         int sign) const {
  for (int i = 0; i < dirs.size(); i++) {
    const Acts::Vector3& dir = dirs.at(i);
    double dirX = dir.x();
    double dirY = dir.y();
    double dirZ = dir.z();

    double denom = 1 + dirZ;
    for (auto sl : sourceLinks) {
      Acts::Vector3 point =
          sl.get().get<SimpleSourceLink>().parametersGlob() - m_shift;
      double pointX = point.x();
      double pointY = point.y();
      double pointZ = point.z();

      double xPrime = (1 - dirX * dirX / denom) * pointX -
                      (dirX * dirY / denom) * pointY - dirX * pointZ;
      double yPrime = -(dirX * dirY / denom) * pointX +
                      (1 - dirY * dirY / denom) * pointY - dirY * pointZ;

      int cellX = std::ceil(xPrime / m_deltaX);
      int cellY = std::ceil(yPrime / m_deltaY);

      votingMap[i][{cellX, cellY}] += sign;
    }
  }
}

void HoughTransformSeeder::fillVotingMap(VotingMap& votingMap,
                                         const std::vector<Acts::Vector3>& dirs,
                                         std::span<SourceLinkRef> sourceLinks,
                                         int sign, int nThreads) const {
#pragma omp parallel for num_threads(nThreads) reduction(mapAdd : votingMap)
  for (int i = 0; i < dirs.size(); i++) {
    const Acts::Vector3& dir = dirs.at(i);
    double dirX = dir.x();
    double dirY = dir.y();
    double dirZ = dir.z();

    double denom = 1 + dirZ;
    for (const auto sl : sourceLinks) {
      Acts::Vector3 point =
          sl.get().get<SimpleSourceLink>().parametersGlob() - m_shift;
      double pointX = point.x();
      double pointY = point.y();
      double pointZ = point.z();

      double xPrime = (1 - dirX * dirX / denom) * pointX -
                      (dirX * dirY / denom) * pointY - dirX * pointZ;
      double yPrime = -(dirX * dirY / denom) * pointX +
                      (1 - dirY * dirY / denom) * pointY - dirY * pointZ;

      int cellX = std::ceil(xPrime / m_deltaX);
      int cellY = std::ceil(yPrime / m_deltaY);

      votingMap[i][{cellX, cellY}] += sign;
    }
  }
}
