#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <cstddef>
#include <map>
#include <span>
#include <vector>

class HoughTransformSeeder {
 public:
  using Cell = std::pair<int, int>;
  using VotingMap = std::map<int, std::map<Cell, int>>;
  using SourceLinkRef = std::reference_wrapper<const Acts::SourceLink>;

  struct Config {
    double boundBoxHalfX;
    double boundBoxHalfY;
    double boundBoxHalfZ;

    int nDirections;
    double phiMin;
    double phiMax;
    double thetaMin;
    double thetaMax;

    int nCellsX;
    int nCellsY;

    std::size_t minSeedSize;
    std::size_t maxSeedSize;

    std::size_t nLSIterations;
  };

  struct Options {
    double boundBoxCenterX;
    double boundBoxCenterY;
    double boundBoxCenterZ;
  };

  struct HTSeed {
    Acts::Vector3 lineRefPoint;
    Acts::Vector3 lineDir;
    std::vector<Acts::SourceLink> sourceLinks;
  };

  HoughTransformSeeder(const Config& cfg);

  std::vector<HTSeed> findSeeds(std::span<SourceLinkRef> sourceLinks,
                                const Options& opt);

 private:
  Config m_cfg;

  void fillVotingMap(VotingMap& votingMap, std::span<SourceLinkRef> points,
                     bool add);

  std::pair<Acts::Vector3, Acts::Vector3> findMaxVotedLine(
      const HoughTransformSeeder::VotingMap& votingMap);

  double computeDistance(const Acts::Vector3& point, const Acts::Vector3& dir,
                         const Acts::Vector3& meas);

  std::vector<int> findLineSourceLinks(std::span<SourceLinkRef> sourceLinks,
                                       const std::vector<bool>& activeIdxs,
                                       const Acts::Vector3& point,
                                       const Acts::Vector3& dir, double dist);

  double orthogonalLeastSquares(
      const std::vector<SourceLinkRef>& sourceLinks, Acts::Vector3& a,
      Acts::Vector3& b);

  Acts::Vector3 m_shift;

  std::vector<Acts::Vector3> m_dirs;

  double m_deltaX;
  double m_deltaY;
  double m_minDelta;
};
