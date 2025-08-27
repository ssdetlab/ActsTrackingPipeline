#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"

#include <cstddef>
#include <span>
#include <unordered_map>
#include <vector>

struct PairHash {
  std::size_t operator()(const std::pair<int, int>& p) const noexcept {
    return std::hash<long long>()(((long long)p.first << 32) ^
                                  (long long)p.second);
  }
};

class HoughTransformSeeder {
 public:
  using Cell = std::pair<int, int>;
  using GridMap = std::unordered_map<Cell, int, PairHash>;
  using VotingMap = std::unordered_map<int, GridMap>;
  using Index = std::tuple<int, int, int, int>;

  using SourceLinkRef = std::reference_wrapper<const Acts::SourceLink>;

  struct Config {
    double boundBoxHalfX;
    double boundBoxHalfY;
    double boundBoxHalfZ;

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

    int firstLayerId;
    int lastLayerId;
    int nLayers;
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

  void fillVotingMap(VotingMap& votingMap,
                     const std::vector<Acts::Vector3>& dirs,
                     std::span<SourceLinkRef> points, int sign, int nThreads);

  void fillVotingMap(VotingMap& votingMap,
                     const std::vector<Acts::Vector3>& dirs,
                     std::span<SourceLinkRef> points, int sign);

  std::pair<Acts::Vector3, Acts::Vector3> findMaxVotedLine(
      const HoughTransformSeeder::VotingMap& votingMap,
                     const std::vector<Acts::Vector3>& dirs
  );

  double computeDistance(const Acts::Vector3& point, const Acts::Vector3& dir,
                         const Acts::Vector3& meas);

  std::vector<int> findLineSourceLinks(std::span<SourceLinkRef> sourceLinks,
                                       const std::vector<bool>& activeIdxs,
                                       const Acts::Vector3& point,
                                       const Acts::Vector3& dir, double dist);

  double orthogonalLeastSquares(const std::vector<SourceLinkRef>& sourceLinks,
                                Acts::Vector3& a, Acts::Vector3& b);

  Acts::Vector3 m_shift;

  double m_deltaX;
  double m_deltaY;
  double m_minDelta;
};
