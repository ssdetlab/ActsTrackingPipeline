#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "TrackingPipeline/Utilities/TupleHash.hpp"

class HoughTransformSeeder {
 public:
  using Cell =
      std::tuple<std::uint16_t, std::uint16_t, std::uint16_t, std::uint16_t>;
  using VotingMap = std::unordered_map<Cell, std::uint8_t, detail::TupleHash>;

  struct Config {
    std::size_t primaryIdx;
    std::size_t longIdx;
    std::size_t shortIdx;
  };

  struct Options {
    double boundBoxCenterPrimary;
    double boundBoxCenterLong;
    double boundBoxCenterShort;

    double boundBoxHalfPrimary;
    double boundBoxHalfLong;
    double boundBoxHalfShort;

    std::size_t nCellsThetaShort;
    std::size_t nCellsRhoShort;

    std::size_t nCellsThetaLong;
    std::size_t nCellsRhoLong;

    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface *>
        surfaceMap;

    int minXCount;
  };

  struct HTSeed {
    Acts::Vector3 lineRefPoint;
    Acts::Vector3 lineDir;
    std::vector<std::size_t> sourceLinkIdxs;
  };

  explicit HoughTransformSeeder(const Config &cfg);

  std::vector<HTSeed> findSeeds(
      const Acts::GeometryContext &gctx,
      const std::vector<Acts::SourceLink> &sourceLinks,
      const std::vector<std::size_t> &sourceLinksIndices, const Options &opt);

 private:
  Config m_cfg;

  void fillVotingMap(VotingMap &votingMap,
                     const std::vector<Acts::SourceLink> &sourceLinks,
                     const std::vector<std::size_t> &sourceLinksIndices,
                     const Options &opt, const Acts::Vector3 &shift,
                     double deltaThetaShort, double deltaRhoShort,
                     double deltaThetaLong, double deltaRhoLong,
                     double maxRhoLong, double maxRhoShort);

  std::vector<std::size_t> findLineSourceLinks(
      const std::vector<Acts::SourceLink> &sourceLinks,
      const std::vector<std::size_t> &sourceLinksIndices,
      const Acts::Vector3 &pointBL, const Acts::Vector3 &dirBL,
      const Acts::Vector3 &pointTL, const Acts::Vector3 &dirTL,
      const Acts::Vector3 &pointBR, const Acts::Vector3 &dirBR,
      const Acts::Vector3 &pointTR, const Acts::Vector3 &dirTR,
      const Acts::Vector3 &shift);

  std::unordered_map<Acts::GeometryIdentifier, double> m_geoPosMap;
};
