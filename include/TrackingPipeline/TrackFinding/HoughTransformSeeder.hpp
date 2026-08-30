#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <cstddef>
#include <unordered_map>
#include <vector>

#include "TrackingPipeline/Utilities/TupleHash.hpp"

/// @brief Hough-Transform-based seeder
///
/// A class performing Hough Transform-based seeding. Initial guess for the
/// track poisition and direction is made based on the Hough space line
/// parameters.
class HoughTransformSeeder {
 public:
  /// Short hands for the Hough space cell and voting map
  using Cell =
      std::tuple<std::uint16_t, std::uint16_t, std::uint16_t, std::uint16_t>;
  using VotingMap = std::unordered_map<Cell, std::uint8_t, detail::TupleHash>;

  /// @brief Configuration struct
  struct Config {
    /// Space indices
    std::size_t primaryIdx;
    std::size_t longIdx;
    std::size_t shortIdx;
  };

  /// @brief Options struct
  struct Options {
    /// Position of the real space bounding box
    double boundBoxCenterPrimary;
    double boundBoxCenterLong;
    double boundBoxCenterShort;
    /// Dimensions of the real space bounding box
    double boundBoxHalfPrimary;
    double boundBoxHalfLong;
    double boundBoxHalfShort;
    /// Hough space binning in the long coordinate
    std::size_t nCellsThetaShort;
    std::size_t nCellsRhoShort;
    /// Hough space binning in the short coordinate
    std::size_t nCellsThetaLong;
    std::size_t nCellsRhoLong;
    /// Surface map
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface *>
        surfaceMap;
    /// Minimum number of Hough space intersections cut
    int minXCount;
  };

  /// @brief Hough Transform seed struct
  struct HTSeed {
    /// Passing point of the seed
    Acts::Vector3 lineRefPoint;
    /// Direction of the seed
    Acts::Vector3 lineDir;
    /// Seed source link indices
    std::vector<std::size_t> sourceLinkIdxs;
    /// Number of cell intersection counts
    std::size_t xCount;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  explicit HoughTransformSeeder(const Config &cfg);

  /// @brief perform seed finding
  ///
  /// @param gctx geometry context
  /// @param sourceLinks source link pull to find the seeds in
  /// @param sourceLinksIndices collection of source link indices
  /// @param opt seed finding options
  ///
  /// @return collection of Hough Transform seeds
  std::vector<HTSeed> findSeeds(
      const Acts::GeometryContext &gctx,
      const std::vector<Acts::SourceLink> &sourceLinks,
      const std::vector<std::size_t> &sourceLinksIndices, const Options &opt);

 private:
  /// @brief fill the voting map with Hough Space intersections
  ///
  /// @param votingMap voting map to be filled
  /// @sourceLinks source links collection
  /// @sourceLinksIndices source links indices collection
  /// @opt seed finding options
  /// @shift real space bounding box center
  /// @deltaThetaShort short axis theta bin size
  /// @deltaRhoShort short axis rho bin size
  /// @deltaThetaLong long axis theta bin size
  /// @deltaRhoLong long axis rho bin size
  /// @maxRhoShort short axis rho cap
  /// @maxRhoLong long axis rho cap
  void fillVotingMap(VotingMap &votingMap,
                     const std::vector<Acts::SourceLink> &sourceLinks,
                     const std::vector<std::size_t> &sourceLinksIndices,
                     const Options &opt, const Acts::Vector3 &shift,
                     double deltaThetaShort, double deltaRhoShort,
                     double deltaThetaLong, double deltaRhoLong,
                     double maxRhoLong, double maxRhoShort);

  /// @brief find source links falling within the real space tunnel
  ///
  /// @sourceLinks source links collection
  /// @sourceLinksIndices source links indices collection
  /// @pointBL real space bottom left tunnel line point
  /// @dirBL real space bottom left tunnel line direction
  /// @pointTL real space top left tunnel line point
  /// @dirTL real space top left tunnel line direction
  /// @pointBR real space bottom right tunnel line point
  /// @dirBR real space bottom right tunnel line direction
  /// @pointTR real space top right tunnel line point
  /// @dirTR real space top right tunnel line direction
  /// @shift real space bounding box center
  ///
  /// @return indices of the source links within the real space tunnel
  std::vector<std::size_t> findLineSourceLinks(
      const std::vector<Acts::SourceLink> &sourceLinks,
      const std::vector<std::size_t> &sourceLinksIndices,
      const Acts::Vector3 &pointBL, const Acts::Vector3 &dirBL,
      const Acts::Vector3 &pointTL, const Acts::Vector3 &dirTL,
      const Acts::Vector3 &pointBR, const Acts::Vector3 &dirBR,
      const Acts::Vector3 &pointTR, const Acts::Vector3 &dirTR,
      const Acts::Vector3 &shift);

  /// Configuration
  Config m_cfg;

  /// Map of primary axis surfaces positions in the real space bounding box
  std::unordered_map<Acts::GeometryIdentifier, double> m_geoPosMap;
};
