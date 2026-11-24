#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "TrackingPipeline/Utilities/TupleHash.hpp"

/// @brief Hough-Transform-based seeder
///
/// A class performing Hough Transform-based seeding. The real space lines are
/// parametrized by their slopes and intercepts in the primary-long and the
/// primary-short planes. The intercepts are taken at the primary position of
/// the real space bounding box center. A measurement maps to a straight line
/// in the slope-intercept Hough space, and every pair of measurements on
/// different primary positions votes for the Hough space cell containing the
/// intersection of their Hough space lines. Initial guess for the track
/// poisition and direction is made based on the Hough space cell center.
class HoughTransformSeeder {
 public:
  /// Short hands for the Hough space cell and voting map
  ///
  /// The cell indices are (slope long, intercept long, slope short,
  /// intercept short)
  using Cell =
      std::tuple<std::uint16_t, std::uint16_t, std::uint16_t, std::uint16_t>;
  using VotingMap = std::unordered_map<Cell, std::size_t, detail::TupleHash>;

  /// @brief Configuration struct
  struct Config {
    /// Space indices
    std::size_t primaryIdx;
    std::size_t longIdx;
    std::size_t shortIdx;
  };

  /// @brief Options struct
  ///
  /// The Hough space ranges are half-open intervals [min, max), the pairs
  /// voting outside of the ranges are discarded.
  struct Options {
    /// Position of the real space bounding box, the slopes and intercepts are
    /// defined relative to it
    double boundBoxCenterPrimary;
    double boundBoxCenterLong;
    double boundBoxCenterShort;
    /// Hough space slope range and binning in the long coordinate
    double minSlopeLong;
    double maxSlopeLong;
    std::size_t nCellsSlopeLong;
    /// Hough space intercept range and binning in the long coordinate
    double minInterceptLong;
    double maxInterceptLong;
    std::size_t nCellsInterceptLong;
    /// Hough space slope range and binning in the short coordinate
    double minSlopeShort;
    double maxSlopeShort;
    std::size_t nCellsSlopeShort;
    /// Hough space intercept range and binning in the short coordinate
    double minInterceptShort;
    double maxInterceptShort;
    std::size_t nCellsInterceptShort;
    /// Minimum number of Hough space intersections cut
    std::size_t minXCount;
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
  /// @param sourceLinks source links collection
  /// @param sourceLinksIndices source links indices collection
  /// @param opt seed finding options
  /// @param shift real space bounding box center
  /// @param deltaSlopeLong long axis slope bin size
  /// @param deltaInterceptLong long axis intercept bin size
  /// @param deltaSlopeShort short axis slope bin size
  /// @param deltaInterceptShort short axis intercept bin size
  void fillVotingMap(VotingMap &votingMap,
                     const std::vector<Acts::SourceLink> &sourceLinks,
                     const std::vector<std::size_t> &sourceLinksIndices,
                     const Options &opt, const Acts::Vector3 &shift,
                     double deltaSlopeLong, double deltaInterceptLong,
                     double deltaSlopeShort, double deltaInterceptShort);

  /// @brief find source links falling within the real space tunnel
  ///
  /// The tunnel is the set of the real space lines in a Hough space cell. A
  /// source link is within the tunnel if its Hough space line crosses the
  /// cell, i.e. if its transverse coordinates lie within the range spanned by
  /// the cell lines at its primary position.
  ///
  /// @param sourceLinks source links collection
  /// @param sourceLinksIndices source links indices collection
  /// @param shift real space bounding box center
  /// @param minSlopeLong cell lower long slope edge
  /// @param maxSlopeLong cell upper long slope edge
  /// @param minInterceptLong cell lower long intercept edge
  /// @param maxInterceptLong cell upper long intercept edge
  /// @param minSlopeShort cell lower short slope edge
  /// @param maxSlopeShort cell upper short slope edge
  /// @param minInterceptShort cell lower short intercept edge
  /// @param maxInterceptShort cell upper short intercept edge
  ///
  /// @return indices of the source links within the real space tunnel
  std::vector<std::size_t> findLineSourceLinks(
      const std::vector<Acts::SourceLink> &sourceLinks,
      const std::vector<std::size_t> &sourceLinksIndices,
      const Acts::Vector3 &shift, double minSlopeLong, double maxSlopeLong,
      double minInterceptLong, double maxInterceptLong, double minSlopeShort,
      double maxSlopeShort, double minInterceptShort, double maxInterceptShort);

  /// Configuration
  Config m_cfg;
};
