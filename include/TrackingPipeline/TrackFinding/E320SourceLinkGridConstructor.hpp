#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Grid.hpp"

#include <unordered_map>

/// TODO: Fast version for single-chip-layer geometries?

/// @brief Class that constructs the source link grid
/// and sort event source links into it
///
/// Class takes the provied layer surfaces and constructs
/// source link grids based on their geometrical parameters.
/// The class is designed to work with layer-representing
/// surfaces combining the senstive surfaces of the detector
/// on the same binning direction, e.g. z-axis.
///
/// @note The convetion of double-digit geometry Ids with
/// the first digit representing the layer number is upheld
class E320SourceLinkGridConstructor {
 public:
  using AxisType =
      Acts::Axis<Acts::AxisType::Equidistant, Acts::AxisBoundaryType::Open>;
  using GridType =
      Acts::Grid<std::vector<Acts::SourceLink>, AxisType, AxisType>;

  /// @brief Nested configuration struct
  struct Config {
    /// Number of bins
    std::pair<int, int> bins;
    /// Merged layers to base the grid on
    std::unordered_map<Acts::GeometryIdentifier, const Acts::Surface*> layers;
    /// Surface accessor for source link transformation
    Acts::SourceLinkSurfaceAccessor surfaceAccessor;
  };

  /// @brief Constructor
  E320SourceLinkGridConstructor(const Config& config);

  /// @brief Bin the source links
  std::unordered_map<Acts::GeometryIdentifier, GridType> constructGrid(
      const Acts::GeometryContext& gctx,
      std::vector<Acts::SourceLink> sourceLinks);

  /// @brief Readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;
};
