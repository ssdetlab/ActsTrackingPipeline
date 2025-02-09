#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <map>

class LayerPathWidthProvider {
 public:
  LayerPathWidthProvider(
      std::map<std::int32_t, std::pair<double, double>> widths)
      : m_widths(widths) {}

  std::map<std::int32_t, std::pair<double, double>> m_widths;

  std::pair<double, double> operator()(
      const Acts::GeometryContext& /*gctx*/,
      const Acts::GeometryIdentifier& geoId) const {
    std::int32_t staveId = static_cast<std::int32_t>(geoId.sensitive() - 1);
    return m_widths.at(staveId);
  }
};
