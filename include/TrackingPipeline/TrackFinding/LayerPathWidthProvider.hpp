#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <map>

class LayerPathWidthProvider {
 public:
  LayerPathWidthProvider(std::map<int, std::pair<double, double>> widths)
      : m_widths(widths) {}

  std::map<int, std::pair<double, double>> m_widths;

  std::pair<double, double> operator()(
      const Acts::GeometryContext& /*gctx*/,
      const Acts::GeometryIdentifier& geoId) const {
    int staveId = static_cast<int>(geoId.sensitive() - 1);
    return m_widths.at(staveId);
  }
};
