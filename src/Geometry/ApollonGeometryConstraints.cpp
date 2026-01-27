#include "TrackingPipeline/Geometry/ApollonGeometryConstraints.hpp"

namespace ApollonGeometry {
// Static member definition (required by C++ standard if odr-used)
std::unique_ptr<const GeometryOptions> GeometryOptions::m_instance;
}  // namespace ApollonGeometry
