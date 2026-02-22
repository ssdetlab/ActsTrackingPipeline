#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"

namespace E320Geometry {
// Static member definition (required by C++ standard if odr-used)
std::unique_ptr<const GeometryOptions> GeometryOptions::m_instance;
}  // namespace E320Geometry