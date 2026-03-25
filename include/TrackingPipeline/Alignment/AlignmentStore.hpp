#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"

#include <unordered_map>

struct AlignmentStore {
  std::unordered_map<Acts::GeometryIdentifier, Acts::Transform3> store;
  Acts::ActsDynamicMatrix covariance;
};
