#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/BinningType.hpp"

#include "TrackingPipeline/Geometry/detail/BinningValueUtils.hpp"

struct SurfaceParameters {
  using AxisPars = std::tuple<Acts::BinningValue, double, double>;

  SurfaceParameters(AxisPars tPrimary, AxisPars tLong, AxisPars tShort, int id)
      : toWorldAnglePrimary(std::get<2>(tPrimary)),
        toWorldAngleLong(std::get<2>(tLong)),
        toWorldAngleShort(std::get<2>(tShort)),
        geoId(id) {
    toWorldTranslation[detail::binningValueToIndex(std::get<0>(tPrimary))] =
        std::get<1>(tPrimary);
    toWorldTranslation[detail::binningValueToIndex(std::get<0>(tLong))] =
        std::get<1>(tLong);
    toWorldTranslation[detail::binningValueToIndex(std::get<0>(tShort))] =
        std::get<1>(tShort);
  }
  SurfaceParameters() = delete;

  Acts::Vector3 toWorldTranslation;

  double toWorldAnglePrimary;
  double toWorldAngleLong;
  double toWorldAngleShort;

  int geoId;
};
