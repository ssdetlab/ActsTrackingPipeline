#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <typeindex>
#include <unordered_map>

template <typename trajectory_t>
class MixedSourceLinkCalibrator {
 public:
  using Calibrator = Acts::Delegate<void(
      const Acts::GeometryContext&, const Acts::CalibrationContext&,
      const Acts::SourceLink&, typename trajectory_t::TrackStateProxy)>;

  MixedSourceLinkCalibrator() = default;
  ~MixedSourceLinkCalibrator() = default;

  template <auto Callable, typename T>
  void connect() {
    std::type_index typeIdx = typeid(T);
    m_calibratorMap.insert({typeIdx, Calibrator()});
    m_calibratorMap.at(typeIdx).template connect<Callable>();
  }

  void operator()(const Acts::GeometryContext& gctx,
                  const Acts::CalibrationContext& cctx,
                  const Acts::SourceLink& sourceLink,
                  typename trajectory_t::TrackStateProxy trackState) const {
    m_calibratorMap.at(sourceLink.type())(gctx, cctx, sourceLink, trackState);
  }

 private:
  std::unordered_map<std::type_index, Calibrator> m_calibratorMap;
};
