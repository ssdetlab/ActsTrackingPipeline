#pragma once

#include "Acts/Detector/Detector.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Delegate.hpp"

#include <typeindex>
#include <unordered_map>

class MixedSourceLinkSurfaceAccessor {
 public:
  using Accessor =
      Acts::Delegate<const Acts::Surface*(const Acts::SourceLink&)>;

  MixedSourceLinkSurfaceAccessor() = default;
  ~MixedSourceLinkSurfaceAccessor() = default;

  template <auto Callable, typename T, typename I>
  void connect(const I* instance) {
    std::type_index typeIdx = typeid(T);
    m_accessorMap.insert({typeIdx, Accessor()});
    m_accessorMap.at(typeIdx).template connect<Callable>(instance);
  }

  const Acts::Surface* operator()(const Acts::SourceLink& sourceLink) const {
    return m_accessorMap.at(sourceLink.type())(sourceLink);
  }

 private:
  std::unordered_map<std::type_index, Accessor> m_accessorMap;
};
