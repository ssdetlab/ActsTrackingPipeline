#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/BinningType.hpp"

class ConstantBoundedMagField : public Acts::MagneticFieldProvider {
 public:
  struct Cache {
    Cache(const Acts::MagneticFieldContext& /*mctx*/) {}
  };

  explicit ConstantBoundedMagField(const Acts::Vector3& field,
                                   const Acts::Extent& fieldExtent)
      : m_field(field), m_extent(fieldExtent) {}

  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& position,
      MagneticFieldProvider::Cache& cache) const override {
    (void)cache;

    bool containsX = (position.x() > m_extent.min(Acts::BinningValue::binX)) &&
                     (position.x() < m_extent.max(Acts::BinningValue::binX));
    bool containsY = (position.y() > m_extent.min(Acts::BinningValue::binY)) &&
                     (position.y() < m_extent.max(Acts::BinningValue::binY));
    bool containsZ = (position.z() > m_extent.min(Acts::BinningValue::binZ)) &&
                     (position.z() < m_extent.max(Acts::BinningValue::binZ));

    return (containsX && containsY && containsZ)
               ? Acts::Result<Acts::Vector3>::success(m_field)
               : Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
  }

  Acts::Result<Acts::Vector3> getFieldGradient(
      const Acts::Vector3& position, Acts::ActsMatrix<3, 3>& derivative,
      MagneticFieldProvider::Cache& cache) const override {
    (void)derivative;
    (void)cache;
    return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
  }

  Acts::MagneticFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override {
    return Acts::MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

 private:
  /// Uagnetic field vector
  Acts::Vector3 m_field;

  /// Magnetic field extent
  Acts::Extent m_extent;
};
