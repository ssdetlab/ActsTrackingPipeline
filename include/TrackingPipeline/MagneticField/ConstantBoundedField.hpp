#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include <Acts/Utilities/BinningType.hpp>

/// @brief Constant magnetic field with bounded region
///
/// This class provides a constant magnetic field within a
/// bounded region. The magnetic field is zero outside the
/// bounded region.
class ConstantBoundedField : public Acts::MagneticFieldProvider {
 public:
  /// @brief Cache for the magnetic field provider
  ///
  /// No specific cache is needed for the constant
  /// magnetic field
  struct Cache {
    /// @brief constructor with context
    Cache(const Acts::MagneticFieldContext& /*mcfg*/) {}
  };

  /// @brief Constructor with magnetic field vector
  ///
  /// @param BField magnetic field vector
  /// @param BFieldExtent magnetic field extent
  explicit ConstantBoundedField(const Acts::Vector3& field,
                                const Acts::Extent& fieldExtent)
      : m_field(field), m_extent(fieldExtent) {}

  /// @brief Get the magnetic field at a given position
  ///
  /// @param position Vector3 position in global coordinate system
  /// @param cache Cache for the magnetic field provider
  /// @return magnetic field vector
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

  /// @brief Get the magnetic field gradient at a given position
  ///
  /// @param position Vector3 position in global coordinate system
  /// @param derivative ActsMatrix<3, 3> to store the gradient
  /// @param cache Cache for the magnetic field provider
  /// @return magnetic field gradient vector
  Acts::Result<Acts::Vector3> getFieldGradient(
      const Acts::Vector3& position, Acts::ActsMatrix<3, 3>& derivative,
      MagneticFieldProvider::Cache& cache) const override {
    (void)derivative;
    (void)cache;
    return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
  }

  /// @brief Get the magnetic field cache
  ///
  /// @param mctx Magnetic field context
  /// @return magnetic field cache
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
