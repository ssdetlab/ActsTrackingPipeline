#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

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
  explicit ConstantBoundedField(const Acts::Vector3& BField,
                                const Acts::Extent& BFieldExtent)
      : m_BField(BField), m_Extent(BFieldExtent) {}

  /// @brief Get the magnetic field at a given position
  ///
  /// @param position Vector3 position in global coordinate system
  /// @param cache Cache for the magnetic field provider
  /// @return magnetic field vector
  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& position,
      MagneticFieldProvider::Cache& cache) const override {
    (void)cache;

    return (m_Extent.contains(position))
               ? Acts::Result<Acts::Vector3>::success(m_BField)
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
    return (m_Extent.contains(position))
               ? Acts::Result<Acts::Vector3>::success(m_BField)
               : Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
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
  /// magnetic field vector
  Acts::Vector3 m_BField;
  /// magnetic field extent
  Acts::Extent m_Extent;
};
