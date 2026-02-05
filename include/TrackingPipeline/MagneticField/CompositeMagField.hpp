#pragma once

#include "Acts/Geometry/Extent.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <memory>

/// @brief Constant magnetic field with bounded region
///
/// This class provides a constant magnetic field within a
/// bounded region. The magnetic field is zero outside the
/// bounded region.
class CompositeMagField : public Acts::MagneticFieldProvider {
 public:
  using FieldComponent =
      std::pair<Acts::Extent,
                const std::shared_ptr<Acts::MagneticFieldProvider>>;
  using FieldComponents = std::vector<FieldComponent>;

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
  /// @param params magnetic field parameters for Bx,By,Bz
  /// @param BFieldExtent magnetic field extent
  CompositeMagField(const FieldComponents& fieldComponents);

  ~CompositeMagField() override;

  /// @brief Get the magnetic field cache
  ///
  /// @param mctx Magnetic field context
  /// @return magnetic field cache
  Acts::MagneticFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override;

  /// @brief Get the magnetic field at a given position
  ///
  /// @param position Vector3 position in global coordinate system
  /// @param cache Cache for the magnetic field provider
  /// @return magnetic field vector
  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& position,
      MagneticFieldProvider::Cache& cache) const override;

  /// @brief Get the magnetic field gradient at a given position
  ///
  /// @param position Vector3 position in global coordinate system
  /// @param derivative ActsMatrix<3, 3> to store the gradient
  /// @param cache Cache for the magnetic field provider
  /// @return magnetic field gradient vector
  Acts::Result<Acts::Vector3> getFieldGradient(
      const Acts::Vector3& position, Acts::ActsMatrix<3, 3>& derivative,
      MagneticFieldProvider::Cache& cache) const override;

 private:
  FieldComponents m_fieldComponents;
};
