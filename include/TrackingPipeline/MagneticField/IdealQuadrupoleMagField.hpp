#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <cstddef>

#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

/// @brief Hard-edge quadupole magnetic field provider
class IdealQuadrupoleMagField : public Acts::MagneticFieldProvider {
 public:
  /// @brief Magnetic field cache struct
  struct Cache {
    /// Field gradient
    double m_gradient;

    /// @brief Constructor
    ///
    /// @param grad field gradient
    explicit Cache(double grad) : m_gradient(grad) {}

    /// @brief Constructor
    ///
    /// @param id field volume ID
    /// @mctx magnetic field context
    /// @defaultGrad default field gradient value
    Cache(std::size_t id, const Acts::MagneticFieldContext& mctx,
          double defaultGrad) {
      if (!mctx.hasValue()) {
        m_gradient = defaultGrad;
        return;
      }
      const auto& store =
          mctx.get<std::shared_ptr<MagneticFieldStore>&>()->store;
      if (store.contains(id)) {
        m_gradient = store.at(id).as<Cache>().m_gradient;
      } else {
        m_gradient = defaultGrad;
      }
    }
  };

  /// @brief Constructor
  ///
  /// @param id field volume ID
  /// @param gradient default field gradient
  /// @param origin quad center in global coordinates
  /// @param rotation quad rotation in global coordinates
  IdealQuadrupoleMagField(std::size_t id, double gradient,
                          const Acts::Vector3& origin,
                          const Acts::RotationMatrix3& rotation);

  /// @brief Get field value in a point
  ///
  /// @param position position of the field inquiry
  /// @param cache current field cache
  ///
  /// @return field value in the result wrapper
  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& position,
      MagneticFieldProvider::Cache& cache) const override;

  /// @brief Get field gradient in a point
  ///
  /// @param position position of the field inquiry
  /// @param derivative field jacobian in global coordinates
  /// @param cache current field cache
  ///
  /// @return field gradient in the result wrapper
  Acts::Result<Acts::Vector3> getFieldGradient(
      const Acts::Vector3& position, Acts::ActsMatrix<3, 3>& derivative,
      MagneticFieldProvider::Cache& cache) const override;

  /// @brief Initialize field cache
  ///
  /// @param mctx current magnetic field context
  ///
  /// @return field cache instance
  Acts::MagneticFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override;

 private:
  /// Magnetic volume ID
  std::size_t m_id;

  /// Field gradient
  double m_gradient;

  /// Quad center in global coordinates
  Acts::Vector3 m_origin;

  /// Quad rotation in global coordinates
  Acts::RotationMatrix3 m_rotation;
};
