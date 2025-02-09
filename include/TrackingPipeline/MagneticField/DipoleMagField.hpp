#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"

class DipoleMagField : public Acts::MagneticFieldProvider {
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
  /// @param params magnetic field parameters for Bx,By,Bz
  /// @param BFieldExtent magnetic field extent
  explicit DipoleMagField(
      const std::tuple<Acts::Vector2, Acts::Vector4, Acts::Vector4> params,
      const double fieldStrength,
      const Acts::RotationMatrix3& rotation = Acts::RotationMatrix3::Identity(),
      const Acts::Vector3& origin = Acts::Vector3::Zero())
      : m_xParams(std::get<0>(params)),
        m_yParams(std::get<1>(params)),
        m_zParams(std::get<2>(params)),
        m_fieldStrength(fieldStrength),
        m_rotation(rotation),
        m_origin(origin) {}

  const double decayFunction(const double x, const Acts::Vector4 params) const {
    return 1 / ((1 + exp((params[0] - x) / params[2])) *
                (1 + exp((x - params[1]) / params[3])));
  }

  /// @brief Get the magnetic field at a given position
  ///
  /// @param position Vector3 position in global coordinate system
  /// @param cache Cache for the magnetic field provider
  /// @return magnetic field vector
  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& position,
      MagneticFieldProvider::Cache& cache) const override {
    (void)cache;
    Acts::Vector3 localB = Acts::Vector3::Zero();

    Acts::Vector3 local = position - m_origin;

    local = m_rotation * local;

    double a1 = (local.x() > m_xParams[0] && local.x() < m_xParams[1])
                    ? double(1.0)
                    : double(0.0);
    double a2 = decayFunction(local.y(), m_yParams);
    double a3 = decayFunction(local.z(), m_zParams);

    localB[0] = m_fieldStrength * a1 * a2 * a3;

    Acts::Vector3 globalB = m_rotation.inverse() * localB;

    return Acts::Result<Acts::Vector3>::success(globalB);
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
  Acts::Vector2 m_xParams;
  Acts::Vector4 m_yParams;
  Acts::Vector4 m_zParams;
  double m_fieldStrength;
  Acts::RotationMatrix3 m_rotation;
  Acts::Vector3 m_origin;
};
