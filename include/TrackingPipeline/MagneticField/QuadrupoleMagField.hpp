#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include "G4MagneticField.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

class QuadrupoleMagField : public Acts::MagneticFieldProvider {
 public:
  /// @brief Cache for the magnetic field provider
  ///
  /// No specific cache is needed for the constant
  /// magnetic field
  struct Cache {
    /// @brief constructor with context
    Cache(const Acts::MagneticFieldContext& /*mcfg*/) {}
  };

  /// @brief Constructor with magnetic field gradient
  ///
  /// @param gradient magnetic field gradient
  QuadrupoleMagField(Acts::ActsScalar gradient);

  /// @brief Constructor with magnetic field gradient
  /// quadrupole origin and orientation
  ///
  /// @param gradient magnetic field gradient
  /// @param origin quadrupole origin
  /// @param rotation quadrupole orientation
  QuadrupoleMagField(Acts::ActsScalar gradient, const Acts::Vector3& origin,
                     const Acts::RotationMatrix3& rotation);

  ~QuadrupoleMagField() override;

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
      const Acts::Vector3& /*position*/, Acts::ActsMatrix<3, 3>& /*derivative*/,
      MagneticFieldProvider::Cache& /*cache*/) const override {
    Acts::Vector3 loaclGrad(m_gradient, m_gradient, 0);

    const Acts::Vector3 globalGrad = m_rotation.inverse() * loaclGrad;

    return Acts::Result<Acts::Vector3>::success(globalGrad);
  }

  /// @brief Get the magnetic field cache
  ///
  /// @param mctx Magnetic field context
  /// @return magnetic field cache
  Acts::MagneticFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override;

 private:
  Acts::ActsScalar m_gradient = 0.0;
  Acts::Vector3 m_origin = Acts::Vector3(0.0, 0.0, 0.0);
  Acts::RotationMatrix3 m_rotation = Acts::RotationMatrix3::Identity();
};
