#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <cstddef>

#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

/// @brief Hard-edge uniform magnetic field provider
class ConstantMagField : public Acts::MagneticFieldProvider {
 public:
  /// @brief Magnetic field cache struct
  struct Cache {
    /// Field vector
    Acts::Vector3 m_field;

    /// @brief Constructor
    ///
    /// @param strength field strength in T
    /// @param dirIdx field direction index in global coordinates
    Cache(double strength, std::size_t dirIdx) {
      m_field = Acts::Vector3::Zero();
      m_field(dirIdx) = strength;
    }

    /// @brief Constructor
    ///
    /// @param id field volume ID
    /// @mctx magnetic field context
    /// @defaultField default field vector
    Cache(std::size_t id, const Acts::MagneticFieldContext& mctx,
          const Acts::Vector3& defaultField) {
      if (!mctx.hasValue()) {
        m_field = defaultField;
        return;
      }
      const auto& store = mctx.get<std::shared_ptr<MagneticFieldStore>&>();
      if (store->store.contains(id)) {
        m_field = store->store.at(id).as<Cache>().m_field;
      } else {
        m_field = defaultField;
      }
    }
  };

  /// @brief Constructor
  ///
  /// @param id field volume ID
  /// @param field default field vector
  ConstantMagField(std::size_t id, const Acts::Vector3& field);

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

  /// Field vector
  Acts::Vector3 m_field;
};
