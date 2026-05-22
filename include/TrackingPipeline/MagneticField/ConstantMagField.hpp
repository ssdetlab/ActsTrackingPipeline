#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <cstddef>

#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

class ConstantMagField : public Acts::MagneticFieldProvider {
 public:
  struct Cache {
    Acts::Vector3 m_field;
    Cache(double strength, std::size_t dirIdx) {
      m_field = Acts::Vector3::Zero();
      m_field(dirIdx) = strength;
    }

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

  ConstantMagField(std::size_t id, const Acts::Vector3& field);

  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& position,
      MagneticFieldProvider::Cache& cache) const override;

  Acts::Result<Acts::Vector3> getFieldGradient(
      const Acts::Vector3& position, Acts::ActsMatrix<3, 3>& derivative,
      MagneticFieldProvider::Cache& cache) const override;

  Acts::MagneticFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override;

 private:
  std::size_t m_id;
  Acts::Vector3 m_field;
};
