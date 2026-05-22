#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <cstddef>

#include "TrackingPipeline/MagneticField/MagneticFieldStore.hpp"

class IdealQuadrupoleMagField : public Acts::MagneticFieldProvider {
 public:
  struct Cache {
    double m_gradient;

    explicit Cache(double grad) : m_gradient(grad) {}

    Cache(std::size_t id, const Acts::MagneticFieldContext& mctx,
          double defaultGrad) {
      if (!mctx.hasValue()) {
        m_gradient = defaultGrad;
        return;
      }
      const auto& store = mctx.get<std::shared_ptr<MagneticFieldStore>&>();
      if (store->store.contains(id)) {
        m_gradient = store->store.at(id).as<Cache>().m_gradient;
      } else {
        m_gradient = defaultGrad;
      }
    }
  };

  IdealQuadrupoleMagField(std::size_t id, double gradient,
                          const Acts::Vector3& origin,
                          const Acts::RotationMatrix3& rotation);

  ~IdealQuadrupoleMagField() override;

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
  double m_gradient;
  Acts::Vector3 m_origin;
  Acts::RotationMatrix3 m_rotation;
};
