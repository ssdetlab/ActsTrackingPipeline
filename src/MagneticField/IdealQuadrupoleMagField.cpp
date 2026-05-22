#include "TrackingPipeline/MagneticField/IdealQuadrupoleMagField.hpp"

#include "Acts/Definitions/Algebra.hpp"

IdealQuadrupoleMagField::IdealQuadrupoleMagField(
    std::size_t id, double gradient, const Acts::Vector3& origin,
    const Acts::RotationMatrix3& rotation)
    : m_id(id), m_gradient(gradient), m_origin(origin), m_rotation(rotation) {};

IdealQuadrupoleMagField::~IdealQuadrupoleMagField() = default;

Acts::Result<Acts::Vector3> IdealQuadrupoleMagField::getField(
    const Acts::Vector3& position, MagneticFieldProvider::Cache& cache) const {
  const auto& quadCache = cache.as<Cache>();
  double gradient = quadCache.m_gradient;

  Acts::Vector3 global(position.x() - m_origin.x(), position.y() - m_origin.y(),
                       position.z() - m_origin.z());

  const Acts::Vector3 local = m_rotation * global;
  const Acts::Vector3 localB(gradient * local.y(), gradient * local.x(), 0);
  const Acts::Vector3 globalB = m_rotation.inverse() * localB;

  return Acts::Result<Acts::Vector3>::success(globalB);
}

Acts::Result<Acts::Vector3> IdealQuadrupoleMagField::getFieldGradient(
    const Acts::Vector3& /*position*/, Acts::ActsMatrix<3, 3>& /*derivative*/,
    MagneticFieldProvider::Cache& /*cache*/) const {
  return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
}

Acts::MagneticFieldProvider::Cache IdealQuadrupoleMagField::makeCache(
    const Acts::MagneticFieldContext& mctx) const {
  return Acts::MagneticFieldProvider::Cache(std::in_place_type<Cache>, m_id,
                                            mctx, m_gradient);
}
