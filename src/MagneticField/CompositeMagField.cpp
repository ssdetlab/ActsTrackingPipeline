#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"

#include "Acts/Definitions/Algebra.hpp"

CompositeMagField::CompositeMagField(std::size_t id,
                                     const FieldComponents& fieldComponents)
    : m_id(id), m_fieldComponents(fieldComponents) {};

CompositeMagField::~CompositeMagField() = default;

Acts::MagneticFieldProvider::Cache CompositeMagField::makeCache(
    const Acts::MagneticFieldContext& mctx) const {
  return Acts::MagneticFieldProvider::Cache(std::in_place_type<Cache>, m_id,
                                            m_fieldComponents, mctx);
}

Acts::Result<Acts::Vector3> CompositeMagField::getField(
    const Acts::Vector3& position, MagneticFieldProvider::Cache& cache) const {
  auto& componentCaches = cache.as<Cache>().m_componentCaches;
  auto& componentExtents = cache.as<Cache>().m_componentExtents;

  Acts::Vector3 fieldValue = Acts::Vector3::Zero();
  for (auto& [id, extent] : componentExtents) {
    if (extent.contains(position)) {
      fieldValue =
          m_fieldComponents.at(id)
              .fieldProvider->getField(position, componentCaches.at(id))
              .value();
      break;
    }
  }
  return Acts::Result<Acts::Vector3>::success(fieldValue);
}

Acts::Result<Acts::Vector3> CompositeMagField::getFieldGradient(
    const Acts::Vector3& /*position*/, Acts::ActsMatrix<3, 3>& /*derivative*/,
    MagneticFieldProvider::Cache& /*cache*/) const {
  return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
}
