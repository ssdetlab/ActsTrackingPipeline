#include "TrackingPipeline/MagneticField/CompositeMagField.hpp"

#include "Acts/Definitions/Algebra.hpp"

CompositeMagField::CompositeMagField(const FieldComponents& fieldComponents)
    : m_fieldComponents(fieldComponents) {};

CompositeMagField::~CompositeMagField() = default;

Acts::MagneticFieldProvider::Cache CompositeMagField::makeCache(
    const Acts::MagneticFieldContext& mctx) const {
  return Acts::MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
}

Acts::Result<Acts::Vector3> CompositeMagField::getField(
    const Acts::Vector3& position, MagneticFieldProvider::Cache& cache) const {
  Acts::Vector3 fieldValue = Acts::Vector3::Zero();
  for (const auto& [extent, field] : m_fieldComponents) {
    if (extent.contains(position)) {
      fieldValue = field->getField(position, cache).value();
      break;
    }
  }
  return Acts::Result<Acts::Vector3>::success(fieldValue);
}

Acts::Result<Acts::Vector3> CompositeMagField::getFieldGradient(
    const Acts::Vector3& position, Acts::ActsMatrix<3, 3>& derivative,
    MagneticFieldProvider::Cache& cache) const {
  (void)position;
  (void)derivative;
  (void)cache;
  return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
}
