#include "TrackingPipeline/MagneticField/ConstantMagField.hpp"

#include <cstddef>

ConstantMagField::ConstantMagField(std::size_t id, const Acts::Vector3& field)
    : m_id(id), m_field(field) {}

Acts::Result<Acts::Vector3> ConstantMagField::getField(
    const Acts::Vector3& position, MagneticFieldProvider::Cache& cache) const {
  const auto& fieldCache = cache.as<Cache>();
  // std::cout << "DIPOLE " << m_id << " GET FIELD FIELD "
  //           << fieldCache.field.transpose() << "\n";
  return Acts::Result<Acts::Vector3>::success(fieldCache.field);
}

Acts::Result<Acts::Vector3> ConstantMagField::getFieldGradient(
    const Acts::Vector3& /*position*/, Acts::ActsMatrix<3, 3>& /*derivative*/,
    MagneticFieldProvider::Cache& /*cache*/) const {
  return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
}

Acts::MagneticFieldProvider::Cache ConstantMagField::makeCache(
    const Acts::MagneticFieldContext& mctx) const {
  // std::cout << "DIPOLE " << m_id << " MAKE CACHE\n";
  return Acts::MagneticFieldProvider::Cache(std::in_place_type<Cache>, m_id,
                                            mctx, m_field);
}
