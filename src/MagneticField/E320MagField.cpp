#include "TrackingPipeline/MagneticField/E320MagField.hpp"

#include "Acts/Definitions/Algebra.hpp"

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"

E320MagField::E320MagField(const double fieldS) : m_fieldS(fieldS) {};

E320MagField::~E320MagField() = default;

Acts::MagneticFieldProvider::Cache E320MagField::makeCache(
    const Acts::MagneticFieldContext& mctx) const {
  return Acts::MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
}

/// @brief helper function for exponential decay of dipole field
///
/// @param x position in space
/// @param params location and rate of decay
/// @return exponential damping factor
double E320MagField::decayFunction(const double x,
                                         const Acts::Vector4 params) const {
  return 1 / ((1 + exp((params[0] - x) / params[2])) *
              (1 + exp((x - params[1]) / params[3])));
}

Acts::Vector3 E320MagField::getDipole(
    const Acts::Vector3& pos,
    const std::tuple<Acts::Vector2, Acts::Vector4, Acts::Vector4>& dipoleParams)
    const {
  Acts::Vector3 fieldValue = Acts::Vector3::Zero();
  auto& [xParams, yParams, zParams] = dipoleParams;

  double a1 =
      (pos[0] > xParams[0] && pos[0] < xParams[1]) ? double(1.0) : double(0.0);

  double a2 = decayFunction(pos[1], yParams);

  double a3 = (pos[0] > -1000 && pos[0] < 1000) ? decayFunction(pos[1], zParams)
                                                : double(0.0);

  fieldValue[0] = m_fieldS * a1 * a2 * a3;
  return fieldValue;
}

Acts::Result<Acts::Vector3> E320MagField::getField(
    const Acts::Vector3& position, MagneticFieldProvider::Cache& cache) const {
  E320Geometry::GeometryOptions gOpt;
  double z = position.z();
  Acts::Vector3 fieldValue = Acts::Vector3::Zero();

  if (gOpt.quad1Translation.z() - gOpt.quad1Bounds[2] < z &&
      z <= gOpt.quad1Translation.z() + gOpt.quad1Bounds[2]) {
    auto grad1 = gOpt.quadrupolesParams[0];
    fieldValue = {grad1 * position[1], grad1 * position[0], 0};
    return Acts::Result<Acts::Vector3>::success(fieldValue);
  } else if (gOpt.quad2Translation.z() - gOpt.quad2Bounds[2] < z &&
             z <= gOpt.quad2Translation.z() + gOpt.quad2Bounds[2]) {
    auto grad2 = gOpt.quadrupolesParams[1];
    fieldValue = {grad2 * position[1], grad2 * position[0], 0};
    return Acts::Result<Acts::Vector3>::success(fieldValue);
  } else if (gOpt.quad3Translation.z() - gOpt.quad3Bounds[2] < z &&
             z <= gOpt.quad3Translation.z() + gOpt.quad3Bounds[2]) {
    auto grad3 = gOpt.quadrupolesParams[2];
    fieldValue = {grad3 * position[1], grad3 * position[0], 0};
    return Acts::Result<Acts::Vector3>::success(fieldValue);
  } else if (gOpt.dipoleTranslation.z() - gOpt.dipoleBounds[2] < z &&
             z <= gOpt.dipoleTranslation.z() + gOpt.dipoleBounds[2]) {
    fieldValue = getDipole(position, gOpt.dipoleParams);
    return Acts::Result<Acts::Vector3>::success(fieldValue);
  }
  return Acts::Result<Acts::Vector3>::success(fieldValue);
}

Acts::Result<Acts::Vector3> E320MagField::getFieldGradient(
    const Acts::Vector3& position, Acts::ActsMatrix<3, 3>& derivative,
    MagneticFieldProvider::Cache& cache) const {
  (void)position;
  (void)derivative;
  (void)cache;
  return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
}
