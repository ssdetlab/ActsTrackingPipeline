#include "TrackingPipeline/MagneticField/QuadrupoleMagField.hpp"

#include "Acts/Definitions/Units.hpp"

#include <iostream>

using namespace Acts::UnitLiterals;

QuadrupoleMagField::QuadrupoleMagField(double gradient)
    : m_gradient(gradient) {};

QuadrupoleMagField::QuadrupoleMagField(double gradient,
                                       const Acts::Vector3& origin,
                                       const Acts::RotationMatrix3& rotation,
                                       double length, double order)
    : m_gradient(gradient),
      m_origin(origin),
      m_rotation(rotation),
      m_order(order) {
  m_width = length / (2 * std::pow(std::log(2), 1.0 / m_order));
};

QuadrupoleMagField::~QuadrupoleMagField() = default;

Acts::MagneticFieldProvider::Cache QuadrupoleMagField::makeCache(
    const Acts::MagneticFieldContext& mctx) const {
  return Acts::MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
}

Acts::Result<Acts::Vector3> QuadrupoleMagField::getField(
    const Acts::Vector3& position, MagneticFieldProvider::Cache& cache) const {
  Acts::Vector3 global(position.x() - m_origin.x(), position.y() - m_origin.y(),
                       position.z() - m_origin.z());

  Acts::Vector3 local = m_rotation * global;
  Acts::Vector3 localB(m_gradient * local.y(), m_gradient * local.x(), 0);
  double factor = std::exp(-std::pow(std::abs(local.z() / m_width), m_order));
  localB *= factor;
  Acts::Vector3 globalB = m_rotation.inverse() * localB;

  return Acts::Result<Acts::Vector3>::success(globalB);
}
