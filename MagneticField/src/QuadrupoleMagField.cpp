#include "ActsLUXEPipeline/QuadrupoleMagField.hpp"

QuadrupoleMagField::QuadrupoleMagField(Acts::ActsScalar gradient) 
    : m_gradient(gradient) {};

QuadrupoleMagField::~QuadrupoleMagField() = default;

Acts::MagneticFieldProvider::Cache QuadrupoleMagField::makeCache(
    const Acts::MagneticFieldContext& mctx) const {
        return Acts::MagneticFieldProvider::Cache(
            std::in_place_type<Cache>, mctx);
}

Acts::Result<Acts::Vector3> QuadrupoleMagField::getField(
    const Acts::Vector3& position, 
    MagneticFieldProvider::Cache& cache) const {
        const Acts::Vector3 fieldValue(
            m_gradient * position.y(),
            m_gradient * position.x(),
            0);

        return Acts::Result<Acts::Vector3>::success(fieldValue);
}

