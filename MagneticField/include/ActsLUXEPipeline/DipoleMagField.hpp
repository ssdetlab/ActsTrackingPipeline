#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Geometry/Extent.hpp"

/// @brief Constant magnetic field with bounded region
///
/// This class provides a constant magnetic field within a
/// bounded region. The magnetic field is zero outside the
/// bounded region.
class DipoleMagField : public Acts::MagneticFieldProvider {
public:
    /// @brief Cache for the magnetic field provider
    ///
    /// No specific cache is needed for the constant
    /// magnetic field
    struct Cache {
        /// @brief constructor with context
        Cache(const Acts::MagneticFieldContext& /*mcfg*/) {}
    };

    /// @brief Constructor with magnetic field vector
    ///
    /// @param params magnetic field parameters for Bx,By,Bz
    /// @param BFieldExtent magnetic field extent
    explicit DipoleMagField(
            const Acts::Vector2& xParams, const Acts::Vector4& yParams, const Acts::Vector4& zParams,
            const Acts::ActsScalar fieldS)
            : m_xParams(xParams),
              m_yParams(yParams),
              m_zParams(zParams),
              m_fieldS(fieldS){}

    const Acts::ActsScalar decayFunction(const Acts::ActsScalar x, const Acts::Vector4 params) const {
        return 1/ ( (1+exp((params[0]-x)/params[2])) * (1+exp((x-params[1])/params[3])));
    }
    /// @brief Get the magnetic field at a given position
    ///
    /// @param position Vector3 position in global coordinate system
    /// @param cache Cache for the magnetic field provider
    /// @return magnetic field vector
    Acts::Result<Acts::Vector3> getField(
            const Acts::Vector3& position,
            MagneticFieldProvider::Cache& cache) const override {
        (void)cache;
        Acts::Vector3 fieldValue = Acts::Vector3::Zero();
        Acts::ActsScalar a1 = (position[0] > m_xParams[0] && position[0] < m_xParams[1]) ?
                              Acts::ActsScalar(1.0) : Acts::ActsScalar(0.0);
        Acts::ActsScalar a2 = decayFunction(position[1],m_yParams);
        Acts::ActsScalar a3 = decayFunction(position[2],m_zParams);
        fieldValue[0] = a1 * a2 * a3;
        return Acts::Result<Acts::Vector3>::success(m_fieldS * fieldValue);
    }

    /// @brief Get the magnetic field gradient at a given position
    ///
    /// @param position Vector3 position in global coordinate system
    /// @param derivative ActsMatrix<3, 3> to store the gradient
    /// @param cache Cache for the magnetic field provider
    /// @return magnetic field gradient vector
    Acts::Result<Acts::Vector3> getFieldGradient(
            const Acts::Vector3& position, Acts::ActsMatrix<3, 3>& derivative,
            MagneticFieldProvider::Cache& cache) const override {
        (void)derivative;
        (void)cache;
        return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
    }

    /// @brief Get the magnetic field cache
    ///
    /// @param mctx Magnetic field context
    /// @return magnetic field cache
    Acts::MagneticFieldProvider::Cache makeCache(
            const Acts::MagneticFieldContext& mctx) const override {
        return Acts::MagneticFieldProvider::Cache(
                std::in_place_type<Cache>, mctx);
    }

private:
    /// magnetic field vector
    Acts::Vector2 m_xParams;
    Acts::Vector4 m_yParams;
    Acts::Vector4 m_zParams;
    Acts::ActsScalar m_fieldS;
};
