#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"

/// @brief Constant magnetic field with bounded region
class E320MagField : public Acts::MagneticFieldProvider {
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
        E320MagField(const Acts::ActsScalar fieldS);
    
        ~E320MagField() override;
    
        /// @brief Get the magnetic field cache
        ///
        /// @param mctx Magnetic field context
        /// @return magnetic field cache
        Acts::MagneticFieldProvider::Cache makeCache(
                const Acts::MagneticFieldContext& mctx) const override;
    
        /// @brief helper function
        const Acts::ActsScalar decayFunction(
                const Acts::ActsScalar x, const Acts::Vector4 params) const;
    
        /// @brief calculate the dipole part of the magnetic field
        ///
        /// @param pos position in space
        /// @param dipoleParams tuple<xParams,yParams,zParams>
        /// @return value of dipole in the x-direction
        const Acts::Vector3 getDipole(
                const Acts::Vector3& pos,
                const std::tuple<Acts::Vector2,
                                Acts::Vector4,
                                Acts::Vector4>& dipoleParams) const;
    
        /// @brief Get the magnetic field at a given position
        ///
        /// @param position Vector3 position in global coordinate system
        /// @param cache Cache for the magnetic field provider
        /// @return magnetic field vector
        Acts::Result<Acts::Vector3> getField(
                const Acts::Vector3& position,
                MagneticFieldProvider::Cache& cache) const override;
    
        /// @brief Get the magnetic field gradient at a given position
        ///
        /// @param position Vector3 position in global coordinate system
        /// @param derivative ActsMatrix<3, 3> to store the gradient
        /// @param cache Cache for the magnetic field provider
        /// @return magnetic field gradient vector
        Acts::Result<Acts::Vector3> getFieldGradient(
                const Acts::Vector3& position, Acts::ActsMatrix<3, 3>& derivative,
                MagneticFieldProvider::Cache& cache) const override;

private:
    /// dipole field strength
    Acts::ActsScalar m_fieldS;
};
