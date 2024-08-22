#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class QuadrupoleMagField : public Acts::MagneticFieldProvider {
    public:
        /// @brief Cache for the magnetic field provider
        ///
        /// No specific cache is needed for the constant 
        /// magnetic field
        struct Cache {
            /// @brief constructor with context
            Cache(const Acts::MagneticFieldContext& /*mcfg*/) {}
        };

        /// @brief Constructor with magnetic field gradient
        ///
        /// @param gradient magnetic field gradient
        QuadrupoleMagField(Acts::ActsScalar gradient);

        ~QuadrupoleMagField() override;

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
            const Acts::Vector3& /*position*/, 
            Acts::ActsMatrix<3, 3>& /*derivative*/,
            MagneticFieldProvider::Cache& /*cache*/) const override {
                Acts::Vector3 loaclGrad(m_gradient, m_gradient, 0);
                return Acts::Result<Acts::Vector3>::success(loaclGrad);
        }

        /// @brief Get the magnetic field cache
        ///
        /// @param mctx Magnetic field context
        /// @return magnetic field cache
        Acts::MagneticFieldProvider::Cache makeCache(
            const Acts::MagneticFieldContext& mctx) const override;

        private:
            Acts::ActsScalar m_gradient = 0.0;
};
