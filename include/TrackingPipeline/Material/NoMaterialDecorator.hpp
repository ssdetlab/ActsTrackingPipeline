#pragma once

#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"

class NoMaterialDecorator : public Acts::IMaterialDecorator {
    public:
        /// Configuration struct for the NoMaterialDecorator
        struct Config {
            /// The surface binnings
            Acts::BinUtility surfaceBinning;
            /// Vetos for the material provider
            std::vector<Acts::GeometryIdentifier> vetos;
        };

        /// Default destructor
        ~NoMaterialDecorator() override = default;

        /// Constructor
        NoMaterialDecorator(const Config& cfg) : m_cfg(cfg) {}

        /// Decorate a surface
        ///
        /// @param surface the non-cost surface that is decorated
        void decorate(Acts::Surface& surface) const override {
            if (std::find(m_cfg.vetos.begin(), 
                    m_cfg.vetos.end(), 
                    surface.geometryId()) != m_cfg.vetos.end()) {
                        return;
            }
            auto surfaceMaterial = 
                std::make_shared<Acts::ProtoSurfaceMaterial>(m_cfg.surfaceBinning);
            surface.assignSurfaceMaterial(surfaceMaterial);
        };

        /// Decorate a TrackingVolume
        ///
        /// @param volume the non-cost volume that is decorated
        void decorate(Acts::TrackingVolume& volume) const override {};

    private:
        /// Private access to the configuration
        Config m_cfg;
};