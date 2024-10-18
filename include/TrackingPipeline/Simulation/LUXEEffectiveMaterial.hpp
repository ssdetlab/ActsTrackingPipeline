#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/MaterialSlab.hpp"

#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"

/// Utility function to build 
/// a particle from the dataset parameters
ActsFatras::Particle makeParticle(
    Acts::PdgParticle pdg, 
    Acts::Vector3 dir, 
    double p, 
    Acts::Vector4 pos4) {
        const auto id = ActsFatras::Barcode().setVertexPrimary(1).setParticle(1);
        return ActsFatras::Particle(id, pdg)
            .setPosition4(pos4)
            .setDirection(dir)
            .setAbsoluteMomentum(p);
}

/// Utility function to build 
/// a mockup material slab used to 
/// approximate scattering and energy loss
Acts::MaterialSlab makeSiliconSlab() {
    using namespace Acts::UnitLiterals;
    return {Acts::Material::fromMolarDensity(
            9.370_cm, 46.52_cm, 28.0855, 14,
            (2.329 / 28.0855) * 1_mol / 1_cm3), 
            100_um};
}