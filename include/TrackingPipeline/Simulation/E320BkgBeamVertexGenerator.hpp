#pragma once

#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/Simulation/IVertexGenerator.hpp"

#include <Acts/Definitions/Algebra.hpp>
#include <random>

namespace E320Sim {

using namespace Acts::UnitLiterals;

/// @brief Class that samples vertex from a ROOT histogram
struct E320BkgBeamVertexGenerator : public IVertexGenerator {
    Acts::ActsScalar y0 = 90.3 - 29.94176/2;
    Acts::ActsScalar y1 = 331.1 + 29.94176/2;
    
    Acts::ActsScalar x0 = 90.3 - 29.94176/2;
    Acts::ActsScalar x1 = 331.1 + 29.94176/2;

    Acts::ActsScalar yPower = 69.8048;
    Acts::ActsScalar yShift = -22165;

    Acts::Vector3 gen(RandomEngine& rng) const override {
        Acts::ActsScalar y0 = m_gOpt.chipY.at(0) - m_gOpt.chipSizeY/2;
        Acts::ActsScalar y1 = m_gOpt.chipY.at(8) + m_gOpt.chipSizeY/2;

        Acts::ActsScalar x0 = m_gOpt.chipX - m_gOpt.chipSizeX/2;
        Acts::ActsScalar x1 = m_gOpt.chipX + m_gOpt.chipSizeX/2;

        std::discrete_distribution<> discrete({0.25, 0.23, 0.24, 0.28});
        std::uniform_real_distribution<> uniform(0, 1);

        // Generate unfiorm x
        Acts::ActsScalar x = x0 + (x1 - x0) * uniform(rng);

        // Generate power-law y
        Acts::ActsScalar gamma = uniform(rng);

        Acts::ActsScalar yTerm1 = gamma * std::pow(y1 - yShift, yPower + 1);
        Acts::ActsScalar yTerm2 = (1 - gamma) * std::pow(y0 - yShift, yPower + 1);

        Acts::ActsScalar y = yShift + std::pow(yTerm1 + yTerm2, 1./(yPower + 1));

        // Generate uniform z
        Acts::ActsScalar z = m_gOpt.staveZ.at(discrete(rng));

        Acts::Vector3 glob(x, y, z);
        glob = m_gOpt.actsToWorld.rotation().inverse() * glob;

        return glob;
    }

    E320Geometry::GeometryOptions m_gOpt;
};

} // namespace E320Sim