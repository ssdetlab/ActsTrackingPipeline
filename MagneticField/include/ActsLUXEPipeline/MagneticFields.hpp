#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsLUXEPipeline/IMagneticField.hpp"

namespace LUXEMagneticField {
using namespace Acts::UnitLiterals;

namespace MagneticFields {

    struct ExampleDipole : public IMagneticField {
        Acts::Vector3 getValue(const std::array<Acts::ActsScalar, 3> &pos) const final {
            Acts::ActsScalar y = pos.at(1);
            std::cout<<pos.at(0)<<" "<<pos.at(1)<<" "<<pos.at(2)<<std::endl;
            if (y < 1450_mm || y > 2650_mm) {
                return Acts::Vector3(0., 0., 0.);
            }
            return Acts::Vector3(0., 0., .95_T);
        }
    };
}
}