#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsLUXEPipeline/IMagneticField.hpp"

namespace LUXEMagneticField {
using namespace Acts::UnitLiterals;

namespace MagneticFields {

    struct ExampleDipole : public IMagneticField {
        Acts::ActsScalar y_Min;
        Acts::ActsScalar y_Max;
        Acts::ActsScalar B_z;

        ExampleDipole(std::pair<Acts::ActsScalar,Acts::ActsScalar> range,
                      Acts::ActsScalar bz) :
            y_Min(range.first),
            y_Max(range.second),
            B_z(bz) {}

        Acts::Vector3 getValue(const std::array<Acts::ActsScalar, 3> &pos) const final {
            Acts::ActsScalar y = pos.at(1);

            if (y < y_Min || y > y_Max) {
                return Acts::Vector3(0., 0., 0.);
            }
            return Acts::Vector3(0., 0., B_z*1_T);
        }
    };
}
}