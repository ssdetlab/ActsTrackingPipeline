#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <string>

namespace LUXEMagneticField {

/// MagneticField implementations only need to implement the `getValue` method.
class IMagneticField {
    public:
        /// This function must be implemented by subclasses.
        virtual Acts::Vector3 getValue(const std::array<Acts::ActsScalar , 3> &pos) const = 0;
};
}