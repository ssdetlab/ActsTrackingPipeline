#pragma once

#include "TVector3.h"
#include "TLorentzVector.h"

#include <vector>

/// @brief Dummy namespace for ROOT dictionary generation
namespace ActsLUXEPipeline {
    using vectorVector3 = std::vector<TVector3>;
    using vectorLorentzVector = std::vector<TLorentzVector>;
} // namespace ActsLUXEPipeline
