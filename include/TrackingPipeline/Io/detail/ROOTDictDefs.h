#pragma once

#include <vector>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector3.h"

/// @brief Dummy namespace for ROOT dictionary generation
namespace TrackingPipelineDummies {
using vectorVector3 = std::vector<TVector3>;
using vectorLorentzVector = std::vector<TLorentzVector>;
}  // namespace TrackingPipelineDummies
