#pragma once

#include <vector>

#include "TLorentzVector.h"
#include "TMatrixDfwd.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"

/// @brief Dummy namespace for ROOT dictionary generation
namespace TrackingPipelineDummies {
using vectorVector2 = std::vector<TVector2>;
using vectorVector3 = std::vector<TVector3>;
using vectorLorentzVector = std::vector<TLorentzVector>;
using vectorMatrixD = std::vector<TMatrixD>;
using vectorVectorD = std::vector<TVectorD>;
}  // namespace TrackingPipelineDummies
