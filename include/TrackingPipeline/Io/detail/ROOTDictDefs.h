#pragma once

#include <vector>

#include "EpicsFrame.h"
#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TVector3.h"

/// @brief Dummy namespace for ROOT dictionary generation
namespace TrackingPipelineDummies {
using vectorVector3 = std::vector<TVector3>;
using vectorVector2 = std::vector<TVector2>;
using vectorLorentzVector = std::vector<TLorentzVector>;
using vectorMatrixD = std::vector<TMatrixD>;
using EPICSFrame = EpicsFrame;
}  // namespace TrackingPipelineDummies
