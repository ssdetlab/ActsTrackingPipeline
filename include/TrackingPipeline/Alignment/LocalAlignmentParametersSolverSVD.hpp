#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"

#include <cstddef>

/// @brief Local alignment parameters solver constraining alignment
/// delta through Hessian mode suppression
///
/// The solver takes the Acts Alignment Hessian and performs its SVD
/// decomposition. The SVD modes are zeroed-out based on their singular values
/// magnitude relative to the greates singular value in the decomposition. The
/// Hessian with low singular value modes removed is inverted and the alignment
/// step is solved for.
class LocalAlignmentParametersSolverSVD {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Tolerance realtive to the max singular value
    double maxSingularValueTol;
    /// Tolerance realtive to the singular value gap
    double singularValueGapTol;
    /// Alignment mask
    ActsAlignment::AlignmentMask alignmentMask;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  /// @param level logging level
  LocalAlignmentParametersSolverSVD(const Config& cfg,
                                    Acts::Logging::Level level);

  /// @brief calculate alignment parameters
  ///
  /// @param gctx current geometry context
  /// @param alignRes alignment result to fill
  /// @param sumChi2Derivative first alignment chi2 derivative
  /// @param sumChi2SecondDerivative alignment chi2 Hessian
  void calculateAlignmentParameters(
      const Acts::GeometryContext& gctx,
      ActsAlignment::AlignmentResult& alignRes,
      const Acts::ActsDynamicVector& sumChi2Derivative,
      const Acts::ActsDynamicMatrix& sumChi2SecondDerivative) const;

 protected:
  /// @brief access to logging instance
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  /// Configuration
  Config m_cfg;

  /// Vector of active alignment dofs indices
  std::vector<std::size_t> m_activeIdxs;

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};
