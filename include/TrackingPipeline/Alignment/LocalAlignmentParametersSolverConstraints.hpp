#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <cstddef>

/// @brief Local alignment parameters solver constraining alignment
/// delta to have zero average transverse shift and zero average rotation
/// about the surface's normal
///
/// The solver takes the Acts Alignment Hessian and chi2 derivative and
/// constraints them such that the resulting delta of the alignment parameters
/// for the alignment surfaces have zero mean transverse shift and rotation
/// about normal
class LocalAlignmentParametersSolverConstraints {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Alignment mask
    ActsAlignment::AlignmentMask alignmentMask;
  };

  /// @brief Constructor
  ///
  /// @param cfg configuration struct
  /// @param level logging level
  LocalAlignmentParametersSolverConstraints(const Config& cfg,
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

  /// Number of alignment step constraints
  const std::size_t m_nConstraints = 3;

  /// Vector of active alignment dofs indices
  std::vector<std::size_t> m_activeIdxs;

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};
