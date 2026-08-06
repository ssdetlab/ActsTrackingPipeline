#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"

#include <cstddef>

/// @brief Global alignment parameters solver performing alignment delta
/// projection onto the rigid body space
///
/// The solver takes the Acts Alignment Hessian and chi2 derivative and projects
/// them onto the rigid space for solving the minimization step
class GlobalAlignmentParametersSolverProjection {
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
  GlobalAlignmentParametersSolverProjection(const Config& cfg,
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
  /// Number of rigid degrees of freedom
  const std::size_t m_nRigidDof = 6;

  /// Configuration
  Config m_cfg;

  /// Vector of active alignment dofs indices
  std::vector<std::size_t> m_activeIdxs;

  /// Flags indicating active rigid degrees of freedom
  bool m_doCenter0;
  bool m_doCenter1;
  bool m_doCenter2;
  bool m_doAngle0;
  bool m_doAngle1;
  bool m_doAngle2;

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};
