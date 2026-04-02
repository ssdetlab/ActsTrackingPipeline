#pragma once

#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"

#include <cstddef>

class LocalAlignmentParametersSolverSVD {
 public:
  struct Config {
    /// Tolerance realtive to the max singular value
    double maxSingularValueTol;
    /// Tolerance realtive to the singular value gap
    double singularValueGapTol;
    /// Alignment mask
    ActsAlignment::AlignmentMask alignmentMask;
  };

  LocalAlignmentParametersSolverSVD(const Config& cfg,
                                    Acts::Logging::Level level);

  void calculateAlignmentParameters(
      const Acts::GeometryContext& /*gctx*/,
      ActsAlignment::AlignmentResult& alignRes,
      const Acts::ActsDynamicVector& sumChi2Derivative,
      const Acts::ActsDynamicMatrix& sumChi2SecondDerivative) const;

 protected:
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  Config m_cfg;

  std::vector<std::size_t> m_activeIdxs;

  std::unique_ptr<const Acts::Logger> m_logger;
};
