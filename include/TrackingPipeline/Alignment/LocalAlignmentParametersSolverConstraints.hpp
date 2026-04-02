#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/AlignmentMask.hpp"

#include <cstddef>

class LocalAlignmentParametersSolverConstraints {
 public:
  struct Config {
    /// Alignment mask
    ActsAlignment::AlignmentMask alignmentMask;
  };

  LocalAlignmentParametersSolverConstraints(const Config& cfg,
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

  const std::size_t m_nConstraints = 3;

  std::vector<std::size_t> m_activeIdxs;

  std::unique_ptr<const Acts::Logger> m_logger;
};
