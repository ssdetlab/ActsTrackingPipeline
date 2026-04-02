#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"

#include <cstddef>

class GlobalAlignmentParametersSolver {
 public:
  struct Config {
    /// Alignment mask
    ActsAlignment::AlignmentMask alignmentMask;
  };

  GlobalAlignmentParametersSolver(const Config& cfg,
                                  Acts::Logging::Level level);

  void calculateAlignmentParameters(
      const Acts::GeometryContext& gctx,
      ActsAlignment::AlignmentResult& alignRes,
      const Acts::ActsDynamicVector& sumChi2Derivative,
      const Acts::ActsDynamicMatrix& sumChi2SecondDerivative) const;

 protected:
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  const std::size_t m_nRigidDof = 6;

  Config m_cfg;

  std::vector<std::size_t> m_activeIdxs;

  bool m_doCenter0;
  bool m_doCenter1;
  bool m_doCenter2;
  bool m_doAngle0;
  bool m_doAngle1;
  bool m_doAngle2;

  std::unique_ptr<const Acts::Logger> m_logger;
};
