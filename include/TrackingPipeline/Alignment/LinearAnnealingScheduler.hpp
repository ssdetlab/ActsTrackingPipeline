#pragma once

#include <cstddef>

#include "TrackingPipeline/Alignment/IAnnealingScheduler.hpp"

class LinearAnnealingScheduler : public IAnnealingScheduler {
 public:
  struct Config {
    double alphaStart;
    double alphaEnd;
    std::size_t nIt;
  };

  explicit LinearAnnealingScheduler(const Config& cfg);

  double getAnnealingFactor(std::size_t it) const override;

 private:
  Config m_cfg;

  double m_a;
  double m_b;
};
