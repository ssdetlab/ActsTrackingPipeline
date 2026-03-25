#pragma once

#include <cstddef>

class IAnnealingScheduler {
 public:
  virtual double getAnnealingFactor(std::size_t it) const = 0;
};
