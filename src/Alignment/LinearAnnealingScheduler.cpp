#include "TrackingPipeline/Alignment/LinearAnnealingScheduler.hpp"

LinearAnnealingScheduler::LinearAnnealingScheduler(const Config& cfg)
    : m_cfg(cfg) {
  m_b = m_cfg.alphaStart;
  m_a = (m_cfg.nIt != 1) ? (m_cfg.alphaEnd - m_b) / (m_cfg.nIt - 1) : 0;
}

double LinearAnnealingScheduler::getAnnealingFactor(std::size_t it) const {
  return m_a * it + m_b;
}
