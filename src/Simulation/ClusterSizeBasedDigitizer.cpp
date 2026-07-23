#include "TrackingPipeline/Simulation/ClusterSizeBasedDigitizer.hpp"

#include <cstddef>
#include <random>
#include <unordered_set>

ClusterSizeBasedDigitizer::ClusterSizeBasedDigitizer(const Config& config)
    : m_cfg(config) {
  std::unordered_set<std::size_t> clSizes;
  std::size_t maxClSize = 0;
  for (const auto& [clSize, probStdDev] : m_cfg.clSizeProbsStdDevs) {
    clSizes.insert(clSize);
    if (clSize > maxClSize) {
      maxClSize = clSize;
    }
  }
  for (std::size_t i = 0; i <= maxClSize; i++) {
    if (m_cfg.clSizeProbsStdDevs.contains(i)) {
      m_clSizeProbs.push_back(std::get<0>(m_cfg.clSizeProbsStdDevs.at(i)));
    } else {
      m_clSizeProbs.push_back(0);
    }
  }

  m_discrete = std::discrete_distribution<std::size_t>(m_clSizeProbs.begin(),
                                                       m_clSizeProbs.end());
}

std::pair<Acts::Vector2, Acts::SquareMatrix2>
ClusterSizeBasedDigitizer::digitize(RandomEngine& rng,
                                    const Acts::GeometryIdentifier& /*geoId*/,
                                    const Acts::Vector2& pos) const {
  std::size_t clSize = m_discrete(rng);

  double resolutionX = std::get<1>(m_cfg.clSizeProbsStdDevs.at(clSize));
  double resolutionY = std::get<2>(m_cfg.clSizeProbsStdDevs.at(clSize));
  Acts::Vector2 stdDev = {resolutionX, resolutionY};
  Acts::Vector2 digLocal =
      pos + stdDev.cwiseProduct(Acts::Vector2(m_normal(rng), m_normal(rng)));
  return {digLocal, stdDev.cwiseProduct(stdDev).asDiagonal()};
}
