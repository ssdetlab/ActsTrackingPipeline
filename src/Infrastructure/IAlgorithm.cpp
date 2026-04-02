#include "TrackingPipeline/Infrastructure/IAlgorithm.hpp"

#include "Acts/Utilities/Logger.hpp"

IAlgorithm::IAlgorithm(const std::string& name, Acts::Logging::Level level)
    : m_name(name), m_logger(Acts::getDefaultLogger(m_name, level)) {}

std::string IAlgorithm::name() const {
  return m_name;
}
