#include "ActsLUXEPipeline/IAlgorithm.hpp"

#include "Acts/Utilities/Logger.hpp"

#include <utility>


IAlgorithm::IAlgorithm(std::string name, Acts::Logging::Level level)
    : m_name(std::move(name)),
    m_logger(Acts::getDefaultLogger(m_name, level)) {}

std::string IAlgorithm::name() const {
    return m_name;
}
