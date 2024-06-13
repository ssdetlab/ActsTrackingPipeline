#pragma once

#include "ActsLUXEPipeline/ILookupProvider.hpp"

#include <fstream>

class CsvLookupTableProvider : public ILookupProvider {
    public:
        /// @brief Configuration struct
        struct Config {
            /// The input file path
            std::string filePath = "";
        };

        /// @brief Constructor
        CsvLookupTableProvider(const Config& config)
            : m_cfg(config) {
                if (m_cfg.filePath.empty()) {
                    throw std::invalid_argument("No file path provided");
                }
            };

        /// @brief Virtual destructor
        ~CsvLookupTableProvider() = default;

        /// @brief Interface function to get the lookup table
        /// for par_1 and par_2 values
        std::unordered_map<
            Acts::ActsScalar,Acts::ActsScalar> getLookup() const override {
                std::unordered_map<Acts::ActsScalar, Acts::ActsScalar> lookupTable;
                std::ifstream file;
                file.open(m_cfg.filePath);
                if (!file.is_open()) {
                    throw std::runtime_error("Could not open file: " + m_cfg.filePath);
                }
                std::string line;
                while (std::getline(file, line)) {
                    std::stringstream ss(line);
                    std::string par1;
                    std::string par2;
                    std::getline(ss, par1, ',');
                    std::getline(ss, par2, ',');
                    lookupTable[std::stod(par1)] = std::stod(par2);
                }
                file.close();
                return lookupTable;
        };

    private:
        /// The configuration
        Config m_cfg;
};
