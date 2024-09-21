#pragma once

#include <fstream>

class CsvLookupTableProvider {
    public:
        using TrackParameters = std::tuple<
            Acts::ActsScalar, 
            Acts::ActsScalar, 
            Acts::Vector3, 
            Acts::Vector3, 
            Acts::Vector3>;

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

                std::ifstream file;
                file.open(m_cfg.filePath);

                if (!file.is_open()) {
                    throw std::runtime_error("Could not open file: " + m_cfg.filePath);
                }

                std::string line;
                while (std::getline(file, line)) {
                    auto tokens = split(line, ",");

                    m_lookupTable[{std::stod(tokens.at(0)), std::stod(tokens.at(1))}] = {
                        std::stod(tokens.at(2)), 
                        std::stod(tokens.at(3)), 
                        Acts::Vector3(
                            std::stod(tokens.at(4)), 
                            std::stod(tokens.at(5)), 
                            std::stod(tokens.at(6))), 
                        Acts::Vector3(
                            std::stod(tokens.at(7)), 
                            std::stod(tokens.at(8)), 
                            std::stod(tokens.at(9))), 
                        Acts::Vector3(
                            std::stod(tokens.at(10)), 
                            std::stod(tokens.at(11)), 
                            std::stod(tokens.at(12)))};
                }

                file.close();
            };

        /// @brief Virtual destructor
        ~CsvLookupTableProvider() = default;

        /// @brief Find closest value in lookup table
        std::pair<Acts::ActsScalar,Acts::ActsScalar> findClosestValue(
            const Acts::ActsScalar& par_1, 
            const Acts::ActsScalar& par_2) const {
                auto start = std::chrono::system_clock::now();

                std::pair<Acts::ActsScalar,Acts::ActsScalar> closestValue;
                Acts::ActsScalar minDistance = std::numeric_limits<Acts::ActsScalar>::max();
                for (const auto& [key, value] : m_lookupTable) {
                    Acts::ActsScalar distance = std::sqrt(
                        std::pow(par_1 - key.first, 2) + 
                        std::pow(par_2 - key.second, 2));
                    if (distance < minDistance) {
                        minDistance = distance;
                        closestValue = {key.first, key.second};
                    }
                }

                auto end = std::chrono::system_clock::now();
                return closestValue;
            };

        /// @brief Interface function to get the lookup table
        /// for par_1 and par_2 values
        TrackParameters operator()(
            const Acts::GeometryContext& /*gctx*/, 
            const Acts::Vector3& pivot) const {
                auto closestValue = findClosestValue(pivot.x(), pivot.z());
                return m_lookupTable.at(closestValue);
        };

        std::vector<std::string> split(std::string& s, const std::string& delimiter) {
            std::vector<std::string> tokens;
            size_t pos = 0;
            std::string token;
            while ((pos = s.find(delimiter)) != std::string::npos) {
                token = s.substr(0, pos);
                tokens.push_back(token);
                s.erase(0, pos + delimiter.length());
            }
            tokens.push_back(s);
        
            return tokens;
        }

    private:
        /// The configuration
        Config m_cfg;

        /// The lookup table
        std::map<
            std::pair<Acts::ActsScalar, Acts::ActsScalar>, 
            TrackParameters> m_lookupTable;
};
