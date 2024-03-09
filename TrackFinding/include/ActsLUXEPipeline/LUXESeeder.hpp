#pragma once

#include "Acts/EventData/SourceLink.hpp"

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
#include "ActsLUXEPipeline/LUXEMeasurement.hpp"

struct Seed {
    LUXEMeasurement::SimMeasurements measurements;
};

using Seeds = std::vector<Seed>;

class IdealSeeder : public IAlgorithm {
    public:
        struct Config {
            std::string inputSourceLinks = "SourceLink";
            double roadWidth = 20;
        };

        IdealSeeder(Config config, Acts::Logging::Level level)
            : IAlgorithm("IdealSeeder", level),
            m_cfg(std::move(config)) {
                m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
        }
        ~IdealSeeder() = default;

        ProcessCode execute(const AlgorithmContext& ctx) const override {
            auto input = m_inputSourceLinks(ctx);

            std::sort(input.begin(), input.end(),
                [](const auto& a, const auto& b) {
                    return a.trackId < b.trackId;
                }
            );

            Seeds seeds;
            LUXEMeasurement::SimMeasurements measurements{input.at(0)};
            for (auto it = input.begin() + 1; it != input.end(); ++it) {
                SimpleSourceLink dummy = it->sourceLink.get<SimpleSourceLink>();
                // dummy.covariance
                if (it->trackId == (it - 1)->trackId) {
                    measurements.push_back(*it);
                }
                else {
                    seeds.push_back(Seed{measurements});
                    measurements.clear();
                    measurements.push_back(*it);
                }
            }
            seeds.push_back(Seed{measurements});

            double avgSeedSize = 0;
            for (auto& seed : seeds) {
                avgSeedSize += seed.measurements.size();
            }
            avgSeedSize /= seeds.size();

            return ProcessCode::SUCCESS;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }
    private:
        Config m_cfg;

        ReadDataHandle<LUXEMeasurement::SimMeasurements> m_inputSourceLinks{
            this, "InputSourceLinks"};

};
