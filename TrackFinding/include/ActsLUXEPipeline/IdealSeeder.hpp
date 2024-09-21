#pragma once

#include "ActsLUXEPipeline/IAlgorithm.hpp"
#include "ActsLUXEPipeline/DataHandle.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"

#include "Acts/EventData/SourceLink.hpp"

/// @brief The ideal seeder for the LUXE simulation
/// takes the the SimMeasurements and converts them 
/// into seeds
class IdealSeeder {
    public:
        /// @brief Delegate to estimate the IP parameters
        /// and the momentum direction at the first tracking layer
        ///
        /// @arg The geometry context to use
        /// @arg The global position of the pivot source link
        ///
        /// @return Particle charge, the IP momentum magnitude, the IP vertex position,
        /// the IP momentum direction, the momentum direction at the
        /// first tracking layer
        using TrackEstimator =
            Acts::Delegate<std::tuple<
                Acts::ActsScalar, Acts::ActsScalar, Acts::Vector3, Acts::Vector3, Acts::Vector3>(
                    const Acts::GeometryContext&, const Acts::Vector3&)>;

        /// @brief Delegate to transform the source link to the
        /// appropriate global frame.
        ///
        /// @arg The geometry context to use
        /// @arg The source link to calibrate
        ///
        /// @return The global position of the source link measurement
        using SourceLinkCalibrator =
            Acts::Delegate<Acts::Vector3(const Acts::GeometryContext&, const Acts::SourceLink&)>;

        /// @brief The nested configuration struct
        struct Config {
            /// Parameters estimator
            TrackEstimator trackEstimator;
            /// SourceLink calibrator
            SourceLinkCalibrator sourceLinkCalibrator;
            /// First layer IDs
            std::vector<Acts::GeometryIdentifier> firstLayerIds;
            /// Minimum number of source links
            int minSourceLinks = 3;
            /// Maximum number of source links
            int maxSourceLinks = 10;
        };

        /// @brief Constructor
        IdealSeeder(Config config) : m_cfg(std::move(config)) {}
        ~IdealSeeder() = default;

        /// @brief The execute method        
        Seeds getSeeds(
            const Acts::GeometryContext& gctx, 
            const SimMeasurements& measurements) const {
                using namespace Acts::UnitLiterals;
    
                // Create the seeds
                Seeds seeds;
                std::vector<Acts::SourceLink> sourceLinks;

                // Create IP covariance matrix from
                // reasonable standard deviations
                Acts::BoundVector ipStdDev;
                ipStdDev[Acts::eBoundLoc0] = 100_um;
                ipStdDev[Acts::eBoundLoc1] = 100_um;
                ipStdDev[Acts::eBoundTime] = 25_ns;
                ipStdDev[Acts::eBoundPhi] = 2_degree;
                ipStdDev[Acts::eBoundTheta] = 2_degree;
                ipStdDev[Acts::eBoundQOverP] = 1 / 100_GeV;
                Acts::BoundSquareMatrix ipCov =
                    ipStdDev.cwiseProduct(ipStdDev).asDiagonal();

                // Insert the first source link
                sourceLinks.push_back(measurements.front().sourceLink);
                for (auto it = measurements.begin() + 1; it != measurements.end(); ++it) {
                    if (it->trackId == (it - 1)->trackId && (it + 1) != measurements.end()) {
                        // Add the source link to the list
                        // if the hit is from the same track
                        sourceLinks.push_back(it->sourceLink);
                    }
                    else {
                        if ((it + 1) == measurements.end()) {
                            // Add the last source link
                            sourceLinks.push_back(it->sourceLink);
                        }

                        if (sourceLinks.size() < m_cfg.minSourceLinks ||
                            sourceLinks.size() > m_cfg.maxSourceLinks) {
                                // Reset the source links
                                sourceLinks.clear();
                                sourceLinks.push_back(it->sourceLink);

                                continue;
                        }

                        bool pivotFound = false;
                        Acts::SourceLink pivot = sourceLinks.front();
                        for (auto& sl : sourceLinks) {
                            auto ssl = sl.get<SimpleSourceLink>();
                            auto sourceLinkId = ssl.geometryId();
                            
                            if (std::find(
                                    m_cfg.firstLayerIds.begin(), 
                                    m_cfg.firstLayerIds.end(), 
                                    sourceLinkId) != m_cfg.firstLayerIds.end()) {
                                        pivot = sl;
                                        pivotFound = true;
                                        break;
                            }
                        }

                        if (pivotFound) {
                            if (m_cfg.trackEstimator.connected()) {
                                Acts::Vector3 globalPos = m_cfg.sourceLinkCalibrator(gctx, pivot);
            
                                // Get the IP parameters
                                auto [q, ipP, ipVertex, ipDir, flDir] =
                                    m_cfg.trackEstimator(gctx, globalPos);
            
                                Acts::Vector4 mPos4 = {ipVertex.x(), ipVertex.y(), ipVertex.z(), 0};
            
                                Acts::ActsScalar p = ipP; 
                                Acts::ActsScalar theta = std::acos(ipDir.z()/p);
                                Acts::ActsScalar phi = std::atan2(ipDir.y(), ipDir.x());
                
                                Acts::CurvilinearTrackParameters ipParameters(
                                    mPos4, phi, theta,
                                    1_e / p, ipCov, 
                                    Acts::ParticleHypothesis::electron());
    
                                // Add the seed to the list
                                seeds.push_back(Seed{
                                    sourceLinks, ipParameters, (it - 1)->trackId});
                            }
                            else {
                                // Ip parameter is the same for all hits
                                // with the same track id
                                Acts::CurvilinearTrackParameters ipParameters =
                                    (it - 1)->ipParameters;
            
                                // Add the seed to the list
                                seeds.push_back(Seed{
                                    sourceLinks, ipParameters, (it - 1)->trackId});
                            }
                        }
                        // Reset the source links
                        sourceLinks.clear();
                        sourceLinks.push_back(it->sourceLink);
                    }
                }
    
                return seeds;
        }

        /// Get readonly access to the config parameters
        const Config& config() const { return m_cfg; }

    private:
    
        Config m_cfg;
};
