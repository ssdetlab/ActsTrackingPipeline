#pragma once

#include "ActsLUXEPipeline/ROOTDataReader.hpp"
#include "ActsLUXEPipeline/E320GeometryConstraints.hpp"
#include "ActsLUXEPipeline/SimpleSourceLink.hpp"
#include "ActsLUXEPipeline/DataContainers.hpp"

#include "Acts/Definitions/Units.hpp"

namespace E320ROOTReader {

using namespace Acts::UnitLiterals;

/// @brief Global to local conversion
/// for the LUXE geometry
Acts::Vector2 convertToLoc(
    const Acts::Vector3& glob, 
    const Acts::GeometryIdentifier geoId,
    const E320Geometry::GeometryOptions& gOpt) {
        int nStave = geoId.sensitive()/10 - 1;
        int nChip = geoId.sensitive()%10 - 1;
        Acts::Vector2 loc = Acts::Vector2(
            (glob.y() - gOpt.chipY.at(nChip)),
            (glob.x() - gOpt.chipX));
        return loc;
}

/// @brief The ROOT file reader for the LUXE simulation
/// that knows about the true hits and the true momenta
///
/// @note Covariance is implemented as a diagonal matrix
/// of ALPIDE intrinsic resolutions
class E320ROOTSimDataReader : 
    public ROOTDataReader<SimMeasurements> {
        public:
            /// @brief The configuration struct
            struct Config 
                : public ROOTDataReader<SimMeasurements>::Config{
                    /// The geometry options
                    E320Geometry::GeometryOptions gOpt;
                    /// Vertex cuts
                    Acts::Extent vertexPosExtent;
                    /// Energy cuts
                    std::pair<Acts::ActsScalar, Acts::ActsScalar> energyCuts = 
                        {8 * Acts::UnitConstants::GeV, 100 * Acts::UnitConstants::GeV};
            };

            E320ROOTSimDataReader(const Config &config, Acts::Logging::Level level) 
                : ROOTDataReader(config, level), m_cfg(config) {
                    m_actsToWorld = 
                        m_cfg.gOpt.actsToWorld.rotation().inverse();
                    
                    m_energyCuts = m_cfg.energyCuts;
                }
    
            std::string name() const override { return "E320ROOTSimDataReader"; }

        private:
            Config m_cfg;

            Acts::RotationMatrix3 m_actsToWorld;

            std::pair<Acts::ActsScalar, Acts::ActsScalar> m_energyCuts;

            inline void prepareMeasurements(
                const AlgorithmContext &context, 
                SimMeasurements* measurements) const override {
                    // Check if the event number is correct
                    auto eventId = m_intColumns.at("eventId");
                    if (eventId != context.eventNumber) {
                        return;
                    }
    
                    // Get the columns with the event data
                    std::int32_t geoIdval; 
                    std::vector<std::int32_t>* trackId;
                    std::vector<std::int32_t>* parentTrackId;
                    std::vector<TVector3>* hits;
                    std::vector<TVector3>* vertices;
                    std::vector<TLorentzVector>* mom;
                    std::vector<TLorentzVector>* momIP;
                    try {
                        geoIdval = m_intColumns.at("geoId");
                        trackId = m_vectorIntColumns.at("tru_trackId");
                        parentTrackId = m_vectorIntColumns.at("tru_parenttrackId");
                        hits = m_vector3Columns.at("tru_hit");
                        vertices = m_vector3Columns.at("tru_vertex");
                        mom = m_lorentzColumns.at("tru_p");
                        momIP = m_lorentzColumns.at("tru_p_ip");
                    } 
                    catch (const std::out_of_range& e) {
                        throw std::runtime_error("Missing columns in the ROOT file");
                    }

                    // This is per the LUXE convention
                    Acts::GeometryIdentifier geoId;
                    geoId.setSensitive(geoIdval + 11);
                    
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

                    // Create the measurements
                    Acts::ActsScalar me = 0.511 * Acts::UnitConstants::MeV;
                    for (int idx = 0; idx < hits->size(); idx++) {
                        auto hitMom = mom->at(idx); 

                        // Apply the cuts
                        if (hitMom.E() < m_energyCuts.first || 
                            hitMom.E() > m_energyCuts.second) {
                            continue;
                        }

                        // Convert the vertex
                        Acts::Vector3 trueVertex3 = 
                            {vertices->at(idx).X() * Acts::UnitConstants::mm, 
                            vertices->at(idx).Y() * Acts::UnitConstants::mm, 
                            vertices->at(idx).Z() * Acts::UnitConstants::mm};
                        trueVertex3 = 
                            m_actsToWorld * trueVertex3;

                        if (!m_cfg.vertexPosExtent.contains(trueVertex3)) {
                            continue;
                        }

                        // Convert the true hit to the local coordinates
                        Acts::Vector3 trueHitGlob = 
                            {hits->at(idx).X() * Acts::UnitConstants::mm, 
                            hits->at(idx).Y() * Acts::UnitConstants::mm, 
                            hits->at(idx).Z() * Acts::UnitConstants::mm};

                        const Acts::Vector2 trueHitLoc = 
                            convertToLoc(trueHitGlob, geoId, m_cfg.gOpt);

                        // KF accepts 4D vectors
                        Acts::Vector4 trueVertex = 
                            {trueVertex3.x(), 
                            trueVertex3.y(), 
                            trueVertex3.z(), 0};

                        // Set up the truth parameters
                        Acts::BoundVector truthPars = Acts::BoundVector::Zero();
                        truthPars[Acts::eBoundLoc0] = trueHitLoc[Acts::eBoundLoc0];
                        truthPars[Acts::eBoundLoc1] = trueHitLoc[Acts::eBoundLoc1];

                        // Infer as the momentum at the first hit
                        Acts::Vector3 trueP = 
                            {hitMom.Px(), hitMom.Py(), hitMom.Pz()};

                        Acts::Vector3 dir = trueP.normalized();

                        Acts::Vector3 dirRotated = 
                            m_actsToWorld * dir;

                        truthPars[Acts::eBoundPhi] = Acts::VectorHelpers::phi(dirRotated);
                        truthPars[Acts::eBoundTheta] = Acts::VectorHelpers::theta(dirRotated);
                        truthPars[Acts::eBoundQOverP] = 
                            1_e/(hitMom.P() * Acts::UnitConstants::GeV);
                        truthPars[Acts::eBoundTime] = hitMom.T();

                        // Set up IP information
                        TLorentzVector ipMom = momIP->at(idx);

                        Acts::Vector3 ipMom3 = 
                            {ipMom.Px(), ipMom.Py(), ipMom.Pz()};
                        Acts::Vector3 dirIP = ipMom3.normalized();

                        Acts::Vector3 dirIPRotated = 
                            m_actsToWorld * dirIP;

                        Acts::CurvilinearTrackParameters ipParameters(
                            trueVertex, 
                            Acts::VectorHelpers::phi(dirIPRotated),
                            Acts::VectorHelpers::theta(dirIPRotated),
                            1_e/(ipMom.P() * Acts::UnitConstants::GeV),
                            ipCov,
                            Acts::ParticleHypothesis::electron());

                        // Covariance is filled with the intrinsic resolution
                        Acts::Vector2 stddev(5  * Acts::UnitConstants::um,
                            5  * Acts::UnitConstants::um);
                        Acts::SquareMatrix2 cov = stddev.cwiseProduct(stddev).asDiagonal();
        
                        // Create the source link
                        SimpleSourceLink ssl(trueHitLoc, cov, geoId, eventId);
                        Acts::SourceLink sl{ssl};

                        // Fill the measurement
                        SimMeasurement measurement{
                            .sourceLink = sl,
                            .truthParameters = truthPars,
                            .ipParameters = ipParameters,
                            .trackId = parentTrackId->at(idx)
                        };

                        measurements->push_back(measurement);
                    }
        };
};

auto defaultSimConfig() {
    E320ROOTSimDataReader::Config config;
    config.treeName = "clusters";
    config.vector3Keys = {"tru_hit", "tru_vertex"};
    config.lorentzKeys = {"tru_p", "tru_p_ip"};
    config.vectorIntKeys = {"tru_trackId", "tru_parenttrackId"};
    config.intKeys = {"eventId", "geoId"};
    return config;
}

}  // namespace E320ROOTReader
