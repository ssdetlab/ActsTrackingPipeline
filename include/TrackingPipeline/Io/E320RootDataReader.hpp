#pragma once

#include "TrackingPipeline/Io/RootDataReader.hpp"
#include "TrackingPipeline/Geometry/E320GeometryConstraints.hpp"
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/EventData/DataContainers.hpp"

#include "Acts/Definitions/Units.hpp"

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

namespace E320Io {

/// @brief The ROOT file reader for the LUXE simulation
/// that knows about the true trueHits and the true momenta
///
/// @note Covariance is implemented as a diagonal matrix
/// of ALPIDE intrinsic resolutions
class E320RootSimDataReader : public RootSimDataReader {
    public:
        /// @brief The configuration struct
        struct Config 
            : public RootSimDataReader::Config{
                /// The geometry options
                E320Geometry::GeometryOptions gOpt;
        };

        E320RootSimDataReader(const Config &config, Acts::Logging::Level level) 
            : RootSimDataReader(config, level), m_cfg(config) {
                m_actsToWorld = 
                    m_cfg.gOpt.actsToWorld.rotation().inverse();
            }

        std::string name() const override { return "E320RootSimDataReader"; }

    private:
        Config m_cfg;

        Acts::RotationMatrix3 m_actsToWorld;

        // Prepare the measurements
        // for the Sequencer pipeline
        inline void prepareMeasurements(
            const AlgorithmContext &context, 
            std::vector<Acts::SourceLink>* sourceLinks,
            SimClusters* clusters) const override {
                // Check if the event number is correct
                auto eventId = m_intColumns.at("eventId");
                if (eventId != context.eventNumber) {
                    return;
                }

                // Columns with the measurable quantities
                std::int32_t geoIdval; 
                std::int32_t sizeX;
                std::int32_t sizeY;
                TVector3* geoCenter;

                // Columns with the truth quantities
                std::vector<std::int32_t>* trackId;
                std::vector<std::int32_t>* parentTrackId;
                std::vector<std::int32_t>* runId;
                std::vector<TVector3>* trueHits;
                std::vector<TVector3>* vertices;
                std::vector<TLorentzVector>* mom;
                std::vector<TLorentzVector>* momIP;
                try {
                    //--------------------------------
                    // Measurable quantities

                    // Geometry ID of the chip where
                    // the cluster occured
                    geoIdval = m_intColumns.at("geoId");

                    // The cluster geometrical center
                    // in the global coordinates
                    geoCenter = m_vector3Columns.at("rglobal_geo");

                    // The cluster size in X
                    sizeX = m_intColumns.at("xsize");

                    // The cluster size in Y
                    sizeY = m_intColumns.at("ysize");

                    //--------------------------------
                    // Truth quantities

                    // Track IDs of the particles 
                    // that created the cluster
                    trackId = m_vIntColumns.at("tru_trackId");

                    // Parent track IDs of the particles
                    // that created the cluster
                    parentTrackId = m_vIntColumns.at("tru_parenttrackId");

                    // Ptarmigan run ID of the particles
                    // that created the cluster
                    runId = m_vIntColumns.at("tru_runId");

                    // The true hit positions of the particles
                    // that created the cluster
                    trueHits = m_vVector3Columns.at("tru_hit");

                    // The true vertex positions of the particles
                    // that created the cluster
                    vertices = m_vVector3Columns.at("tru_vertex");

                    // The true momenta of the particles
                    // that created the cluster
                    mom = m_vLorentzColumns.at("tru_p");

                    // The true momenta of the particles
                    // that created the cluster at the IP
                    momIP = m_vLorentzColumns.at("tru_p_ip");
                } 
                catch (const std::out_of_range& e) {
                    throw std::runtime_error("Missing columns in the ROOT file");
                }

                //-------------------------------
                // Measurable quantities

                // Apply the Geometry ID convention
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

                // Convert the cluster geometrical center
                // to the local coordinates
                Acts::Vector3 hitGlob = 
                    {geoCenter->X() * Acts::UnitConstants::mm, 
                    geoCenter->Y() * Acts::UnitConstants::mm, 
                    geoCenter->Z() * Acts::UnitConstants::mm};
                
                const Acts::Vector2 hitLoc =
                    convertToLoc(hitGlob, geoId, m_cfg.gOpt);

                // Estimate error from the cluster size
                double pixSizeX = 27_um;
                double pixSizeY = 29_um;

                Acts::Vector2 stddev(sizeX * pixSizeX, sizeY * pixSizeY);
                Acts::SquareMatrix2 cov = stddev.cwiseProduct(stddev).asDiagonal();

                // Fill the measurement
                SimpleSourceLink ssl(hitLoc, cov, geoId, eventId, sourceLinks->size());

                //-------------------------------
                // Truth quantities

                // Fill the cluster data
                SimCluster cluster{
                    .sourceLink = ssl
                };

                if(!m_cfg.clusterFilter->operator()(
                    context.geoContext, cluster)) {
                        return;
                }

                sourceLinks->push_back(Acts::SourceLink(ssl));
                Acts::ActsScalar me = 0.511 * Acts::UnitConstants::MeV;
                for (int idx = 0; idx < trueHits->size(); idx++) {
                    auto hitMom = mom->at(idx); 

                    // Convert the vertex
                    Acts::Vector3 trueVertex3 = 
                        {vertices->at(idx).X() * Acts::UnitConstants::mm, 
                        vertices->at(idx).Y() * Acts::UnitConstants::mm, 
                        vertices->at(idx).Z() * Acts::UnitConstants::mm};
                    trueVertex3 = 
                        m_actsToWorld * trueVertex3;

                    // KF accepts 4D vectors
                    Acts::Vector4 trueVertex = 
                        {trueVertex3.x(), 
                        trueVertex3.y(), 
                        trueVertex3.z(), 0};

                    // Convert the true hit to the local coordinates
                    Acts::Vector3 trueHitGlob = 
                        {trueHits->at(idx).X() * Acts::UnitConstants::mm, 
                        trueHits->at(idx).Y() * Acts::UnitConstants::mm, 
                        trueHits->at(idx).Z() * Acts::UnitConstants::mm};

                    const Acts::Vector2 trueHitLoc = 
                        convertToLoc(trueHitGlob, geoId, m_cfg.gOpt);

                    // Set up the truth parameters
                    Acts::BoundVector truthPars = Acts::BoundVector::Zero();
                    truthPars[Acts::eBoundLoc0] = trueHitLoc[Acts::eBoundLoc0];
                    truthPars[Acts::eBoundLoc1] = trueHitLoc[Acts::eBoundLoc1];

                    // Get the momentum at the first hit
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

                    SimpleSourceLink truthSsl(
                        trueHitLoc, cov, geoId, eventId, clusters->size());

                    // Fill the cluster data
                    cluster.truthHits.push_back(
                        SimHit(
                            Acts::SourceLink(truthSsl),
                            truthPars,
                            ipParameters,
                            trackId->at(idx),
                            parentTrackId->at(idx),
                            runId->at(idx)
                        ));
                    cluster.isSignal = m_intColumns.at("isSignal");
                    cluster.index = ssl.index();
                }
                
                clusters->push_back(cluster);
        };
};

auto defaultSimConfig() {
    E320RootSimDataReader::Config config;
    config.treeName = "clusters";
    config.vVector3Keys = {"tru_hit", "tru_vertex"};
    config.vector3Keys = {"rglobal_geo"};
    config.vLorentzKeys = {"tru_p", "tru_p_ip"};
    config.vIntKeys = {"tru_trackId", "tru_parenttrackId", "tru_runId"};
    config.intKeys = {"eventId", "geoId", "xsize", "ysize", "isSignal"};
    return config;
}

} // namespace E320Io
