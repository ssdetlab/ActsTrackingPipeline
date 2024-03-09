#pragma once

#include "Acts/Definitions/Units.hpp"
#include "ActsLUXEPipeline/ROOTDataReader.hpp"
#include "ActsLUXEPipeline/LUXEGeometryConstraints.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
#include "ActsLUXEPipeline/LUXEMeasurement.hpp"

namespace LUXEROOTReader {

/// @brief Global to local conversion
/// for the LUXE geometry
Acts::Vector2 convertToLoc(
    const Acts::Vector3& glob, 
    const Acts::GeometryIdentifier geoId,
    const LUXEGeometry::GeometryOptions& gOpt) {
        // TODO: geoIds for the electron arm
        int nStave = geoId.sensitive()/10;
        int nChip = geoId.sensitive()%10;
        if (nStave % 2 == 0) {
            return Acts::Vector2(
                (glob.x() - gOpt.chipTranslationXEven.at(nChip)),
                (glob.y() - gOpt.chipTranslationY));
        }
        else {
            return Acts::Vector2(
                (glob.x() - gOpt.chipTranslationXOdd.at(nChip)),
                (glob.y() - gOpt.chipTranslationY));
        }
}

/// @brief The ROOT file reader for the LUXE simulation
/// that knows about the true hits and the true momenta
///
/// @note Covariance is implemented as a diagonal matrix
/// of ALPIDE intrinsic resolutions
class LUXEROOTSimDataReader : 
    public ROOTDataReader<LUXEMeasurement::SimMeasurements> {
        public:
            struct Config 
                : public ROOTDataReader<LUXEMeasurement::SimMeasurements>::Config{
                    LUXEGeometry::GeometryOptions gOpt;
            };

            LUXEROOTSimDataReader(const Config &config, Acts::Logging::Level level) 
                : ROOTDataReader(config, level), m_cfg(config) {}
    
            std::string name() const override { return "LUXEROOTSimDataReader"; }

        private:
            Config m_cfg;

            inline void prepareMeasurements(
                const AlgorithmContext &context, 
                LUXEMeasurement::SimMeasurements* measurements) const override {
                    auto eventId = m_intColumns.at("eventId");
                    if (eventId != context.eventNumber) {
                        return;
                    }
    
                    std::int32_t geoIdval; 
                    std::vector<std::int32_t>* trackId;
                    std::vector<TVector3>* hits;
                    std::vector<TLorentzVector>* dirs;
                    try {
                        geoIdval = m_intColumns.at("geoId");
                        trackId = m_vectorIntColumns.at("tru_trackId");
                        hits = m_vector3Columns.at("tru_hit");
                        dirs = m_lorentzColumns.at("tru_p");
                    } catch (const std::out_of_range& e) {
                        throw std::runtime_error("Missing columns in the ROOT file");
                    }

                    Acts::GeometryIdentifier geoId;
                    geoId.setSensitive(geoIdval);
    
                    for (int idx = 0; idx < hits->size(); idx++) {
                        const Acts::Vector3 trueHitGlob = 
                            {hits->at(idx).X() * Acts::UnitConstants::mm, 
                            hits->at(idx).Y() * Acts::UnitConstants::mm, 
                            hits->at(idx).Z() * Acts::UnitConstants::mm};
        
                        const Acts::Vector2 trueHitLoc = 
                            convertToLoc(trueHitGlob, geoId, m_cfg.gOpt);
        
                        Acts::BoundVector parameters = Acts::BoundVector::Zero();
                        parameters[Acts::eBoundLoc0] = trueHitLoc[Acts::eBoundLoc0];
                        parameters[Acts::eBoundLoc1] = trueHitLoc[Acts::eBoundLoc1];
        
                        const Acts::Vector3 trueP = 
                            {dirs->at(idx).Px(), dirs->at(idx).Py(), dirs->at(idx).Pz()};
        
                        const Acts::Vector3 direction = trueP.normalized();
        
                        parameters[Acts::eBoundPhi] = Acts::VectorHelpers::phi(direction);
                        parameters[Acts::eBoundTheta] = Acts::VectorHelpers::theta(direction);
                        parameters[Acts::eBoundQOverP] = 
                            1/(dirs->at(idx).E() * Acts::UnitConstants::GeV);
                        parameters[Acts::eBoundTime] = dirs->at(idx).T();
    
                        Acts::Vector2 stddev(5  * Acts::UnitConstants::um,
                            5  * Acts::UnitConstants::um);
                        Acts::SquareMatrix2 cov = stddev.cwiseProduct(stddev).asDiagonal();
        
                        SimpleSourceLink ssl(trueHitLoc, cov, geoId, eventId);
                        Acts::SourceLink sl{ssl};
    
                        LUXEMeasurement::SimMeasurement measurement{
                            .sourceLink = sl,
                            .truthParameters = parameters,
                            .trackId = trackId->at(idx)
                        };

                        measurements->push_back(measurement);
                    }
        };
};

auto defaultSimConfig() {
    LUXEROOTSimDataReader::Config config;
    config.treeName = "clusters";
    config.vector3Keys = {"tru_hit"};
    config.lorentzKeys = {"tru_p"};
    config.vectorIntKeys = {"tru_trackId"};
    config.intKeys = {"eventId", "geoId"};
    return config;
}

}  // namespace LUXEROOTReader
