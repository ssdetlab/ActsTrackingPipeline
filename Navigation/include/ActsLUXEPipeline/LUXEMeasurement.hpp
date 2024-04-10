#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
#include "ActsLUXEPipeline/LUXEDataContainers.hpp"
#include "ActsLUXEPipeline/LUXEEffectiveMaterial.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"

#include <vector>

#include <memory>
#include <iostream>
#include <fstream>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

/// Propagator action to create measurements.
namespace LUXENavigator {

// Store only essential components
    struct MeasurementComponents {
        int eId;
        int gId;
        double l0;
        double l1;
        double phi;
        double theta;
        double qOverP;
        double time;
        double trueVertexX;
        double trueVertexY;
        double trueVertexZ;
        double trueVertexT;

        MeasurementComponents() :
                eId(),
                gId(),
                l0(),
                l1(),
                phi(),
                theta(),
                qOverP(),
                time(),
                trueVertexX(),
                trueVertexY(),
                trueVertexZ(),
                trueVertexT() {}

        MeasurementComponents(LUXEDataContainer::SimMeasurement m) :
                eId(m.sourceLink.get<SimpleSourceLink>().eventId),
                gId(m.sourceLink.get<SimpleSourceLink>().geometryId().sensitive()),
                l0(m.sourceLink.get<SimpleSourceLink>().parameters[0]),
                l1(m.sourceLink.get<SimpleSourceLink>().parameters[1]),
                phi(m.truthParameters[Acts::eBoundPhi]),
                theta(m.truthParameters[Acts::eBoundTheta]),
                qOverP(m.truthParameters[Acts::eBoundQOverP]),
                time(m.truthParameters[Acts::eBoundTime]),
                trueVertexX(m.trueVertex[0]),
                trueVertexY(m.trueVertex[1]),
                trueVertexZ(m.trueVertex[2]),
                trueVertexT(m.trueVertex[3]){}


        template<class Archive>
        void serialize(Archive& ar, const unsigned int version) {
            ar & eId;
            ar & gId;
            ar & l0;
            ar & l1;
            ar & phi;
            ar & theta;
            ar & qOverP;
            ar & time;
            ar & trueVertexX;
            ar & trueVertexY;
            ar & trueVertexZ;
            ar & trueVertexT;
        }
    };

    void saveMeasurementsToFile(const LUXEDataContainer::SimMeasurements& measurements, const std::string& filename) {
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }
        boost::archive::binary_oarchive oa(ofs);
        size_t size=measurements.size();
        oa << size; // Write the size of the vector

        for (const auto& measurement : measurements) {
            MeasurementComponents m_comp(measurement);
            oa << m_comp; // Serialize each Measurements object
        }
        ofs.flush();
    }

// Function to load vector of Measurements from a file
    LUXEDataContainer::SimMeasurements loadMeasurementsFromFile(const std::string& filename) {
        LUXEDataContainer::SimMeasurements measurements;
        std::ifstream ifs(filename, std::ios::binary);
        if (!ifs) {
            std::cerr << "Error opening file for reading: " << filename << std::endl;
            return measurements;
        }
        boost::archive::binary_iarchive ia(ifs);
        size_t vectorSize;
        ia >> vectorSize; // Read the size of the vector
        measurements.reserve(vectorSize);
        Acts::SquareMatrix2 covariance = Acts::SquareMatrix2::Identity();
        for (size_t i = 0; i < vectorSize; ++i) {
            MeasurementComponents mC;
            ia >> mC;
            Acts::ActsVector<2> params{mC.l0,mC.l1};
            Acts::GeometryIdentifier geometryId;
            geometryId.setSensitive(mC.gId);
            SimpleSourceLink ssl(params, covariance, geometryId, mC.eId);
            Acts::SourceLink sl{ssl};
            Acts::Vector4 trueVertex{mC.trueVertexX,mC.trueVertexY,
                                     mC.trueVertexZ,mC.trueVertexT};
            Acts::BoundVector truthParams;
            truthParams[Acts::eBoundLoc0] = mC.l0;
            truthParams[Acts::eBoundLoc1] = mC.l1;
            truthParams[Acts::eBoundPhi] = mC.phi;
            truthParams[Acts::eBoundTheta] = mC.theta;
            truthParams[Acts::eBoundQOverP] = mC.qOverP;
            truthParams[Acts::eBoundTime] = mC.time;
            LUXEDataContainer::SimMeasurement m{sl, truthParams, trueVertex, mC.eId};
            measurements.push_back(m);
        }
        std::cout<<"Loaded "<< vectorSize << " measurements from: "<<filename<<std::endl;
        return measurements;
    }

} // LUXENavigator