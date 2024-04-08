#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsLUXEPipeline/LUXESimpleSourceLink.hpp"
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

    enum class MeasurementType {
        eLoc0,
        eLoc1,
        eLoc01,
    };
    struct MeasurementResolution {
        MeasurementType type = MeasurementType::eLoc0;
        // Depending on the type, only the first value is used.
        std::array<double, 2> stddev = {15 * Acts::UnitConstants::um, 15 * Acts::UnitConstants::um};
    };

/// Measurement resolution configuration for a full detector geometry.
    using MeasurementResolutionMap =
            Acts::GeometryHierarchyMap<MeasurementResolution>;

/// Result struct for generated measurements and outliers.
    struct Measurement {
        unsigned int eventId;
        std::vector<SimpleSourceLink> sourceLinks;
        std::vector<Acts::Vector3> fullTrack;
        std::vector<Acts::BoundVector> truthParameters;
        std::vector<Acts::Vector3> globalPosition;
    };

// Store only essential components
    struct MeasurementComponents {
        int eId;
        int gId;
        double l0;
        double l1;
        double x;
        double y;
        double z;

        MeasurementComponents() :
                eId(),
                gId(),
                l0(),
                l1(),
                x(),
                y(),
                z() {}

        MeasurementComponents(Measurement m, int index) :
                eId(m.sourceLinks[index].eventId),
                gId(m.sourceLinks[index].geometryId().sensitive()),
                l0(m.sourceLinks[index].parameters[0]),
                l1(m.sourceLinks[index].parameters[1]),
                x(m.globalPosition[index][0]),
                y(m.globalPosition[index][1]),
                z(m.globalPosition[index][2]) {}


        template<class Archive>
        void serialize(Archive& ar, const unsigned int version) {
            ar & eId;
            ar & gId;
            ar & l0;
            ar & l1;
            ar & x;
            ar & y;
            ar & z;
        }
    };

    void saveMeasurementsToFile(const std::vector<Measurement>& measurements, const std::string& filename) {
        std::ofstream ofs(filename, std::ios::binary);
        if (!ofs) {
            std::cerr << "Error opening file for writing: " << filename << std::endl;
            return;
        }
        boost::archive::binary_oarchive oa(ofs);
        size_t size=0;
        for (auto m : measurements) {
            size+=m.sourceLinks.size();
        }
        oa << size; // Write the size of the vector

        for (const auto& measurement : measurements) {
            for (int i = 0; i<measurement.sourceLinks.size(); i++) {
                MeasurementComponents m_comp(measurement,i);
                oa << m_comp; // Serialize each Measurements object
            }
        }
        ofs.flush();
    }

// Function to load vector of Measurements from a file
    std::vector<Measurement> loadMeasurementsFromFile(const std::string& filename) {
        std::vector<Measurement> measurements;
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
            Measurement m;
            MeasurementComponents mC;
            ia >> mC;
            Acts::ActsVector<2> params{mC.l0,mC.l1};
            Acts::GeometryIdentifier geometryId;
            geometryId.setSensitive(mC.gId);
            SimpleSourceLink ssl(params, covariance, geometryId, mC.eId);
            m.sourceLinks = {ssl};
            m.globalPosition = {Acts::Vector3{mC.x,mC.y,mC.z}};
            measurements.push_back(m);
        }
        std::cout<<"Loaded "<< vectorSize << " measurements from: "<<filename<<std::endl;
        return measurements;
    }

} // LUXENavigator