#include "ActsLUXEPipeline/LUXEBinnedMagneticField.hpp"

namespace LUXEMagneticField {

Acts::InterpolatedBFieldMap<vGrid> buildBinnedBField(
    const Acts::MagneticFieldProvider& mFieldVal,
    const posTransform& transformPos,
    const fieldTransform& transformBField,
    const vGridOptions& gridOpt,
    const Acts::MagneticFieldContext& mctx) {
        vAxis x(gridOpt.xBins);
        vAxis y(gridOpt.yBins);
        vAxis z(gridOpt.zBins);

        auto mCache = mFieldVal.makeCache(mctx);

        // magnetic field known on grid in (x,y,z)
        vGrid g(std::make_tuple(std::move(x), std::move(y), std::move(z)));
    
        // set grid values
        for (std::size_t i = 1; i <= g.numLocalBins().at(0) + 1; i++) {
            for (std::size_t j = 1; j <= g.numLocalBins().at(1) + 1; j++) {
                for (std::size_t k = 1; k <= g.numLocalBins().at(2) + 1; k++) {
                    typename vGrid::index_t indices = {{i, j, k}};
                    const auto &llCorner = g.lowerLeftBinEdge(indices);
                    Acts::Vector3 pos(llCorner[0], llCorner[1], llCorner[2]);
                    g.atLocalBins(indices) = mFieldVal.getField(pos, mCache).value();
                }
            }
        }
    
        // create BField service
        Acts::InterpolatedBFieldMap<vGrid> 
            bField{{transformPos, transformBField, std::move(g)}};
    
        return bField;
}

Acts::InterpolatedBFieldMap<eGrid> buildBinnedBField(
    const Acts::MagneticFieldProvider& mFieldVal,
    const posTransform& transformPos,
    const fieldTransform& transformBField,
    const eGridOptions& gridOpt,
    const Acts::MagneticFieldContext& mctx) {
        auto [xMin, xMax, xBins] = gridOpt.xBins;
        auto [yMin, yMax, yBins] = gridOpt.yBins;
        auto [zMin, zMax, zBins] = gridOpt.zBins;

        eAxis x(xMin, xMax, xBins);
        eAxis y(yMin, yMax, yBins);
        eAxis z(zMin, zMax, zBins);

        auto mCache = mFieldVal.makeCache(mctx);

        // magnetic field known on grid in (x,y,z)
        eGrid g(std::make_tuple(std::move(x), std::move(y), std::move(z)));
    
        // set grid values
        for (std::size_t i = 1; i <= g.numLocalBins().at(0) + 1; i++) {
            for (std::size_t j = 1; j <= g.numLocalBins().at(1) + 1; j++) {
                for (std::size_t k = 1; k <= g.numLocalBins().at(2) + 1; k++) {
                    typename eGrid::index_t indices = {{i, j, k}};
                    const auto &llCorner = g.lowerLeftBinEdge(indices);
                    Acts::Vector3 pos(llCorner[0], llCorner[1], llCorner[2]);
                    g.atLocalBins(indices) = mFieldVal.getField(pos, mCache).value();
                }
            }
        }
    
        // create BField service
        Acts::InterpolatedBFieldMap<eGrid> 
            bField{{transformPos, transformBField, std::move(g)}};
    
        return bField;
}

} // namespace LUXEMagneticField
