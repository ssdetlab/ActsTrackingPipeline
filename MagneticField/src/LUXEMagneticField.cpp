#include "ActsLUXEPipeline/LUXEMagneticField.hpp"
#include "ActsLUXEPipeline/MagneticFields.hpp"

namespace LUXEMagneticField {
    
BField_t buildLUXEBField(const transformationPos_t& transformPos,
                         const transformationBField_t& transformBField,
                         const GridOptions gridOpt,
                         const IMagneticField& magneticField) {

    Acts::MagneticFieldContext mfContext = Acts::MagneticFieldContext();

    axis_t x(gridOpt.xBins);
    axis_t y(gridOpt.yBins);
    axis_t z(gridOpt.zBins);

    // magnetic field known on grid in (x,y,z)
    Grid_t g(std::make_tuple(std::move(x), std::move(y), std::move(z)));

    // set grid values
    for (std::size_t i = 1; i <= g.numLocalBins().at(0) + 1; i++) {
        for (std::size_t j = 1; j <= g.numLocalBins().at(1) + 1; j++) {
            for (std::size_t k = 1; k <= g.numLocalBins().at(2) + 1; k++) {
                Grid_t::index_t indices = {{i, j, k}};
                const auto &llCorner = g.lowerLeftBinEdge(indices);
                g.atLocalBins(indices) = magneticField.getValue(llCorner);
            }
        }
    }

    // create BField service
    BField_t bField{{ transformPos, transformBField, std::move(g) }};

    return bField;
}
} // namespace LUXEMagneticField
