#include "ActsLUXEPipeline/LUXEMagneticField.hpp"

namespace LUXEMagneticField {

BField_t buildLUXEBField(const transformationPos_t& transformPos,
                         const transformationBField_t& transformBField,
                         const bFieldValue_t& bFieldValue,
                         const std::vector<unsigned int> bins) {
    Acts::MagneticFieldContext mfContext = Acts::MagneticFieldContext();

    // magnetic field known on grid in (x,y,z)
    // TODO: make non-equidistant binning an option
    Acts::detail::EquidistantAxis x(0.0, 5.0, bins[0]);
    Acts::detail::EquidistantAxis y(0.0, 5.0, bins[1]);
    Acts::detail::EquidistantAxis z(0.0, 5.0, bins[2]);

    Grid_t g(std::make_tuple(std::move(x), std::move(y), std::move(z)));

    // set grid values
    for (std::size_t i = 1; i <= g.numLocalBins().at(0) + 1; ++i) {
        for (std::size_t j = 1; j <= g.numLocalBins().at(1) + 1; ++j) {
            for (std::size_t k = 1; k <= g.numLocalBins().at(2) + 1; ++k) {
                Grid_t::index_t indices = {{i, j, k}};
                const auto &llCorner = g.lowerLeftBinEdge(indices);
                g.atLocalBins(indices) = bFieldValue(llCorner);
            }
        }
    }

    // create BField service
    BField_t bField{{ transformPos, transformBField, std::move(g) }};

    return bField;
}

} // namespace LUXEPipeline
