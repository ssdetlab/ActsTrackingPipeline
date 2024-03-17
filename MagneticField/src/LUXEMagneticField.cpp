#include "ActsLUXEPipeline/LUXEMagneticField.hpp"

namespace LUXEMagneticField {

auto exampleDipole = [](const std::array<double, 3> &v) {
    double x = v.at(0);
    double y = v.at(1);
    double z = v.at(2);
    double r = std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
    double r5 = std::pow(r,5);

    return Acts::Vector3(-0.001 / r5, -0.001 , 0 * (std::pow(z, 2) - std::pow(r, 2)) / r5);
};

BField_t buildLUXEBField(const transformationPos_t& transformPos,
                         const transformationBField_t& transformBField,
                         const GridOptions gridOpt) {
    Acts::MagneticFieldContext mfContext = Acts::MagneticFieldContext();

    // magnetic field known on grid in (x,y,z)
    Acts::detail::EquidistantAxis x(gridOpt.limits[0].first, gridOpt.limits[0].second, gridOpt.bins[0]);
    Acts::detail::EquidistantAxis y(gridOpt.limits[1].first, gridOpt.limits[1].second, gridOpt.bins[1]);
    Acts::detail::EquidistantAxis z(gridOpt.limits[2].first, gridOpt.limits[2].second, gridOpt.bins[2]);

    Grid_t g(std::make_tuple(std::move(x), std::move(y), std::move(z)));

    // set grid values
    for (std::size_t i = 1; i <= g.numLocalBins().at(0) + 1; ++i) {
        for (std::size_t j = 1; j <= g.numLocalBins().at(1) + 1; ++j) {
            for (std::size_t k = 1; k <= g.numLocalBins().at(2) + 1; ++k) {
                Grid_t::index_t indices = {{i, j, k}};
                const auto &llCorner = g.lowerLeftBinEdge(indices);
                g.atLocalBins(indices) = LUXEMagneticField::exampleDipole(llCorner);
            }
        }
    }

    // create BField service
    BField_t bField{{ transformPos, transformBField, std::move(g) }};

    return bField;
}
} // namespace LUXEMagneticField
