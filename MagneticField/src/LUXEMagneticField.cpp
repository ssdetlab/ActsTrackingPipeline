#include "ActsLUXEPipeline/LUXEMagneticField.hpp"

namespace LUXEMagneticField {

auto exampleDipole = [](const std::array<double, 3> &v) {
    double x = v.at(0);
    double y = v.at(1);
    double z = v.at(2);
    double r = std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
    double r5 = std::pow(r,5);
    // linear in r and z so interpolation should be exact
    if (z<3962) {
        return Acts::Vector3(0,0,0);
    } else {
        return Acts::Vector3(0 * x * z / r5, -0.001 ,
                             (0 * std::pow(z, 2) - std::pow(r, 2)) / r5);
    }
};

BField_t buildLUXEBField(const transformationPos_t& transformPos,
                         const transformationBField_t& transformBField,
                         const std::vector<unsigned int> bins) {
    Acts::MagneticFieldContext mfContext = Acts::MagneticFieldContext();

    // magnetic field known on grid in (x,y,z)
    Acts::detail::EquidistantAxis x(-3000.0, 10000.0, bins[0]);
    Acts::detail::EquidistantAxis y(-3000.0, 10000.0, bins[1]);
    Acts::detail::EquidistantAxis z(-3000.0, 10000.0, bins[2]);

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
} // namespace LUXEPipeline
