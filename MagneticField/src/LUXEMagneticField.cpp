#include "ActsLUXEPipeline/LUXEMagneticField.hpp"

namespace LUXEMagneticField {

BField_t buildLUXEBField() {
    Acts::MagneticFieldContext mfContext = Acts::MagneticFieldContext();

    // definition of BField
    struct BField {
        static Acts::Vector3 value(const std::array<double, 3> &v) {
            double x = v.at(0);
            double y = v.at(1);
            double z = v.at(2);
            double r = std::sqrt(std::pow(x,2) + std::pow(y,2) + std::pow(z,2));
            double r5 = std::pow(r,5);
            // linear in r and z so interpolation should be exact
            return Acts::Vector3(3 * x * z / r5, 3 * y * z / r5,
                                 (3 * std::pow(z, 2) - std::pow(r, 2)) / r5);
        }
    };

    // map (x,y,z) -> (x,y,z)
    auto transformPos =[](const Acts::Vector3& pos) {
        return pos;
    };

    // map (Bx,By,Bz) -> (Bx,By,Bz)
    auto transformBField =[](const Acts::Vector3& field, const Acts::Vector3&) {
        return field;
    };

    // magnetic field known on grid in (x,y,z)
    Acts::detail::EquidistantAxis x(0.0, 5.0, 5u);
    Acts::detail::EquidistantAxis y(0.0, 5.0, 5u);
    Acts::detail::EquidistantAxis z(0.0, 5.0, 5u);

    Grid_t g(std::make_tuple(std::move(x), std::move(y), std::move(z)));

    // set grid values
    for (std::size_t i = 1; i <= g.numLocalBins().at(0) + 1; ++i) {
        for (std::size_t j = 1; j <= g.numLocalBins().at(1) + 1; ++j) {
            for (std::size_t k = 1; k <= g.numLocalBins().at(2) + 1; ++k) {
                Grid_t::index_t indices = {{i, j, k}};
                const auto &llCorner = g.lowerLeftBinEdge(indices);
                g.atLocalBins(indices) = BField::value(llCorner);
            }
        }
    }

    // create BField service
    BField_t bField{{ transformPos, transformBField, std::move(g) }};

    return bField;
}

} // namespace LUXEMagneticField
