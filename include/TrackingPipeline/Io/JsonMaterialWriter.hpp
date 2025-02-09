#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <map>
#include <memory>
#include <string>
#include <utility>

#include "TrackingPipeline/Material/IMaterialWriter.hpp"

namespace Acts {

class ISurfaceMaterial;
class IVolumeMaterial;

using SurfaceMaterialMap =
    std::map<GeometryIdentifier, std::shared_ptr<const ISurfaceMaterial>>;

using VolumeMaterialMap =
    std::map<GeometryIdentifier, std::shared_ptr<const IVolumeMaterial>>;

using DetectorMaterialMaps = std::pair<SurfaceMaterialMap, VolumeMaterialMap>;

}  // namespace Acts

/// @class Json Material writer
///
/// @brief Writes out Detector material maps
/// using the Json Geometry converter
class JsonMaterialWriter : public IMaterialWriter {
 public:
  struct Config {
    /// The config class of the converter
    Acts::MaterialMapJsonConverter::Config converterCfg;
    /// Output file name
    std::string fileName = "material";
  };

  /// Constructor
  ///
  /// @param config The configuration struct of the writer
  /// @param level The log level
  JsonMaterialWriter(const Config& config, Acts::Logging::Level level);

  /// Virtual destructor
  ~JsonMaterialWriter() override;

  /// Write out the material map
  ///
  /// @param detMaterial is the SurfaceMaterial and VolumeMaterial maps
  void writeMaterial(const Acts::DetectorMaterialMaps& detMaterial) override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  const Acts::Logger& logger() const { return *m_logger; }

  /// The logger instance
  std::unique_ptr<const Acts::Logger> m_logger{nullptr};

  /// The config of the writer
  Config m_cfg;

  /// The material converter
  std::unique_ptr<Acts::MaterialMapJsonConverter> m_converter{nullptr};
};
