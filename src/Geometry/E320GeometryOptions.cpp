#include "TrackingPipeline/Geometry/E320GeometryOptions.hpp"

#include "TrackingPipeline/Geometry/detail/BinningValueUtils.hpp"
#include "toml++/impl/parser.hpp"
#include "toml++/toml.hpp"

E320::GeometryOptions::GeometryOptions() {
  const std::string pathToCfg =
      "/home/romanurmanov/work/TrackingPipeline/ActsTrackingPipeline/conf/"
      "geometry/"
      "E320PrototypeGeometryConfig.toml";
  auto cfg = toml::parse_file(pathToCfg);

  auto getEntryDouble = [&cfg](const std::string& section,
                               const std::string& subsection) {
    return cfg[section][subsection].value<double>().value();
  };
  auto getEntrySizeT = [&cfg](const std::string& section,
                              const std::string& subsection) {
    return cfg[section][subsection].value<std::size_t>().value();
  };
  auto getEntryInt = [&cfg](const std::string& section,
                            const std::string& subsection) {
    return cfg[section][subsection].value<int>().value();
  };

  // --------------------------------------------------------------
  // General parameters

  // Detector binning directions
  primaryBinValue =
      detail::indexToBinningValue(getEntrySizeT("General", "primaryBinValue"));
  longBinValue =
      detail::indexToBinningValue(getEntrySizeT("General", "longBinValue"));
  shortBinValue =
      detail::indexToBinningValue(getEntrySizeT("General", "shortBinValue"));

  primaryIdx = detail::binningValueToIndex(primaryBinValue);
  longIdx = detail::binningValueToIndex(longBinValue);
  shortIdx = detail::binningValueToIndex(shortBinValue);

  primaryDir = detail::binningValueToDirection(primaryBinValue);
  longDir = detail::binningValueToDirection(longBinValue);
  shortDir = detail::binningValueToDirection(shortBinValue);

  // Volume geometry ID prefactors
  magVolumeIdPrefactor = getEntrySizeT("General", "magVolumeIdPrefactor");
  gapVolumeIdPrefactor = getEntrySizeT("General", "gapVolumeIdPrefactor");

  worldHalfLong = getEntryDouble("General", "worldHalfLong") * 1_m;
  worldHalfShort = getEntryDouble("General", "worldHalfShort") * 1_m;

  // Rotation of the plane surfaces into the global frame
  toWorldAngleX = M_PI_2 * getEntryDouble("General", "toWorldAngleX");
  toWorldAngleY = M_PI_2 * getEntryDouble("General", "toWorldAngleY");
  toWorldAngleZ = M_PI_2 * getEntryDouble("General", "toWorldAngleZ");

  // --------------------------------------------------------------
  // Parameters of the BPMs

  bpmHalfX = getEntryDouble("BPMs", "bpmHalfX") * 1_mm;
  bpmHalfY = getEntryDouble("BPMs", "bpmHalfY") * 1_mm;
  bpmThickness = getEntryDouble("BPMs", "bpmThickness") * 1_mm;

  bpm0CenterPrimary = getEntryDouble("BPMs", "bpm0CenterPrimary") * 1_mm;
  bpm1CenterPrimary = getEntryDouble("BPMs", "bpm1CenterPrimary") * 1_mm;
  bpm2CenterPrimary = getEntryDouble("BPMs", "bpm2CenterPrimary") * 1_mm;
  bpm3CenterPrimary = getEntryDouble("BPMs", "bpm3CenterPrimary") * 1_mm;

  bpm0CenterLong = getEntryDouble("BPMs", "bpm0CenterLong") * 1_mm;
  bpm1CenterLong = getEntryDouble("BPMs", "bpm1CenterLong") * 1_mm;
  bpm2CenterLong = getEntryDouble("BPMs", "bpm2CenterLong") * 1_mm;
  bpm3CenterLong = getEntryDouble("BPMs", "bpm3CenterLong") * 1_mm;

  bpm0CenterShort = getEntryDouble("BPMs", "bpm0CenterShort") * 1_mm;
  bpm1CenterShort = getEntryDouble("BPMs", "bpm1CenterShort") * 1_mm;
  bpm2CenterShort = getEntryDouble("BPMs", "bpm2CenterShort") * 1_mm;
  bpm3CenterShort = getEntryDouble("BPMs", "bpm3CenterShort") * 1_mm;

  int bpm0GeometryId = getEntryInt("BPMs", "bpm0GeometryId");
  int bpm1GeometryId = getEntryInt("BPMs", "bpm1GeometryId");
  int bpm2GeometryId = getEntryInt("BPMs", "bpm2GeometryId");
  int bpm3GeometryId = getEntryInt("BPMs", "bpm3GeometryId");

  bpm0Parameters =
      SurfaceParameters{{primaryBinValue, bpm0CenterPrimary, toWorldAngleX},
                        {longBinValue, bpm0CenterLong, toWorldAngleY},
                        {shortBinValue, bpm0CenterShort, toWorldAngleZ},
                        bpm0GeometryId};
  bpm1Parameters =
      SurfaceParameters{{primaryBinValue, bpm1CenterPrimary, toWorldAngleX},
                        {longBinValue, bpm1CenterLong, toWorldAngleY},
                        {shortBinValue, bpm1CenterShort, toWorldAngleZ},
                        bpm1GeometryId};
  bpm2Parameters =
      SurfaceParameters{{primaryBinValue, bpm2CenterPrimary, toWorldAngleX},
                        {longBinValue, bpm2CenterLong, toWorldAngleY},
                        {shortBinValue, bpm2CenterShort, toWorldAngleZ},
                        bpm2GeometryId};
  bpm3Parameters =
      SurfaceParameters{{primaryBinValue, bpm3CenterPrimary, toWorldAngleX},
                        {longBinValue, bpm3CenterLong, toWorldAngleY},
                        {shortBinValue, bpm3CenterShort, toWorldAngleZ},
                        bpm3GeometryId};

  // --------------------------------------------------------------
  // Parameters of the Be window

  beWindowHalfX = getEntryDouble("BeWindow", "beWindowHalfX") * 1_mm;
  beWindowHalfY = getEntryDouble("BeWindow", "beWindowHalfY") * 1_mm;
  beWindowThickness = getEntryDouble("BeWindow", "beWindowThickness") * 1_mm;
  beWindowCenterPrimary =
      getEntryDouble("BeWindow", "beWindowCenterPrimary") * 1_mm;
  beWindowCenterLong = getEntryDouble("BeWindow", "beWindowCenterLong") * 1_mm;
  beWindowCenterShort =
      getEntryDouble("BeWindow", "beWindowCenterShort") * 1_mm;
  int beWindowGeometryId = getEntryInt("BeWindow", "beWindowGeometryId");

  beWindowParameters =
      SurfaceParameters{{primaryBinValue, beWindowCenterPrimary, toWorldAngleX},
                        {longBinValue, beWindowCenterLong, toWorldAngleY},
                        {shortBinValue, beWindowCenterShort, toWorldAngleZ},
                        beWindowGeometryId};

  // --------------------------------------------------------------
  // Parameters of the IP surface

  ipSurfaceHalfX = getEntryDouble("IPSurface", "ipSurfaceHalfX") * 1_mm;
  ipSurfaceHalfY = getEntryDouble("IPSurface", "ipSurfaceHalfY") * 1_mm;
  ipSurfaceThickness = getEntryDouble("IPSurface", "ipSurfaceThickness") * 1_mm;
  ipSurfaceCenterPrimary =
      getEntryDouble("IPSurface", "ipSurfaceCenterPrimary") * 1_mm;
  ipSurfaceCenterLong =
      getEntryDouble("IPSurface", "ipSurfaceCenterLong") * 1_mm;
  ipSurfaceCenterShort =
      getEntryDouble("IPSurface", "ipSurfaceCenterShort") * 1_mm;
  int ipSurfaceGeometryId = getEntryInt("IPSurface", "ipSurfaceGeometryId");

  ipSurfaceParameters = SurfaceParameters{
      {primaryBinValue, ipSurfaceCenterPrimary, toWorldAngleX},
      {longBinValue, ipSurfaceCenterLong, toWorldAngleY},
      {shortBinValue, ipSurfaceCenterShort, toWorldAngleZ},
      ipSurfaceGeometryId};

  // --------------------------------------------------------------
  // Parameters of the PDC window

  pdcWindowHalfX = getEntryDouble("PDCWindowSurface", "pdcWindowHalfX") * 1_mm;
  pdcWindowHalfY = getEntryDouble("PDCWindowSurface", "pdcWindowHalfY") * 1_mm;

  pdcWindowThickness =
      getEntryDouble("PDCWindowSurface", "pdcWindowThickness") * 1_mm;

  std::size_t pdcWindowMaterialNBinsX =
      getEntrySizeT("PDCWindowSurface", "pdcWindowMaterialBinningX");
  pdcWindowMaterialBinningX =
      Acts::BinUtility(pdcWindowMaterialNBinsX, -pdcWindowHalfX, pdcWindowHalfX,
                       Acts::closed, Acts::BinningValue::binX);
  std::size_t pdcWindowMaterialNBinsY =
      getEntrySizeT("PDCWindowSurface", "pdcWindowMaterialBinningY");
  pdcWindowMaterialBinningY =
      Acts::BinUtility(pdcWindowMaterialNBinsY, -pdcWindowHalfY, pdcWindowHalfY,
                       Acts::closed, Acts::BinningValue::binY);

  pdcWindowCenterPrimary =
      getEntryDouble("PDCWindowSurface", "pdcWindowCenterPrimary") * 1_mm;
  pdcWindowCenterLong =
      getEntryDouble("PDCWindowSurface", "pdcWindowCenterLong") * 1_mm;
  pdcWindowCenterShort =
      getEntryDouble("PDCWindowSurface", "pdcWindowCenterShort") * 1_mm;
  int pdcWindowGeometryId =
      getEntryInt("PDCWindowSurface", "pdcWindowGeometryId");

  pdcWindowParameters = SurfaceParameters{
      {primaryBinValue, pdcWindowCenterPrimary, toWorldAngleX},
      {longBinValue, pdcWindowCenterLong, toWorldAngleY},
      {shortBinValue, pdcWindowCenterShort, toWorldAngleZ},
      pdcWindowGeometryId};

  // --------------------------------------------------------------
  // Parameters of the tracking chambers

  // Chip size in chips's local coordinates
  chipHalfX = getEntryDouble("TrackingChamber", "chipHalfX") * 1_mm;
  chipHalfY = getEntryDouble("TrackingChamber", "chipHalfY") * 1_mm;

  pixelHalfX = getEntryDouble("TrackingChamber", "pixelHalfX") * 1_um;
  pixelHalfY = getEntryDouble("TrackingChamber", "pixelHalfY") * 1_um;

  pixelThickness = getEntryDouble("TrackingChamber", "pixelThickness") * 1_um;

  // Allpix2 measurement resolutions
  const auto* allpixResolutions =
      cfg["TrackingChamber"]["allpixResolution"].as_array();
  for (auto it = allpixResolutions->begin(); it != allpixResolutions->end();
       it++) {
    const auto& entry = *it->as_table();
    allpixErrors.insert(
        {entry["clusterSize"].value<std::size_t>().value(),
         {
             entry["resolutionX"].value<double>().value() * 1_um,
             entry["resolutionY"].value<double>().value() * 1_um,
         }});
  }

  // Volume spacing around the chips
  chipVolumeHalfSpacing =
      getEntryDouble("TrackingChamber", "chipVolumeHalfSpacing") * 1_mm;

  // Distance between the chips
  interChipDistance =
      cfg["TrackingChamber"]["interChipDistance"].value<double>().value() *
      1_mm;
  tcWindowToFirstChipDistance =
      getEntryDouble("TrackingChamber", "tcWindowToFirstChipDistance") * 1_mm;
  tcWindowToLastChipDistance =
      getEntryDouble("TrackingChamber", "tcWindowToLastChipDistance") * 1_mm;

  // Transverse volume parameters
  tcHalfLong = chipHalfX;
  tcHalfShort = chipHalfY;

  std::size_t chipMaterialNBinsX =
      getEntrySizeT("TrackingChamber", "chipMaterialBinningX");
  chipMaterialBinningX =
      Acts::BinUtility(chipMaterialNBinsX, -chipHalfX, chipHalfX, Acts::closed,
                       Acts::BinningValue::binX);
  std::size_t chipMaterialNBinsY =
      getEntrySizeT("TrackingChamber", "chipMaterialBinningY");
  chipMaterialBinningY =
      Acts::BinUtility(chipMaterialNBinsY, -chipHalfY, chipHalfY, Acts::closed,
                       Acts::BinningValue::binY);

  // Placement
  ipTcDistance = getEntryDouble("TrackingChamber", "ipTcDistance") * 1_mm;
  tcCenterLong = getEntryDouble("TrackingChamber", "tcCenterLong") * 1_mm;
  tcCenterShort = getEntryDouble("TrackingChamber", "tcCenterShort") * 1_mm;

  int chip0GeometryId = getEntryInt("TrackingChamber", "chip0GeometryId");
  int chip1GeometryId = getEntryInt("TrackingChamber", "chip1GeometryId");
  int chip2GeometryId = getEntryInt("TrackingChamber", "chip2GeometryId");
  int chip3GeometryId = getEntryInt("TrackingChamber", "chip3GeometryId");
  int chip4GeometryId = getEntryInt("TrackingChamber", "chip4GeometryId");

  tcParameters = std::vector<SurfaceParameters>{
      SurfaceParameters{{primaryBinValue, ipTcDistance + 0 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tcCenterLong, toWorldAngleY},
                        {shortBinValue, tcCenterShort, toWorldAngleZ},
                        chip0GeometryId},
      SurfaceParameters{{primaryBinValue, ipTcDistance + 1 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tcCenterLong, toWorldAngleY},
                        {shortBinValue, tcCenterShort, toWorldAngleZ},
                        chip1GeometryId},
      SurfaceParameters{{primaryBinValue, ipTcDistance + 2 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tcCenterLong, toWorldAngleY},
                        {shortBinValue, tcCenterShort, toWorldAngleZ},
                        chip2GeometryId},
      SurfaceParameters{{primaryBinValue, ipTcDistance + 3 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tcCenterLong, toWorldAngleY},
                        {shortBinValue, tcCenterShort, toWorldAngleZ},
                        chip3GeometryId},
      SurfaceParameters{{primaryBinValue, ipTcDistance + 4 * interChipDistance,
                         toWorldAngleX},
                        {longBinValue, tcCenterLong, toWorldAngleY},
                        {shortBinValue, tcCenterShort, toWorldAngleZ},
                        chip4GeometryId}};

  tcHalfPrimary = interChipDistance * (tcParameters.size() - 1) / 2.0;
  tcCenterPrimary = ipTcDistance + tcHalfPrimary;

  // --------------------------------------------------------------
  // Parameters of the composite field

  compositeMagFieldId = getEntrySizeT("Magnets", "compositeMagFieldId");

  // --------------------------------------------------------------
  // Parameters of the quads

  quad1Id = getEntrySizeT("Magnets", "quad1Id");
  quad1HalfPrimary = getEntryDouble("Magnets", "quad1HalfPrimary") * 1_mm;
  quad1HalfLong = getEntryDouble("Magnets", "quad1HalfLong") * 1_mm;
  quad1HalfShort = getEntryDouble("Magnets", "quad1HalfShort") * 1_mm;
  quad1Gradient = getEntryDouble("Magnets", "quad1Gradient") * 1_T / 1_m;

  quad2Id = getEntrySizeT("Magnets", "quad2Id");
  quad2HalfPrimary = getEntryDouble("Magnets", "quad2HalfPrimary") * 1_mm;
  quad2HalfLong = getEntryDouble("Magnets", "quad2HalfLong") * 1_mm;
  quad2HalfShort = getEntryDouble("Magnets", "quad2HalfShort") * 1_mm;
  quad2Gradient = getEntryDouble("Magnets", "quad2Gradient") * 1_T / 1_m;

  quad3Id = getEntrySizeT("Magnets", "quad3Id");
  quad3HalfPrimary = getEntryDouble("Magnets", "quad3HalfPrimary") * 1_mm;
  quad3HalfLong = getEntryDouble("Magnets", "quad3HalfLong") * 1_mm;
  quad3HalfShort = getEntryDouble("Magnets", "quad3HalfShort") * 1_mm;
  quad3Gradient = getEntryDouble("Magnets", "quad3Gradient") * 1_T / 1_m;

  quad0CenterPrimary = getEntryDouble("Magnets", "quad0CenterPrimary") * 1_mm;
  quad0CenterLong = getEntryDouble("Magnets", "quad0CenterLong") * 1_mm;
  quad0CenterShort = getEntryDouble("Magnets", "quad0CenterShort") * 1_mm;

  quad1CenterPrimary = getEntryDouble("Magnets", "quad1CenterPrimary") * 1_mm;
  quad1CenterLong = getEntryDouble("Magnets", "quad1CenterLong") * 1_mm;
  quad1CenterShort = getEntryDouble("Magnets", "quad1CenterShort") * 1_mm;

  quad2CenterPrimary = getEntryDouble("Magnets", "quad2CenterPrimary") * 1_mm;
  quad2CenterLong = getEntryDouble("Magnets", "quad2CenterLong") * 1_mm;
  quad2CenterShort = getEntryDouble("Magnets", "quad2CenterShort") * 1_mm;

  quad3CenterPrimary = getEntryDouble("Magnets", "quad3CenterPrimary") * 1_mm;
  quad3CenterLong = getEntryDouble("Magnets", "quad3CenterLong") * 1_mm;
  quad3CenterShort = getEntryDouble("Magnets", "quad3CenterShort") * 1_mm;

  // --------------------------------------------------------------
  // Parameters of the x-corrector

  xCorrectorId = getEntrySizeT("Magnets", "xCorrectorId");
  xCorrectorHalfPrimary =
      getEntryDouble("Magnets", "xCorrectorHalfPrimary") * 1_mm;
  xCorrectorHalfLong = getEntryDouble("Magnets", "xCorrectorHalfLong") * 1_mm;
  xCorrectorHalfShort = getEntryDouble("Magnets", "xCorrectorHalfShort") * 1_mm;
  xCorrectorFieldStrength =
      getEntryDouble("Magnets", "xCorrectorFieldStrength") * 1_T;

  xCorrectorDirIdx = longIdx;

  xCorrectorFieldPrimary = 0;
  xCorrectorFieldLong = xCorrectorFieldStrength;
  xCorrectorFieldShort = 0;

  xCorrectorCenterPrimary =
      getEntryDouble("Magnets", "xCorrectorCenterPrimary") * 1_mm;
  xCorrectorCenterLong =
      getEntryDouble("Magnets", "xCorrectorCenterLong") * 1_mm;
  xCorrectorCenterShort =
      getEntryDouble("Magnets", "xCorrectorCenterShort") * 1_mm;

  // --------------------------------------------------------------
  // Parameters of the dipole

  dipoleId = getEntrySizeT("Magnets", "dipoleId");
  dipoleHalfPrimary = getEntryDouble("Magnets", "dipoleHalfPrimary") * 1_mm;
  dipoleHalfLong = getEntryDouble("Magnets", "dipoleHalfLong") * 1_mm;
  dipoleHalfShort = getEntryDouble("Magnets", "dipoleHalfShort") * 1_mm;
  dipoleFieldStrength = getEntryDouble("Magnets", "dipoleFieldStrength") * 1_T;

  dipoleDirIdx = shortIdx;

  dipoleFieldPrimary = 0;
  dipoleFieldLong = 0;
  dipoleFieldShort = dipoleFieldStrength;

  dipoleCenterPrimary = getEntryDouble("Magnets", "dipoleCenterPrimary") * 1_mm;
  dipoleCenterLong = getEntryDouble("Magnets", "dipoleCenterLong") * 1_mm;
  dipoleCenterShort = getEntryDouble("Magnets", "dipoleCenterShort") * 1_mm;
}
