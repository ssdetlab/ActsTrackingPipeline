
#include "TrackingPipeline/EventData/SimpleSourceLink.hpp"
#include "TrackingPipeline/TrackFitting/FastGX2Fiter.hpp"

using namespace Acts::UnitLiterals;

FastGX2Fitter::FastGX2Fitter(const Config& cfg) : m_cfg(cfg) {
  m_firstLayerGeoId.setSensitive(m_cfg.firstLayerGeoId);
  m_lastLayerGeoId.setSensitive(m_cfg.lastLayerGeoId);

  Acts::GeometryContext gctx;

  m_primaryInterlayerDistance =
      (m_cfg.surfaceMap.at(m_lastLayerGeoId)->center(gctx)(m_cfg.primaryIdx) -
       m_cfg.surfaceMap.at(m_firstLayerGeoId)->center(gctx)(m_cfg.primaryIdx)) /
      (m_cfg.surfaceMap.size() - 1);
}

double FastGX2Fitter::gx2Fit(const Acts::GeometryContext& gctx,
                             const std::vector<Acts::SourceLink>& sourceLinks,
                             const std::vector<std::size_t>& sourceLinksIndices,
                             Acts::Vector3& pos, Acts::Vector3& dir,
                             Acts::ActsSquareMatrix<6>& cov) {
  Acts::ActsDynamicMatrix G =
      Acts::ActsDynamicMatrix::Zero(2 * sourceLinksIndices.size(), 4);

  Eigen::VectorXd X(2 * sourceLinksIndices.size());
  for (std::size_t i = 0; i < sourceLinksIndices.size(); i++) {
    std::size_t idx = sourceLinksIndices.at(i);
    const auto& ssl = sourceLinks.at(idx).get<SimpleSourceLink>();
    const Acts::Vector3& parameters = ssl.parametersGlob();
    const auto& geoId = ssl.geometryId();
    X(2 * i) = parameters(m_cfg.longIdx);
    X(2 * i + 1) = parameters(m_cfg.shortIdx);

    G(2 * i, 0) = 1;
    G(2 * i, 2) = m_primaryInterlayerDistance * i;
    G(2 * i + 1, 1) = 1;
    G(2 * i + 1, 3) = m_primaryInterlayerDistance * i;
  }
  Eigen::LDLT<Acts::ActsDynamicMatrix> ldltD(
      constructCov(sourceLinks, sourceLinksIndices, dir));
  Acts::SquareMatrix4 B = G.transpose() * ldltD.solve(G);
  Acts::Vector4 estimates = B.ldlt().solve(G.transpose() * ldltD.solve(X));

  pos(m_cfg.primaryIdx) =
      m_cfg.surfaceMap.at(m_firstLayerGeoId)->center(gctx)(m_cfg.primaryIdx);
  pos(m_cfg.longIdx) = estimates(0);
  pos(m_cfg.shortIdx) = estimates(1);

  double tLong = estimates(2);
  double tShort = estimates(3);
  dir(m_cfg.primaryIdx) = 1.0;
  dir(m_cfg.longIdx) = tLong;
  dir(m_cfg.shortIdx) = tShort;
  dir.normalize();

  double denom = std::pow(1 + tLong * tLong + tShort * tShort, 1.5);

  Acts::ActsMatrix<6, 4> jacToGlob;
  jacToGlob.setZero();
  jacToGlob(1, 0) = 1;
  jacToGlob(2, 1) = 1;

  jacToGlob(3, 2) = -tLong / denom;
  jacToGlob(3, 3) = -tShort / denom;

  jacToGlob(4, 2) = (1 + tShort * tShort) / denom;
  jacToGlob(4, 3) = -tLong * tShort / denom;

  jacToGlob(5, 2) = -tLong * tShort / denom;
  jacToGlob(5, 3) = (1 + tLong * tLong) / denom;

  cov = jacToGlob * B.inverse() * jacToGlob.transpose();
  return (X - G * estimates).transpose() * ldltD.solve(X - G * estimates);
}

Acts::ActsDynamicMatrix FastGX2Fitter::constructCov(
    const std::vector<Acts::SourceLink>& sourceLinks,
    const std::vector<std::size_t>& sourceLinksIndices,
    const Acts::Vector3& dir) {
  Acts::ActsDynamicMatrix D = Acts::ActsDynamicMatrix::Zero(
      2 * sourceLinksIndices.size(), 2 * sourceLinksIndices.size());
  double mcsVarFactor =
      std::pow(m_primaryInterlayerDistance * dir.x() * m_cfg.thetaMcpRms, 2);
  for (std::size_t i = 0; i < sourceLinksIndices.size(); i++) {
    std::size_t idx = sourceLinksIndices.at(i);
    const auto& ssl = sourceLinks.at(idx).get<SimpleSourceLink>();

    double mcpVar = mcsVarFactor * i * (1 + 3 * i + 2 * i * i) / 6.0;
    Acts::SquareMatrix2 varMat =
        ssl.covariance() + Acts::SquareMatrix2::Identity() * mcpVar;
    D.block(i * 2, i * 2, 2, 2) = varMat;
    for (std::size_t j = 0; j < i; j++) {
      double covPiPj = mcsVarFactor * j * (j + 1) * (1 - j + 3 * i) / 6.0;
      Acts::SquareMatrix2 covMat = Acts::SquareMatrix2::Identity() * covPiPj;
      D.block(i * 2, j * 2, 2, 2) = covMat;
      D.block(j * 2, i * 2, 2, 2) = covMat;
    }
  }
  return std::move(D);
}
