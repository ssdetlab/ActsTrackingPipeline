#pragma once

#include "Eigen/Eigen"
#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Multivariate normal random variable generator
class NormalRandomVariable {
 public:
  NormalRandomVariable(Eigen::MatrixXd const& covar)
      : NormalRandomVariable(Eigen::VectorXd::Zero(covar.rows()), covar) {}

  NormalRandomVariable(Eigen::VectorXd const& mean,
                       Eigen::DiagonalWrapper<Eigen::MatrixXd> const& covar)
      : m_mean(mean), m_cov(covar) {
    Eigen::MatrixXd cov = covar.toDenseMatrix().matrix();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(cov);
    m_transform = eigenSolver.eigenvectors() *
                  eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  NormalRandomVariable(Eigen::VectorXd const& mean,
                       Eigen::MatrixXd const& covar)
      : m_mean(mean), m_cov(covar) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(covar);
    m_transform = eigenSolver.eigenvectors() *
                  eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
  }

  Eigen::VectorXd gen(RandomEngine& rng) const {
    static std::normal_distribution<> dist;

    return m_mean + m_transform * Eigen::VectorXd{m_mean.size()}.unaryExpr(
                                      [&](auto x) { return dist(rng); });
  }

  Eigen::VectorXd getMean() const { return m_mean; }

  Eigen::MatrixXd getCovariance() const { return m_cov; }

  Eigen::MatrixXd getTransform() const { return m_transform; }

 private:
  Eigen::VectorXd m_mean;
  Eigen::MatrixXd m_cov;
  Eigen::MatrixXd m_transform;
};
