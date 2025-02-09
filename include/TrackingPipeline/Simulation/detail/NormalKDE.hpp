#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "Eigen/Eigen"
#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

/// @brief Multivariate normal random variable generator
template <std::size_t Dim>
class NormalKDE {
 public:
  using Vector = Eigen::Vector<double, Dim>;

  NormalKDE(std::vector<Vector> sample, std::size_t nIterations,
            double sensitivity)
      : m_sample(sample) {
    Vector mu = Vector::Zero();
    Vector mu2 = Vector::Zero();

    for (const auto& point : sample) {
      mu += point;
      mu2 += point.cwiseProduct(point);
    }
    mu /= sample.size() - 1;
    mu2 /= sample.size() - 1;

    Vector scale = (mu2 - mu.cwiseProduct(mu)).cwiseSqrt();

    Vector iqr;
    for (std::size_t i = 0; i < Dim; i++) {
      std::sort(sample.begin(), sample.end(),
                [&](const auto& pointA, const auto& pointB) {
                  return pointA(i) < pointB(i);
                });

      double q1 = sample.at(sample.size() / 4)(i);
      double q3 = sample.at(3 * sample.size() / 4)(i);
      iqr(i) = q3 - q1;
    }

    Vector initialBw =
        0.9 * scale.cwiseMin(iqr / 1.34) * std::pow(sample.size(), -0.2);

    for (std::size_t j = 0; j < nIterations; j++) {
      initialBw = sheatherJonesBW(m_sample, initialBw);
    }

    if (sensitivity == 0) {
      m_bandwidths = {initialBw};
      return;
    }

    std::vector<Vector> pilotDensity;
    pilotDensity.reserve(m_sample.size());
    Vector invBw = initialBw.cwiseInverse();
    for (const auto& x : m_sample) {
      Vector res = Vector::Zero();
      for (const auto& point : sample) {
        res += normalKernel((x - point).cwiseProduct(invBw));
      }
      res = res.cwiseProduct(invBw) / m_sample.size();

      pilotDensity.push_back(res);
    }

    Vector logSum = Vector::Zero();
    for (const auto& pd : pilotDensity) {
      logSum += (pd.array() + 1e-10).log().matrix();
    }

    Vector geoMean = logSum / pilotDensity.size();
    geoMean = geoMean.array().exp();

    m_bandwidths.reserve(m_sample.size());
    geoMean = geoMean.cwiseInverse();
    for (auto& point : pilotDensity) {
      point = point.cwiseProduct(geoMean);
      point = point.array().pow(-sensitivity);
      m_bandwidths.push_back(initialBw.cwiseProduct(point));
    }
  }

  Vector sample(RandomEngine& rng) {
    std::uniform_int_distribution<> index(0, m_sample.size() - 1);

    Vector mean = m_sample.at(index(rng));
    Vector bw = (m_bandwidths.size() == 1) ? m_bandwidths.at(0)
                                           : m_bandwidths.at(index(rng));

    std::normal_distribution<> normal(0, 1);

    return mean + bw.cwiseProduct(Eigen::VectorXd{mean.size()}.unaryExpr(
                      [&](auto x) { return normal(rng); }));
  }

  // double eval(double x) {
  // double res = 0;
  // for (const auto& point : m_sample) {
  // res += normalKernel((x - point)/m_bandwidth);
  // }
  // res /= m_sample.size() * m_bandwidth;

  // return res;
  // }

  // double bandwidth() {
  // return m_bandwidth;
  // }

 private:
  std::vector<Vector> m_sample;
  std::vector<Vector> m_bandwidths;

  Vector normalKernel(const Vector& x) const {
    Vector prod = x.cwiseProduct(x);
    prod *= -0.5;
    prod = prod.array().exp();

    return 1. / std::sqrt(2 * M_PI) * prod;
  }

  Vector normalKernelSecondDerivative(const Vector& x) const {
    Vector prod = x.cwiseProduct(x);
    prod.array() -= 1;
    return prod.cwiseProduct(normalKernel(x));
  }

  Vector RSecondDerivative(const std::vector<Vector>& sample,
                           Vector pilotBandwidth) {
    Vector sum = Vector::Zero();
    for (const auto& point : sample) {
      auto u = pilotBandwidth.cwiseInverse().cwiseProduct(point);
      sum += normalKernelSecondDerivative(u).cwiseProduct(
          normalKernelSecondDerivative(u));
    }
    pilotBandwidth = pilotBandwidth.array().pow(5);
    pilotBandwidth *= sample.size();
    sum = sum.cwiseProduct(pilotBandwidth.cwiseInverse());
    return sum;
  }

  Vector sheatherJonesBW(const std::vector<Vector>& sample,
                         Vector pilotBandwidth) {
    const double RK = 1.0 / (2.0 * std::sqrt(M_PI));
    const double mu2K = 1.0;

    auto rSecondDerivative = RSecondDerivative(sample, pilotBandwidth);
    double n = sample.size();

    rSecondDerivative =
        RK * rSecondDerivative.cwiseInverse() / (n * mu2K * mu2K);
    rSecondDerivative = rSecondDerivative.array().pow(1.0 / 5.0);

    return rSecondDerivative;
  }
};
