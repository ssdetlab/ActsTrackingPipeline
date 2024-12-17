#pragma once

#include "TrackingPipeline/Infrastructure/RandomNumbers.hpp"

#include "Eigen/Eigen"

/// @brief Multivariate normal random variable generator
struct NormalRandomVariable {
    NormalRandomVariable(Eigen::MatrixXd const& covar)
        : NormalRandomVariable(
            Eigen::VectorXd::Zero(covar.rows()), 
            covar) {}

    NormalRandomVariable(Eigen::VectorXd const& mean, Eigen::DiagonalWrapper<Eigen::MatrixXd> const& covar)
        : mean(mean) {
            Eigen::MatrixXd cov = covar.toDenseMatrix().matrix();

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> 
                eigenSolver(cov);
            transform = 
                eigenSolver.eigenvectors() * 
                    eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    NormalRandomVariable(Eigen::VectorXd const& mean, Eigen::MatrixXd const& covar)
        : mean(mean) {
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> 
                eigenSolver(covar);
            transform = 
                eigenSolver.eigenvectors() * 
                    eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();
    }

    Eigen::VectorXd mean;
    Eigen::MatrixXd transform;

    Eigen::VectorXd gen(RandomEngine& rng) const {
        static std::normal_distribution<> dist;

        return mean + 
            transform * Eigen::VectorXd{mean.size()}.unaryExpr(
                [&](auto x) { return dist(rng); });
    }
};
