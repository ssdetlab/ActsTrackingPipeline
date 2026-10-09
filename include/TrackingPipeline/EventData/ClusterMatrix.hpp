#pragma once

#include <cstddef>
#include <cstdint>
#include <stdexcept>

#include <Eigen/Eigen>

#include "TrackingPipeline/Utilities/IsAnyOfConcept.hpp"
#include "TrackingPipeline/Utilities/TypeList.hpp"

/// @brief Types list of allowed pixel position indices
using PixelIndexTypes = TypeList<std::uint8_t, std::uint16_t, std::uint16_t,
                                 std::uint64_t, std::size_t>;

/// @brief Class holding the matrix representation of a pixel cluster
///
/// @tparam T type indexing the pixel positon
template <typename T>
  requires IsAnyOf<T, PixelIndexTypes>
class ClusterMatrix {
 public:
  /// @brief Constructor
  ///
  /// @param pixels cluster pixels
  explicit ClusterMatrix(const std::vector<std::pair<T, T>> &pixels,
                         std::size_t matrixDim)
      : m_size(pixels.size()) {
    if (m_size > matrixDim) {
      throw std::runtime_error(
          "Cluster size is greater than allocated matrix dimensions");
    }
    m_matrix = Eigen::MatrixXi::Zero(matrixDim, matrixDim);

    // Find the top left pixel of the cluster
    std::size_t topLeftX = 1024;
    std::size_t topLeftY = 0;
    for (const auto &[hitX, hitY] : pixels) {
      if (topLeftX > hitX) {
        topLeftX = hitX;
      }
      if (topLeftY < hitY) {
        topLeftY = hitY;
      }
    }

    // Construct the pixel matrix
    for (const auto &[hitX, hitY] : pixels) {
      std::size_t col = hitX - topLeftX;
      std::size_t row = topLeftY - hitY;
      m_matrix(row, col) = 1;
    }

    // Calculate the X, Y extent
    m_lengthX = m_matrix.rowwise().sum().maxCoeff();
    m_lengthY = m_matrix.colwise().sum().maxCoeff();

    // Calculate shape ID
    m_shapeId = 0;
    for (std::size_t i = 0; i < m_matrix.size(); i++) {
      m_shapeId |= m_matrix(i) << i;
    }
  }

  /// @brief get Eigen representation of the pixel matrix
  const Eigen::MatrixXi &matrix() const { return m_matrix; }

  /// @brief get cluster size in number of pixels
  std::size_t size() const { return m_size; }

  /// @brief get cluster extent in local X direction
  std::size_t lengthX() const { return m_lengthX; }

  /// @brief get cluster extent in local Y direction
  std::size_t lengthY() const { return m_lengthY; }

  /// @brief get cluster unique shape identifier
  std::size_t shapeId() const { return m_shapeId; }

  /// @brief equality comparison operator
  ///
  /// @param other cluster matrix to compare to
  ///
  /// @return true if pixel matrices are identical, false otherwise
  bool operator==(const ClusterMatrix &other) const {
    return m_shapeId == other.shapeId();
  }

 private:
  /// Cluster shape identifier
  std::size_t m_shapeId;

  /// Cluster extent in local X direction
  std::size_t m_lengthX;

  /// Cluster extent in local Y direction
  std::size_t m_lengthY;

  /// Cluster size in number of pixels
  std::size_t m_size;

  /// Cluster pixel matrix
  Eigen::MatrixXi m_matrix;
};
