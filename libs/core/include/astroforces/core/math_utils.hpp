/**
 * @file math_utils.hpp
 * @brief Shared vector/matrix math helpers.
 * @author Watosn
 */
#pragma once

#include <array>

#include "astroforces/core/types.hpp"

namespace astroforces::core {

/**
 * @brief Dense 3x3 matrix in row-major storage.
 */
struct Mat3 {
  std::array<double, 9> v{};
  [[nodiscard]] double& operator()(const int r, const int c) { return v[static_cast<std::size_t>(r * 3 + c)]; }
  [[nodiscard]] double operator()(const int r, const int c) const { return v[static_cast<std::size_t>(r * 3 + c)]; }
};

/**
 * @brief Create identity 3x3 matrix.
 */
inline Mat3 mat_identity() {
  Mat3 m{};
  m(0, 0) = 1.0;
  m(1, 1) = 1.0;
  m(2, 2) = 1.0;
  return m;
}

/**
 * @brief Matrix-matrix multiplication.
 */
inline Mat3 mat_mul(const Mat3& a, const Mat3& b) {
  Mat3 c{};
  for (int r = 0; r < 3; ++r) {
    for (int col = 0; col < 3; ++col) {
      c(r, col) = a(r, 0) * b(0, col) + a(r, 1) * b(1, col) + a(r, 2) * b(2, col);
    }
  }
  return c;
}

/**
 * @brief Matrix transpose.
 */
inline Mat3 mat_transpose(const Mat3& m) {
  Mat3 t{};
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      t(r, c) = m(c, r);
    }
  }
  return t;
}

/**
 * @brief Matrix-vector multiplication.
 */
inline Vec3 mat_vec(const Mat3& m, const Vec3& x) {
  return Vec3{
      m(0, 0) * x.x + m(0, 1) * x.y + m(0, 2) * x.z,
      m(1, 0) * x.x + m(1, 1) * x.y + m(1, 2) * x.z,
      m(2, 0) * x.x + m(2, 1) * x.y + m(2, 2) * x.z,
  };
}

/**
 * @brief Vector cross product.
 */
inline Vec3 vec_cross(const Vec3& a, const Vec3& b) {
  return Vec3{
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x,
  };
}

}  // namespace astroforces::core
