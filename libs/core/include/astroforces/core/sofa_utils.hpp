/**
 * @file sofa_utils.hpp
 * @brief SOFA-style Earth orientation building blocks (modern C++ implementation).
 * @author Watosn
 */
#pragma once

#include <cmath>

#include "astroforces/core/constants.hpp"
#include "astroforces/core/math_utils.hpp"

namespace astroforces::core::sofa {

inline Mat3 rot_x(const double a) {
  Mat3 m = mat_identity();
  const double c = std::cos(a);
  const double s = std::sin(a);
  m(1, 1) = c;
  m(1, 2) = s;
  m(2, 1) = -s;
  m(2, 2) = c;
  return m;
}

inline Mat3 rot_y(const double a) {
  Mat3 m = mat_identity();
  const double c = std::cos(a);
  const double s = std::sin(a);
  m(0, 0) = c;
  m(0, 2) = -s;
  m(2, 0) = s;
  m(2, 2) = c;
  return m;
}

inline Mat3 rot_z(const double a) {
  Mat3 m = mat_identity();
  const double c = std::cos(a);
  const double s = std::sin(a);
  m(0, 0) = c;
  m(0, 1) = s;
  m(1, 0) = -s;
  m(1, 1) = c;
  return m;
}

inline Mat3 drot_x_da(const double a) {
  Mat3 m{};
  const double c = std::cos(a);
  const double s = std::sin(a);
  m(1, 1) = -s;
  m(1, 2) = c;
  m(2, 1) = -c;
  m(2, 2) = -s;
  return m;
}

inline Mat3 drot_y_da(const double a) {
  Mat3 m{};
  const double c = std::cos(a);
  const double s = std::sin(a);
  m(0, 0) = -s;
  m(0, 2) = -c;
  m(2, 0) = c;
  m(2, 2) = -s;
  return m;
}

inline Mat3 drot_z_da(const double a) {
  Mat3 m{};
  const double c = std::cos(a);
  const double s = std::sin(a);
  m(0, 0) = -s;
  m(0, 1) = c;
  m(1, 0) = -c;
  m(1, 1) = -s;
  return m;
}

// SOFA iauEra00-compatible expression.
inline double era00(const double jd_ut1) {
  double era = constants::kTwoPi * (0.7790572732640 + 1.00273781191135448 * (jd_ut1 - constants::kJ2000Jd));
  era = std::fmod(era, constants::kTwoPi);
  if (era < 0.0) {
    era += constants::kTwoPi;
  }
  return era;
}

// SOFA iauSp00-compatible linear approximation.
inline double sp00(const double jd_tt) {
  const double t = (jd_tt - constants::kJ2000Jd) / 36525.0;
  return -47e-6 * constants::kArcsecToRad * t;
}

// SOFA iauC2ixys-style construction.
inline Mat3 c2ixys(const double x_rad, const double y_rad, const double s_rad) {
  const double r2 = x_rad * x_rad + y_rad * y_rad;
  const double e = (r2 > 0.0) ? std::atan2(y_rad, x_rad) : 0.0;
  const double d = std::atan(std::sqrt(r2 / std::max(1e-30, 1.0 - r2)));
  return mat_mul(rot_z(-(e + s_rad)), mat_mul(rot_y(d), rot_z(e)));
}

// SOFA iauPom00-style polar-motion matrix.
inline Mat3 pom00(const double xp_rad, const double yp_rad, const double sp_rad) {
  return mat_mul(rot_z(sp_rad), mat_mul(rot_y(-xp_rad), rot_x(-yp_rad)));
}

// SOFA iauC2tcio-style final matrix: RPOM * R3(ERA) * RC2I.
inline Mat3 c2tcio(const Mat3& rc2i, const double era_rad, const Mat3& rpom) {
  return mat_mul(rpom, mat_mul(rot_z(era_rad), rc2i));
}

// SOFA-style compact periodic model for TDB-TT (seconds), suitable for force
// model epoch conversion where microsecond-level effects matter.
inline double dtdb_seconds_approx(const double jd_tt) {
  const double t = (jd_tt - constants::kJ2000Jd) / 36525.0;
  const double g = (357.5277233 + 35999.05034 * t) * constants::kDegToRad;   // Mean anomaly of the Sun.
  const double l = (280.4665 + 36000.7698 * t) * constants::kDegToRad;        // Mean longitude of the Sun.
  return 0.001657 * std::sin(g) + 0.000022 * std::sin(2.0 * g) + 0.000014 * std::sin(3.0 * g) + 0.000005 * std::sin(4.0 * g)
         + 0.000005 * std::sin(l) + 0.000002 * std::sin(2.0 * l);
}

}  // namespace astroforces::core::sofa
