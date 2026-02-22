/**
 * @file solid_earth_tide_freqdep.cpp
 * @brief Frequency-dependent solid Earth tide (IERS2010 6.5 style) helpers.
 * @author Watosn
 */

#include "astroforces/forces/gravity/tides/solid_earth_tide_freqdep.hpp"

#include <cmath>

#include "astroforces/atmo/constants.hpp"

namespace astroforces::forces::tides {
namespace {

constexpr double kTwoPi = 6.283185307179586476925286766559;

double wrap_0_2pi(const double x) {
  const double y = std::fmod(x, kTwoPi);
  return (y < 0.0) ? (y + kTwoPi) : y;
}

double deg_to_rad(const double deg) { return deg * astroforces::core::constants::kDegToRad; }

double poly_arcsec_to_rad(const double base_deg,
                          const double c1,
                          const double c2,
                          const double c3,
                          const double c4,
                          const double T) {
  const double arcsec = c1 * T + c2 * T * T + c3 * T * T * T + c4 * T * T * T * T;
  return wrap_0_2pi(deg_to_rad(base_deg) + deg_to_rad(arcsec / 3600.0));
}

Eigen::Array<double, 6, 1> doodson_angles(const double jd_utc, const double gmst_rad) {
  const double T = (jd_utc - astroforces::core::constants::kJ2000Jd) / 36525.0;

  const double l = poly_arcsec_to_rad(134.96340251, 1717915923.2178, 31.8792, 0.051635, -0.00024470, T);
  const double lp = poly_arcsec_to_rad(357.52910918, 129596581.0481, -0.5532, 0.000136, -0.00001149, T);
  const double f = poly_arcsec_to_rad(93.27209062, 1739527262.8478, -12.7512, -0.001037, 0.00000417, T);
  const double d = poly_arcsec_to_rad(297.85019547, 1602961601.2090, -6.3706, 0.006593, -0.00003169, T);
  const double omega = poly_arcsec_to_rad(125.04455501, -6962890.5431, 7.4722, 0.007702, -0.00005939, T);

  Eigen::Array<double, 6, 1> out{};
  out(4) = -omega;
  out(1) = f + omega;
  out(0) = gmst_rad - out(1);
  out(2) = out(1) - d;
  out(3) = out(1) - l;
  out(5) = out(1) - d - lp;
  for (int i = 0; i < 6; ++i) {
    out(i) = wrap_0_2pi(out(i));
  }
  return out;
}

}  // namespace

void add_solid_earth_tide2_delta(const double jd_utc,
                                 const double gmst_rad,
                                 Eigen::MatrixXd& dC,
                                 Eigen::MatrixXd& dS) {
  if (dC.rows() <= 2 || dC.cols() <= 2 || dS.rows() <= 2 || dS.cols() <= 2) {
    return;
  }

  Eigen::Array<double, 48, 8> tab65a;
  Eigen::Array<double, 21, 8> tab65b;
  Eigen::Array<double, 2, 7> tab65c;

  tab65a <<
      1.0e0, -3.0e0,  0.0e0,  2.0e0,  0.0e0,  0.0e0,   -0.1e0,    0.0e0,
      1.0e0, -3.0e0,  2.0e0,  0.0e0,  0.0e0,  0.0e0,   -0.1e0,    0.0e0,
      1.0e0, -2.0e0,  0.0e0,  1.0e0, -1.0e0,  0.0e0,   -0.1e0,    0.0e0,
      1.0e0, -2.0e0,  0.0e0,  1.0e0,  0.0e0,  0.0e0,   -0.7e0,    0.1e0,
      1.0e0, -2.0e0,  2.0e0, -1.0e0,  0.0e0,  0.0e0,   -0.1e0,    0.0e0,
      1.0e0, -1.0e0,  0.0e0,  0.0e0, -1.0e0,  0.0e0,   -1.3e0,    0.1e0,
      1.0e0, -1.0e0,  0.0e0,  0.0e0,  0.0e0,  0.0e0,   -6.8e0,    0.6e0,
      1.0e0, -1.0e0,  2.0e0,  0.0e0,  0.0e0,  0.0e0,    0.1e0,    0.0e0,
      1.0e0,  0.0e0, -2.0e0,  1.0e0,  0.0e0,  0.0e0,    0.1e0,    0.0e0,
      1.0e0,  0.0e0,  0.0e0, -1.0e0, -1.0e0,  0.0e0,    0.1e0,    0.0e0,
      1.0e0,  0.0e0,  0.0e0, -1.0e0,  0.0e0,  0.0e0,    0.4e0,    0.0e0,
      1.0e0,  0.0e0,  0.0e0,  1.0e0,  0.0e0,  0.0e0,    1.3e0,   -0.1e0,
      1.0e0,  0.0e0,  0.0e0,  1.0e0,  1.0e0,  0.0e0,    0.3e0,    0.0e0,
      1.0e0,  0.0e0,  2.0e0, -1.0e0,  0.0e0,  0.0e0,    0.3e0,    0.0e0,
      1.0e0,  0.0e0,  2.0e0, -1.0e0,  1.0e0,  0.0e0,    0.1e0,    0.0e0,
      1.0e0,  1.0e0, -3.0e0,  0.0e0,  0.0e0,  1.0e0,   -1.9e0,    0.1e0,
      1.0e0,  1.0e0, -2.0e0,  0.0e0, -1.0e0,  0.0e0,    0.5e0,    0.0e0,
      1.0e0,  1.0e0, -2.0e0,  0.0e0,  0.0e0,  0.0e0,  -43.4e0,    2.9e0,
      1.0e0,  1.0e0, -1.0e0,  0.0e0,  0.0e0, -1.0e0,    0.6e0,    0.0e0,
      1.0e0,  1.0e0, -1.0e0,  0.0e0,  0.0e0,  1.0e0,    1.6e0,   -0.1e0,
      1.0e0,  1.0e0,  0.0e0, -2.0e0, -1.0e0,  0.0e0,    0.1e0,    0.0e0,
      1.0e0,  1.0e0,  0.0e0,  0.0e0, -2.0e0,  0.0e0,    0.1e0,    0.0e0,
      1.0e0,  1.0e0,  0.0e0,  0.0e0, -1.0e0,  0.0e0,   -8.8e0,    0.5e0,
      1.0e0,  1.0e0,  0.0e0,  0.0e0,  0.0e0,  0.0e0,  470.9e0,  -30.2e0,
      1.0e0,  1.0e0,  0.0e0,  0.0e0,  1.0e0,  0.0e0,   68.1e0,   -4.6e0,
      1.0e0,  1.0e0,  0.0e0,  0.0e0,  2.0e0,  0.0e0,   -1.6e0,    0.1e0,
      1.0e0,  1.0e0,  1.0e0, -1.0e0,  0.0e0,  0.0e0,    0.1e0,    0.0e0,
      1.0e0,  1.0e0,  1.0e0,  0.0e0, -1.0e0, -1.0e0,   -0.1e0,    0.0e0,
      1.0e0,  1.0e0,  1.0e0,  0.0e0,  0.0e0, -1.0e0,  -20.6e0,   -0.3e0,
      1.0e0,  1.0e0,  1.0e0,  0.0e0,  0.0e0,  1.0e0,    0.3e0,    0.0e0,
      1.0e0,  1.0e0,  1.0e0,  0.0e0,  1.0e0, -1.0e0,   -0.3e0,    0.0e0,
      1.0e0,  1.0e0,  2.0e0, -2.0e0,  0.0e0,  0.0e0,   -0.2e0,    0.0e0,
      1.0e0,  1.0e0,  2.0e0, -2.0e0,  1.0e0,  0.0e0,   -0.1e0,    0.0e0,
      1.0e0,  1.0e0,  2.0e0,  0.0e0,  0.0e0,  0.0e0,   -5.0e0,    0.3e0,
      1.0e0,  1.0e0,  2.0e0,  0.0e0,  1.0e0,  0.0e0,    0.2e0,    0.0e0,
      1.0e0,  1.0e0,  3.0e0,  0.0e0,  0.0e0, -1.0e0,   -0.2e0,    0.0e0,
      1.0e0,  2.0e0, -2.0e0,  1.0e0,  0.0e0,  0.0e0,   -0.5e0,    0.0e0,
      1.0e0,  2.0e0, -2.0e0,  1.0e0,  1.0e0,  0.0e0,   -0.1e0,    0.0e0,
      1.0e0,  2.0e0,  0.0e0, -1.0e0, -1.0e0,  0.0e0,    0.1e0,    0.0e0,
      1.0e0,  2.0e0,  0.0e0, -1.0e0,  0.0e0,  0.0e0,   -2.1e0,    0.1e0,
      1.0e0,  2.0e0,  0.0e0, -1.0e0,  1.0e0,  0.0e0,   -0.4e0,    0.0e0,
      1.0e0,  3.0e0, -2.0e0,  0.0e0,  0.0e0,  0.0e0,   -0.2e0,    0.0e0,
      1.0e0,  3.0e0,  0.0e0, -2.0e0,  0.0e0,  0.0e0,   -0.1e0,    0.0e0,
      1.0e0,  3.0e0,  0.0e0,  0.0e0,  0.0e0,  0.0e0,   -0.6e0,    0.0e0,
      1.0e0,  3.0e0,  0.0e0,  0.0e0,  1.0e0,  0.0e0,   -0.4e0,    0.0e0,
      1.0e0,  3.0e0,  0.0e0,  0.0e0,  2.0e0,  0.0e0,   -0.1e0,    0.0e0,
      1.0e0,  4.0e0,  0.0e0, -1.0e0,  0.0e0,  0.0e0,   -0.1e0,    0.0e0,
      1.0e0,  4.0e0,  0.0e0, -1.0e0,  1.0e0,  0.0e0,   -0.1e0,    0.0e0;

  tab65b <<
      0.0e0,  0.0e0,  0.0e0,  0.0e0,  1.0e0,  0.0e0,  16.6e0,  -6.7e0,
      0.0e0,  0.0e0,  0.0e0,  0.0e0,  2.0e0,  0.0e0,  -0.1e0,   0.1e0,
      0.0e0,  0.0e0,  1.0e0,  0.0e0,  0.0e0, -1.0e0,  -1.2e0,   0.8e0,
      0.0e0,  0.0e0,  2.0e0,  0.0e0,  0.0e0,  0.0e0,  -5.5e0,   4.3e0,
      0.0e0,  0.0e0,  2.0e0,  0.0e0,  1.0e0,  0.0e0,   0.1e0,  -0.1e0,
      0.0e0,  0.0e0,  3.0e0,  0.0e0,  0.0e0, -1.0e0,  -0.3e0,   0.2e0,
      0.0e0,  1.0e0, -2.0e0,  1.0e0,  0.0e0,  0.0e0,  -0.3e0,   0.7e0,
      0.0e0,  1.0e0,  0.0e0, -1.0e0, -1.0e0,  0.0e0,   0.1e0,  -0.2e0,
      0.0e0,  1.0e0,  0.0e0, -1.0e0,  0.0e0,  0.0e0,  -1.2e0,   3.7e0,
      0.0e0,  1.0e0,  0.0e0, -1.0e0,  1.0e0,  0.0e0,   0.1e0,  -0.2e0,
      0.0e0,  1.0e0,  0.0e0,  1.0e0,  0.0e0,  0.0e0,   0.1e0,  -0.2e0,
      0.0e0,  2.0e0, -2.0e0,  0.0e0,  0.0e0,  0.0e0,   0.0e0,   0.6e0,
      0.0e0,  2.0e0,  0.0e0, -2.0e0,  0.0e0,  0.0e0,   0.0e0,   0.3e0,
      0.0e0,  2.0e0,  0.0e0,  0.0e0,  0.0e0,  0.0e0,   0.6e0,   6.3e0,
      0.0e0,  2.0e0,  0.0e0,  0.0e0,  1.0e0,  0.0e0,   0.2e0,   2.6e0,
      0.0e0,  2.0e0,  0.0e0,  0.0e0,  2.0e0,  0.0e0,   0.0e0,   0.2e0,
      0.0e0,  3.0e0, -2.0e0,  1.0e0,  0.0e0,  0.0e0,   0.1e0,   0.2e0,
      0.0e0,  3.0e0,  0.0e0, -1.0e0,  0.0e0,  0.0e0,   0.4e0,   1.1e0,
      0.0e0,  3.0e0,  0.0e0, -1.0e0,  1.0e0,  0.0e0,   0.2e0,   0.5e0,
      0.0e0,  4.0e0, -2.0e0,  0.0e0,  0.0e0,  0.0e0,   0.1e0,   0.2e0,
      0.0e0,  4.0e0,  0.0e0, -2.0e0,  0.0e0,  0.0e0,   0.1e0,   0.1e0;

  tab65c <<
      2.0e0, -1.0e0,  0.0e0,  1.0e0,  0.0e0,  0.0e0,   -0.3e0,
      2.0e0,  0.0e0,  0.0e0,  0.0e0,  0.0e0,  0.0e0,   -1.2e0;

  const auto dood = doodson_angles(jd_utc, gmst_rad);

  for (int i = 0; i < tab65b.rows(); ++i) {
    double theta = 0.0;
    for (int k = 0; k < 6; ++k) {
      theta += tab65b(i, k) * dood(k);
    }
    dC(2, 0) += tab65b(i, 6) * 1e-12 * std::cos(theta) - tab65b(i, 7) * 1e-12 * std::sin(theta);
  }

  for (int i = 0; i < tab65a.rows(); ++i) {
    double theta = 0.0;
    for (int k = 0; k < 6; ++k) {
      theta += tab65a(i, k) * dood(k);
    }
    dC(2, 1) += tab65a(i, 6) * 1e-12 * std::sin(theta) + tab65a(i, 7) * 1e-12 * std::cos(theta);
    dS(2, 1) += tab65a(i, 6) * 1e-12 * std::cos(theta) - tab65a(i, 7) * 1e-12 * std::sin(theta);
  }

  for (int i = 0; i < tab65c.rows(); ++i) {
    double theta = 0.0;
    for (int k = 0; k < 6; ++k) {
      theta += tab65c(i, k) * dood(k);
    }
    dC(2, 2) += tab65c(i, 6) * 1e-12 * std::cos(theta);
    dS(2, 2) += -tab65c(i, 6) * 1e-12 * std::sin(theta);
  }
}

}  // namespace astroforces::forces::tides
