/**
 * @file solid_earth_tide.cpp
 * @brief Solid Earth tide contribution helpers for gravity SH models.
 * @author Watson
 */

#include "astroforces/forces/gravity/tides/solid_earth_tide.hpp"

#include <cmath>

namespace astroforces::forces::tides {
namespace {

struct LegendreCache {
  explicit LegendreCache(int nmax) : nmax(nmax) {
    anm = Eigen::MatrixXd::Zero(nmax + 1, nmax + 1);
    bnm = Eigen::MatrixXd::Zero(nmax + 1, nmax + 1);
    fnm = Eigen::MatrixXd::Zero(nmax + 1, nmax + 1);
    Pnm = Eigen::MatrixXd::Zero(nmax + 1, nmax + 1);
    dPnm = Eigen::MatrixXd::Zero(nmax + 1, nmax + 1);

    for (int n = 0; n <= nmax; ++n) {
      for (int m = 0; m <= n; ++m) {
        if (n - m > 0 && n + m > 0) {
          anm(n, m) = std::sqrt((static_cast<double>(2 * n - 1) * (2 * n + 1)) / ((n - m) * (n + m)));
        }
        if (n >= 2 && n - m > 0 && n + m > 0 && 2 * n - 3 > 0 && n + m - 1 >= 0 && n - m - 1 >= 0) {
          bnm(n, m) = std::sqrt((static_cast<double>(2 * n + 1) * (n + m - 1) * (n - m - 1)) /
                                ((n - m) * (n + m) * (2 * n - 3)));
        }
        if (n >= 1 && 2 * n - 1 > 0) {
          fnm(n, m) = std::sqrt((static_cast<double>(n * n - m * m) * (2 * n + 1)) / (2 * n - 1));
        }
      }
    }
  }

  bool calculate(double x) {
    if (std::abs(x) >= 1.0) {
      return false;
    }
    const double u = std::sqrt(1.0 - x * x);
    Pnm.setZero();
    dPnm.setZero();

    Pnm(0, 0) = 1.0;
    if (nmax == 0) {
      return true;
    }

    Pnm(1, 0) = std::sqrt(3.0) * x;
    Pnm(1, 1) = std::sqrt(3.0) * u;
    dPnm(1, 0) = (x * Pnm(1, 0) - std::sqrt(3.0) * Pnm(0, 0)) / u;
    dPnm(1, 1) = x * Pnm(1, 1) / u;

    for (int n = 2; n <= nmax; ++n) {
      for (int m = 0; m < n; ++m) {
        Pnm(n, m) = anm(n, m) * x * Pnm(n - 1, m) - bnm(n, m) * Pnm(n - 2, m);
        dPnm(n, m) = (n * x * Pnm(n, m) - fnm(n, m) * Pnm(n - 1, m)) / u;
      }
      Pnm(n, n) = u * std::sqrt((2.0 * n + 1.0) / (2.0 * n)) * Pnm(n - 1, n - 1);
      dPnm(n, n) = n * x * Pnm(n, n) / u;
    }

    return true;
  }

  int nmax{0};
  Eigen::MatrixXd anm{};
  Eigen::MatrixXd bnm{};
  Eigen::MatrixXd fnm{};
  Eigen::MatrixXd Pnm{};
  Eigen::MatrixXd dPnm{};
};

}  // namespace

void add_solid_earth_tide1_delta(const astroforces::core::Vec3& body_ecef_m,
                                 double mu_body_m3_s2,
                                 double mu_earth_m3_s2,
                                 double radius_m,
                                 int max_deg,
                                 Eigen::MatrixXd& dC,
                                 Eigen::MatrixXd& dS) {
  if (max_deg < 2) {
    return;
  }

  Eigen::Matrix<double, 5, 5> elasticLove;
  elasticLove << 0, 0, 0, 0, 0,
                 0, 0, 0, 0, 0,
                 0.29525, 0.29470, 0.29801, 0, 0,
                 0.09300, 0.09300, 0.09300, 0.09400, 0,
                 -0.00087, -0.00079, -0.00057, 0, 0;

  const double rb = astroforces::core::norm(body_ecef_m);
  if (!(rb > 0.0)) {
    return;
  }

  const double sin_lat = body_ecef_m.z / rb;
  const double lon = std::atan2(body_ecef_m.y, body_ecef_m.x);

  LegendreCache leg(3);
  if (!leg.calculate(sin_lat)) {
    return;
  }

  const int lim_23 = std::min(max_deg, 3);
  for (int n = 2; n <= lim_23; ++n) {
    for (int m = 0; m <= n; ++m) {
      const double cst = elasticLove(n, m) / (2.0 * n + 1.0);
      const double amp = cst * mu_body_m3_s2 / mu_earth_m3_s2 * std::pow(radius_m / rb, n + 1) * leg.Pnm(n, m);
      dC(n, m) += amp * std::cos(m * lon);
      if (m > 0) {
        dS(n, m) += amp * std::sin(m * lon);
      }
    }
  }

  if (max_deg >= 4) {
    for (int m = 0; m <= 2; ++m) {
      const double cst = elasticLove(4, m) / 5.0;
      const double amp = cst * mu_body_m3_s2 / mu_earth_m3_s2 * std::pow(radius_m / rb, 3) * leg.Pnm(2, m);
      dC(4, m) += amp * std::cos(m * lon);
      if (m > 0) {
        dS(4, m) += amp * std::sin(m * lon);
      }
    }
  }
}

}  // namespace astroforces::forces::tides
