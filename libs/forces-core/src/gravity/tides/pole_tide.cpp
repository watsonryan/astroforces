/**
 * @file pole_tide.cpp
 * @brief Pole tide correction helpers for gravity SH models.
 * @author Watosn
 */

#include "astroforces/forces/gravity/tides/pole_tide.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "astroforces/atmo/constants.hpp"

namespace astroforces::forces::tides {
namespace {

void secular_pole_mas(const double mjd_tt, double& xpv_mas, double& ypv_mas) {
  const double t_years = (mjd_tt - 51544.5) / 365.25;
  xpv_mas = 55.0 + 1.677 * t_years;
  ypv_mas = 320.5 + 3.460 * t_years;
}

std::vector<std::string> split_ws(const std::string& s) {
  std::istringstream iss(s);
  std::vector<std::string> out;
  for (std::string tok; iss >> tok;) {
    out.push_back(tok);
  }
  return out;
}

}  // namespace

std::unique_ptr<OceanPoleTideModel> OceanPoleTideModel::load_from_file(const std::filesystem::path& path, int max_degree) {
  std::ifstream in(path);
  if (!in) {
    return {};
  }
  auto out = std::unique_ptr<OceanPoleTideModel>(new OceanPoleTideModel(std::max(0, max_degree)));
  const int nmax = out->max_degree_;
  out->cnmp_ = Eigen::MatrixXd::Zero(nmax + 1, nmax + 1);
  out->cnmm_ = Eigen::MatrixXd::Zero(nmax + 1, nmax + 1);
  out->snmp_ = Eigen::MatrixXd::Zero(nmax + 1, nmax + 1);
  out->snmm_ = Eigen::MatrixXd::Zero(nmax + 1, nmax + 1);

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    const auto toks = split_ws(line);
    if (toks.size() < 6) {
      continue;
    }
    const int n = std::stoi(toks[0]);
    const int m = std::stoi(toks[1]);
    if (n < 0 || m < 0 || m > n || n > nmax) {
      continue;
    }
    out->cnmp_(n, m) = std::stod(toks[2]);
    out->cnmm_(n, m) = std::stod(toks[3]);
    out->snmp_(n, m) = std::stod(toks[4]);
    out->snmm_(n, m) = std::stod(toks[5]);
  }
  return out;
}

bool OceanPoleTideModel::empty() const noexcept {
  return cnmp_.size() == 0 || cnmm_.size() == 0 || snmp_.size() == 0 || snmm_.size() == 0;
}

void OceanPoleTideModel::add_delta_coefficients(const double m1_arcsec,
                                                const double m2_arcsec,
                                                Eigen::MatrixXd& dC,
                                                Eigen::MatrixXd& dS) const {
  if (empty() || dC.rows() == 0 || dC.cols() == 0 || dS.rows() == 0 || dS.cols() == 0) {
    return;
  }
  constexpr double gamma2_r = 0.6870;
  constexpr double gamma2_i = 0.0036;
  const double coeff1 = (m1_arcsec * gamma2_r + m2_arcsec * gamma2_i) * astroforces::core::constants::kArcsecToRad;
  const double coeff2 = (m2_arcsec * gamma2_r - m1_arcsec * gamma2_i) * astroforces::core::constants::kArcsecToRad;
  const int nmax = std::min<int>(max_degree_, std::min<int>(dC.rows() - 1, dC.cols() - 1));
  dC.topLeftCorner(nmax + 1, nmax + 1) += coeff1 * cnmp_.topLeftCorner(nmax + 1, nmax + 1)
                                         + coeff2 * cnmm_.topLeftCorner(nmax + 1, nmax + 1);
  dS.topLeftCorner(nmax + 1, nmax + 1) += coeff1 * snmp_.topLeftCorner(nmax + 1, nmax + 1)
                                         + coeff2 * snmm_.topLeftCorner(nmax + 1, nmax + 1);
}

void add_pole_solid_tide_delta(const double mjd_tt,
                               const double xp_rad,
                               const double yp_rad,
                               Eigen::MatrixXd& dC,
                               Eigen::MatrixXd& dS) {
  if (dC.rows() <= 2 || dC.cols() <= 1 || dS.rows() <= 2 || dS.cols() <= 1) {
    return;
  }

  double xpv_mas = 0.0;
  double ypv_mas = 0.0;
  secular_pole_mas(mjd_tt, xpv_mas, ypv_mas);

  const double m1 = +(xp_rad / astroforces::core::constants::kArcsecToRad - xpv_mas / 1000.0);
  const double m2 = -(yp_rad / astroforces::core::constants::kArcsecToRad - ypv_mas / 1000.0);

  dC(2, 1) += -1.333e-9 * (m1 + 0.0115 * m2);
  dS(2, 1) += -1.333e-9 * (m2 - 0.0115 * m1);
}

void add_pole_ocean_tide_delta(const double mjd_tt,
                               const double xp_rad,
                               const double yp_rad,
                               Eigen::MatrixXd& dC,
                               Eigen::MatrixXd& dS) {
  if (dC.rows() <= 2 || dC.cols() <= 1 || dS.rows() <= 2 || dS.cols() <= 1) {
    return;
  }

  double xpv_mas = 0.0;
  double ypv_mas = 0.0;
  secular_pole_mas(mjd_tt, xpv_mas, ypv_mas);

  const double m1 = +(xp_rad / astroforces::core::constants::kArcsecToRad - xpv_mas / 1000.0);
  const double m2 = -(yp_rad / astroforces::core::constants::kArcsecToRad - ypv_mas / 1000.0);

  dC(2, 1) += -2.1778e-10 * (m1 - 0.01724 * m2);
  dS(2, 1) += -1.7232e-10 * (m2 - 0.03365 * m1);
}

}  // namespace astroforces::forces::tides
