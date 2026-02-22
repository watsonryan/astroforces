/**
 * @file constituent_tide.cpp
 * @brief Constituent ocean/atmospheric tide SH correction model.
 * @author Watosn
 */

#include "astroforces/forces/gravity/tides/constituent_tide.hpp"

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>

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

struct FundamentalArgs {
  double l_rad{};
  double lp_rad{};
  double f_rad{};
  double d_rad{};
  double omega_rad{};
};

FundamentalArgs fundamental_args(const double jd_utc) {
  const double T = (jd_utc - astroforces::core::constants::kJ2000Jd) / 36525.0;

  FundamentalArgs args{};
  args.l_rad = poly_arcsec_to_rad(134.96340251, 1717915923.2178, 31.8792, 0.051635, -0.00024470, T);
  args.lp_rad = poly_arcsec_to_rad(357.52910918, 129596581.0481, -0.5532, 0.000136, -0.00001149, T);
  args.f_rad = poly_arcsec_to_rad(93.27209062, 1739527262.8478, -12.7512, -0.001037, 0.00000417, T);
  args.d_rad = poly_arcsec_to_rad(297.85019547, 1602961601.2090, -6.3706, 0.006593, -0.00003169, T);
  args.omega_rad = poly_arcsec_to_rad(125.04455501, -6962890.5431, 7.4722, 0.007702, -0.00005939, T);
  return args;
}

std::array<double, 6> doodson_beta(const double jd_utc, const double gmst_rad) {
  const auto fa = fundamental_args(jd_utc);

  std::array<double, 6> beta{};
  beta[4] = -fa.omega_rad;
  beta[1] = fa.f_rad + fa.omega_rad;
  beta[0] = gmst_rad - beta[1];
  beta[2] = beta[1] - fa.d_rad;
  beta[3] = beta[1] - fa.l_rad;
  beta[5] = beta[1] - fa.d_rad - fa.lp_rad;

  for (double& b : beta) {
    b = wrap_0_2pi(b);
  }
  return beta;
}

std::array<double, 6> parse_doodson(const std::string& token) {
  std::array<double, 6> out{};
  const auto dot = token.find('.');
  if (dot == std::string::npos || dot < 3 || token.size() < dot + 4) {
    return out;
  }

  auto d = [&](const std::size_t idx) -> int {
    const char c = token[idx];
    return (c >= '0' && c <= '9') ? static_cast<int>(c - '0') : 0;
  };

  out[0] = static_cast<double>(d(dot - 3));
  out[1] = static_cast<double>(d(dot - 2) - 5);
  out[2] = static_cast<double>(d(dot - 1) - 5);
  out[3] = static_cast<double>(d(dot + 1) - 5);
  out[4] = static_cast<double>(d(dot + 2) - 5);
  out[5] = static_cast<double>(d(dot + 3) - 5);
  return out;
}

}  // namespace

std::unique_ptr<ConstituentTideModel> ConstituentTideModel::load_from_file(const std::filesystem::path& path,
                                                                            const int max_degree) {
  if (path.empty() || max_degree < 0) {
    return {};
  }

  std::ifstream in(path);
  if (!in.is_open()) {
    return {};
  }

  auto model = std::unique_ptr<ConstituentTideModel>(new ConstituentTideModel(max_degree));

  std::string line;
  for (int i = 0; i < 4 && std::getline(in, line); ++i) {
  }

  struct WaveWork {
    std::array<double, 6> doodson{};
    Eigen::MatrixXd CnmP;
    Eigen::MatrixXd CnmM;
    Eigen::MatrixXd SnmP;
    Eigen::MatrixXd SnmM;
  };

  std::unordered_map<std::string, WaveWork> wave_map;

  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }

    std::istringstream iss(line);
    std::string doodson_token;
    std::string wave_name;
    int n = 0;
    int m = 0;
    double cnmp = 0.0;
    double snmp = 0.0;
    double cnmm = 0.0;
    double snmm = 0.0;

    if (!(iss >> doodson_token >> wave_name >> n >> m >> cnmp >> snmp >> cnmm >> snmm)) {
      continue;
    }
    if (n < 0 || m < 0 || m > n || n > max_degree) {
      continue;
    }

    auto it = wave_map.find(wave_name);
    if (it == wave_map.end()) {
      WaveWork work{};
      work.doodson = parse_doodson(doodson_token);
      work.CnmP = Eigen::MatrixXd::Zero(max_degree + 1, max_degree + 1);
      work.CnmM = Eigen::MatrixXd::Zero(max_degree + 1, max_degree + 1);
      work.SnmP = Eigen::MatrixXd::Zero(max_degree + 1, max_degree + 1);
      work.SnmM = Eigen::MatrixXd::Zero(max_degree + 1, max_degree + 1);
      it = wave_map.emplace(wave_name, std::move(work)).first;
    }

    it->second.CnmP(n, m) = cnmp;
    it->second.SnmP(n, m) = snmp;
    it->second.CnmM(n, m) = cnmm;
    it->second.SnmM(n, m) = snmm;
  }

  model->waves_.reserve(wave_map.size());
  for (auto& [name, w] : wave_map) {
    Wave out{};
    out.name = name;
    out.doodson = w.doodson;
    out.C1 = w.CnmP + w.CnmM;
    out.C2 = w.SnmP + w.SnmM;
    out.S1 = w.SnmP - w.SnmM;
    out.S2 = w.CnmP - w.CnmM;
    model->waves_.push_back(std::move(out));
  }

  return model;
}

bool ConstituentTideModel::empty() const noexcept { return waves_.empty(); }

int ConstituentTideModel::max_degree() const noexcept { return max_degree_; }

void ConstituentTideModel::add_delta_coefficients(const double jd_utc,
                                                  const double gmst_rad,
                                                  Eigen::MatrixXd& dC,
                                                  Eigen::MatrixXd& dS) const {
  if (waves_.empty()) {
    return;
  }

  const auto beta = doodson_beta(jd_utc, gmst_rad);

  for (const auto& wave : waves_) {
    double theta = 0.0;
    for (int i = 0; i < 6; ++i) {
      theta += beta[i] * wave.doodson[i];
    }

    const double cos_theta = std::cos(theta);
    const double sin_theta = std::sin(theta);

    dC += (wave.C1 * cos_theta + wave.C2 * sin_theta) * 1e-11;
    dS += (wave.S1 * cos_theta - wave.S2 * sin_theta) * 1e-11;
  }

  if (dS.rows() > 0) {
    dS.col(0).setZero();
  }
  if (dS.cols() > 0) {
    dS.row(0).setZero();
  }
}

}  // namespace astroforces::forces::tides
