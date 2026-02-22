/**
 * @file gravity_sph_model.cpp
 * @brief Full spherical harmonic gravity model implementation.
 * @author Watosn
 */

#include "astroforces/forces/gravity/gravity_sph_model.hpp"

#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include "astroforces/core/transforms.hpp"
#include "astroforces/forces/gravity/tides/aod1b_tide.hpp"
#include "astroforces/forces/gravity/tides/constituent_tide.hpp"
#include "astroforces/forces/gravity/tides/pole_tide.hpp"
#include "astroforces/forces/gravity/tides/solid_earth_tide.hpp"
#include "astroforces/forces/gravity/tides/solid_earth_tide_freqdep.hpp"
#include "jpl_eph/jpl_eph.hpp"

namespace astroforces::forces {

struct GravitySphAccelerationModel::GravityFieldData {
  Eigen::MatrixXd C{};
  Eigen::MatrixXd S{};
  struct TimeTerm {
    double c{0.0};
    double s{0.0};
    double t0_jd{0.0};
    double t1_jd{0.0};
    double period_years{0.0};
  };
  struct TemporalSeries {
    std::vector<TimeTerm> gfct{};
    std::vector<TimeTerm> trnd{};
    std::vector<TimeTerm> acos{};
    std::vector<TimeTerm> asin{};
  };
  std::unordered_map<std::int64_t, TemporalSeries> temporal{};
  bool has_time_variable{false};
  int max_degree{0};
  double mu_earth_m3_s2{astroforces::core::constants::kEarthMuM3S2};
  double radius_m{astroforces::core::constants::kEarthEquatorialRadiusM};
  bool is_tide_free{true};
};

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

astroforces::core::Status map_jpl_error(const jpl::eph::Status& s) {
  switch (s.code) {
    case jpl::eph::ErrorCode::kInvalidArgument:
      return astroforces::core::Status::InvalidInput;
    case jpl::eph::ErrorCode::kIo:
    case jpl::eph::ErrorCode::kCorruptFile:
    case jpl::eph::ErrorCode::kOutOfRange:
    case jpl::eph::ErrorCode::kUnsupported:
      return astroforces::core::Status::DataUnavailable;
    case jpl::eph::ErrorCode::kOk:
    default:
      return astroforces::core::Status::NumericalError;
  }
}

astroforces::core::Vec3 to_vec3(const std::array<double, 6>& pv) { return astroforces::core::Vec3{pv[0], pv[1], pv[2]}; }

std::vector<std::string> split_ws(const std::string& s) {
  std::istringstream iss(s);
  std::vector<std::string> out;
  for (std::string tok; iss >> tok;) {
    out.push_back(tok);
  }
  return out;
}

double icgem_date_token_to_jd(const std::string& token) {
  const std::size_t dot = token.find('.');
  const std::string ymd = (dot == std::string::npos) ? token : token.substr(0, dot);
  const std::string frac = (dot == std::string::npos) ? std::string{} : token.substr(dot + 1);
  if (ymd.size() < 8) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const int year = std::stoi(ymd.substr(0, 4));
  const int month = std::stoi(ymd.substr(4, 2));
  const int day = std::stoi(ymd.substr(6, 2));

  double day_frac = 0.0;
  if (!frac.empty()) {
    day_frac = std::stod("0." + frac);
  }
  day_frac = std::clamp(day_frac, 0.0, 0.999999999);

  int y = year;
  int m = month;
  int d = day;
  if (m <= 2) {
    y -= 1;
    m += 12;
  }
  const int A = y / 100;
  const int B = 2 - A + (A / 4);
  const double jd = std::floor(365.25 * (y + 4716)) + std::floor(30.6001 * (m + 1)) + d + day_frac + B - 1524.5;
  return jd;
}

std::int64_t nm_key(int n, int m) {
  return (static_cast<std::int64_t>(n) << 32) | static_cast<std::uint32_t>(m);
}

bool in_window(double jd, double t0_jd, double t1_jd) {
  return jd >= t0_jd && (jd < t1_jd || std::abs(jd - t1_jd) < 1e-12);
}

const GravitySphAccelerationModel::GravityFieldData::TimeTerm* select_active_term(
    const std::vector<GravitySphAccelerationModel::GravityFieldData::TimeTerm>& terms,
    double jd) {
  const GravitySphAccelerationModel::GravityFieldData::TimeTerm* best = nullptr;
  for (const auto& t : terms) {
    if (!in_window(jd, t.t0_jd, t.t1_jd)) {
      continue;
    }
    if (best == nullptr || t.t0_jd > best->t0_jd) {
      best = &t;
    }
  }
  return best;
}

void apply_time_variable_coefficients(
    const GravitySphAccelerationModel::GravityFieldData& field,
    double jd_utc,
    int nmax,
    Eigen::MatrixXd& Ceff,
    Eigen::MatrixXd& Seff) {
  constexpr double kDaysPerYear = 365.25;
  constexpr double kTwoPi = 6.283185307179586476925286766559;

  for (const auto& [key, ts] : field.temporal) {
    const int n = static_cast<int>(key >> 32);
    const int m = static_cast<int>(static_cast<std::uint32_t>(key & 0xFFFFFFFFULL));
    if (n < 0 || m < 0 || m > n || n > nmax || n >= Ceff.rows() || m >= Ceff.cols()) {
      continue;
    }

    double c = Ceff(n, m);
    double s = Seff(n, m);

    if (const auto* base = select_active_term(ts.gfct, jd_utc); base != nullptr) {
      c = base->c;
      s = base->s;
    }

    if (const auto* trend = select_active_term(ts.trnd, jd_utc); trend != nullptr) {
      const double dt_years = (jd_utc - trend->t0_jd) / kDaysPerYear;
      c += trend->c * dt_years;
      s += trend->s * dt_years;
    }

    for (const auto& ac : ts.acos) {
      if (!in_window(jd_utc, ac.t0_jd, ac.t1_jd) || !(ac.period_years > 0.0)) {
        continue;
      }
      const double dt_years = (jd_utc - ac.t0_jd) / kDaysPerYear;
      const double ph = kTwoPi * dt_years / ac.period_years;
      const double f = std::cos(ph);
      c += ac.c * f;
      s += ac.s * f;
    }

    for (const auto& as : ts.asin) {
      if (!in_window(jd_utc, as.t0_jd, as.t1_jd) || !(as.period_years > 0.0)) {
        continue;
      }
      const double dt_years = (jd_utc - as.t0_jd) / kDaysPerYear;
      const double ph = kTwoPi * dt_years / as.period_years;
      const double f = std::sin(ph);
      c += as.c * f;
      s += as.s * f;
    }

    Ceff(n, m) = c;
    if (m > 0) {
      Seff(n, m) = s;
    }
  }
}

std::shared_ptr<GravitySphAccelerationModel::GravityFieldData> load_gravity_field(
    const std::filesystem::path& file,
    int max_degree,
    bool convert_to_tide_free) {
  if (file.empty() || max_degree < 0) {
    return {};
  }

  std::ifstream in(file);
  if (!in) {
    return {};
  }

  auto data = std::make_shared<GravitySphAccelerationModel::GravityFieldData>();
  data->C = Eigen::MatrixXd::Zero(max_degree + 1, max_degree + 1);
  data->S = Eigen::MatrixXd::Zero(max_degree + 1, max_degree + 1);
  data->C(0, 0) = 1.0;

  bool in_header = true;
  int degree_seen = 0;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    const auto toks = split_ws(line);
    if (toks.empty()) {
      continue;
    }

    if (in_header) {
      if (toks[0] == "end_of_head") {
        in_header = false;
        continue;
      }
      if (toks[0] == "earth_gravity_constant" && toks.size() > 1) {
        data->mu_earth_m3_s2 = std::stod(toks[1]);
      } else if (toks[0] == "radius" && toks.size() > 1) {
        data->radius_m = std::stod(toks[1]);
      } else if (toks[0] == "tide_system" && toks.size() > 1) {
        data->is_tide_free = (toks[1] == "tide_free");
      }
      continue;
    }

    if (toks[0] != "gfc" && toks[0] != "gfct" && toks[0] != "trnd" && toks[0] != "acos" && toks[0] != "asin") {
      continue;
    }
    if (toks.size() < 5) {
      continue;
    }

    const int n = std::stoi(toks[1]);
    const int m = std::stoi(toks[2]);
    if (n < 0 || m < 0 || m > n || n > max_degree) {
      continue;
    }
    const double c = std::stod(toks[3]);
    const double s = std::stod(toks[4]);

    if (toks[0] == "gfc") {
      data->C(n, m) = c;
      if (m > 0) {
        data->S(n, m) = s;
      }
      degree_seen = std::max(degree_seen, n);
      continue;
    }

    if (toks[0] == "gfct" || toks[0] == "trnd") {
      if (toks.size() < 9) {
        continue;
      }
      const double t0_jd = icgem_date_token_to_jd(toks[7]);
      const double t1_jd = icgem_date_token_to_jd(toks[8]);
      if (!std::isfinite(t0_jd) || !std::isfinite(t1_jd) || t1_jd <= t0_jd) {
        continue;
      }
      GravitySphAccelerationModel::GravityFieldData::TimeTerm term{};
      term.c = c;
      term.s = s;
      term.t0_jd = t0_jd;
      term.t1_jd = t1_jd;
      auto& ts = data->temporal[nm_key(n, m)];
      if (toks[0] == "gfct") {
        ts.gfct.push_back(term);
      } else {
        ts.trnd.push_back(term);
      }
      data->has_time_variable = true;
      degree_seen = std::max(degree_seen, n);
      continue;
    }

    if (toks[0] == "acos" || toks[0] == "asin") {
      if (toks.size() < 10) {
        continue;
      }
      const double t0_jd = icgem_date_token_to_jd(toks[7]);
      const double t1_jd = icgem_date_token_to_jd(toks[8]);
      const double period_years = std::stod(toks[9]);
      if (!std::isfinite(t0_jd) || !std::isfinite(t1_jd) || t1_jd <= t0_jd || !(period_years > 0.0)) {
        continue;
      }
      GravitySphAccelerationModel::GravityFieldData::TimeTerm term{};
      term.c = c;
      term.s = s;
      term.t0_jd = t0_jd;
      term.t1_jd = t1_jd;
      term.period_years = period_years;
      auto& ts = data->temporal[nm_key(n, m)];
      if (toks[0] == "acos") {
        ts.acos.push_back(term);
      } else {
        ts.asin.push_back(term);
      }
      data->has_time_variable = true;
      degree_seen = std::max(degree_seen, n);
      continue;
    }
  }

  data->max_degree = degree_seen;
  if (!data->is_tide_free && convert_to_tide_free && data->max_degree >= 2) {
    data->C(2, 0) -= -4.1736e-9;
    data->is_tide_free = true;
  }

  return data;
}

astroforces::core::Vec3 central_accel(const astroforces::core::Vec3& r_m, double mu_earth_m3_s2) {
  const double r2 = astroforces::core::dot(r_m, r_m);
  if (!(r2 > 0.0)) {
    return {};
  }
  const double r = std::sqrt(r2);
  const double r3 = r2 * r;
  return (-mu_earth_m3_s2 / r3) * r_m;
}

astroforces::core::Vec3 accel_sph_noncentral(const astroforces::core::Vec3& r_ecef_m,
                                             const Eigen::MatrixXd& C,
                                             const Eigen::MatrixXd& S,
                                             int max_deg,
                                             double mu_earth_m3_s2,
                                             double radius_m) {
  const double R = astroforces::core::norm(r_ecef_m);
  if (!(R > 0.0)) {
    return {};
  }
  const double sin_lat = r_ecef_m.z / R;
  const double Rxy = std::sqrt(r_ecef_m.x * r_ecef_m.x + r_ecef_m.y * r_ecef_m.y);
  if (!(Rxy > 0.0)) {
    return {};
  }

  const double cos_lon = r_ecef_m.x / Rxy;
  const double sin_lon = r_ecef_m.y / Rxy;

  Eigen::VectorXd cosphi = Eigen::VectorXd::Zero(max_deg + 1);
  Eigen::VectorXd sinphi = Eigen::VectorXd::Zero(max_deg + 1);
  cosphi(0) = 1.0;
  if (max_deg >= 1) {
    cosphi(1) = cos_lon;
    sinphi(1) = sin_lon;
    for (int m = 2; m <= max_deg; ++m) {
      cosphi(m) = cosphi(m - 1) * cos_lon - sinphi(m - 1) * sin_lon;
      sinphi(m) = sinphi(m - 1) * cos_lon + cosphi(m - 1) * sin_lon;
    }
  }

  LegendreCache leg(max_deg);
  if (!leg.calculate(sin_lat)) {
    return {};
  }

  double dVr = 0.0;
  double dVtheta = 0.0;
  double dVlambda = 0.0;

  double ratio_pow = radius_m / R;
  for (int n = 2; n <= max_deg; ++n) {
    ratio_pow *= radius_m / R;
    double dVr_n = 0.0;
    double dVtheta_n = 0.0;
    double dVlambda_n = 0.0;

    for (int m = 0; m <= n; ++m) {
      const double common = C(n, m) * cosphi(m) + S(n, m) * sinphi(m);
      dVr_n += leg.Pnm(n, m) * common;
      dVtheta_n += leg.dPnm(n, m) * common;
      dVlambda_n += static_cast<double>(m) * leg.Pnm(n, m) * (S(n, m) * cosphi(m) - C(n, m) * sinphi(m));
    }

    dVr += -(n + 1.0) * ratio_pow * dVr_n;
    dVtheta += ratio_pow * dVtheta_n;
    dVlambda += ratio_pow * dVlambda_n;
  }

  const double constant = mu_earth_m3_s2 / R;
  dVtheta *= constant;
  dVlambda *= constant;
  dVr *= mu_earth_m3_s2 / (R * R);

  astroforces::core::Vec3 acc{};
  acc.x = dVr * r_ecef_m.x / R + dVtheta * r_ecef_m.z * r_ecef_m.x / (R * R * Rxy) - dVlambda * r_ecef_m.y / (Rxy * Rxy);
  acc.y = dVr * r_ecef_m.y / R + dVtheta * r_ecef_m.z * r_ecef_m.y / (R * R * Rxy) + dVlambda * r_ecef_m.x / (Rxy * Rxy);
  acc.z = dVr * r_ecef_m.z / R - dVtheta * Rxy / (R * R);
  return acc;
}

double gmst_rad_from_jd_utc(double jd_utc) {
  const double T = (jd_utc - 2451545.0) / 36525.0;
  const double gmst_deg = 280.46061837 + 360.98564736629 * (jd_utc - 2451545.0) + 0.000387933 * T * T - (T * T * T) / 38710000.0;
  constexpr double kPi = 3.1415926535897932384626433832795;
  const double gmst_wrapped = std::fmod(gmst_deg, 360.0);
  return (gmst_wrapped < 0.0 ? gmst_wrapped + 360.0 : gmst_wrapped) * kPi / 180.0;
}

astroforces::core::Vec3 rot_z(double theta, const astroforces::core::Vec3& v) {
  const double c = std::cos(theta);
  const double s = std::sin(theta);
  return astroforces::core::Vec3{c * v.x + s * v.y, -s * v.x + c * v.y, v.z};
}

bool finite_vec(const astroforces::core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

}  // namespace

std::unique_ptr<GravitySphAccelerationModel> GravitySphAccelerationModel::Create(const Config& config) {
  auto out = std::unique_ptr<GravitySphAccelerationModel>(new GravitySphAccelerationModel(config));

  if (config.use_sph && !config.gravity_model_file.empty()) {
    out->field_ = load_gravity_field(config.gravity_model_file, std::max(0, config.max_degree), config.convert_to_tide_free);
  }

  if (config.use_solid_earth_tides && (config.use_sun_tide || config.use_moon_tide)) {
    auto opened = jpl::eph::Ephemeris::Open(config.ephemeris_file.string());
    if (opened.has_value()) {
      out->ephemeris_ = opened.value();
      out->workspace_ = std::make_shared<jpl::eph::Workspace>();
    }
  }

  if (config.use_pole_tide_solid || config.use_pole_tide_ocean) {
    if (!config.eop_finals_file.empty()) {
      auto eop_series = std::make_shared<astroforces::core::eop::Series>(
          astroforces::core::eop::Series::load_iers_finals(config.eop_finals_file));
      if (!eop_series->empty()) {
        out->eop_ = std::move(eop_series);
      }
    }
    if (config.use_pole_tide_ocean && !config.ocean_pole_tide_file.empty()) {
      auto model = tides::OceanPoleTideModel::load_from_file(config.ocean_pole_tide_file, std::max(0, config.max_degree));
      if (model && !model->empty()) {
        out->ocean_pole_tide_ = std::shared_ptr<tides::OceanPoleTideModel>(std::move(model));
      }
    }
  }

  if (config.use_ocean_tide && !config.ocean_tide_file.empty()) {
    auto model = tides::ConstituentTideModel::load_from_file(config.ocean_tide_file, std::max(0, config.max_degree));
    if (model && !model->empty()) {
      out->ocean_tide_ = std::shared_ptr<tides::ConstituentTideModel>(std::move(model));
    }
  }

  if (config.use_atmos_tide && !config.atmos_tide_file.empty()) {
    auto model = tides::ConstituentTideModel::load_from_file(config.atmos_tide_file, std::max(0, config.max_degree));
    if (model && !model->empty()) {
      out->atmos_tide_ = std::shared_ptr<tides::ConstituentTideModel>(std::move(model));
    }
  }

  if (config.use_aod && !config.aod_file.empty()) {
    auto model = tides::Aod1bTideModel::load_from_file(config.aod_file, std::max(0, config.max_degree));
    if (model && !model->empty()) {
      out->aod_tide_ = std::shared_ptr<tides::Aod1bTideModel>(std::move(model));
    }
  }

  return out;
}

GravitySphResult GravitySphAccelerationModel::evaluate(const astroforces::core::StateVector& state) const {
  if (state.frame != astroforces::core::Frame::ECI && state.frame != astroforces::core::Frame::ECEF) {
    return GravitySphResult{.status = astroforces::core::Status::InvalidInput};
  }
  if (config_.max_degree < 0) {
    return GravitySphResult{.status = astroforces::core::Status::InvalidInput};
  }

  const double jd_utc = astroforces::core::utc_seconds_to_julian_date_utc(state.epoch.utc_seconds);
  const double gmst = gmst_rad_from_jd_utc(jd_utc);

  astroforces::core::Vec3 r_ecef = state.position_m;
  if (state.frame == astroforces::core::Frame::ECI) {
    if (!config_.use_simple_eci_to_ecef) {
      return GravitySphResult{.status = astroforces::core::Status::InvalidInput};
    }
    r_ecef = rot_z(gmst, state.position_m);
  }

  GravitySphResult out{};

  const double mu_earth = (config_.mu_earth_m3_s2 > 0.0)
                              ? config_.mu_earth_m3_s2
                              : (field_ ? field_->mu_earth_m3_s2 : astroforces::core::constants::kEarthMuM3S2);
  const double radius_m = (config_.earth_equatorial_radius_m > 0.0)
                              ? config_.earth_equatorial_radius_m
                              : (field_ ? field_->radius_m : astroforces::core::constants::kEarthEquatorialRadiusM);

  if (!(mu_earth > 0.0) || !(radius_m > 0.0)) {
    return GravitySphResult{.status = astroforces::core::Status::InvalidInput};
  }

  if (config_.use_central) {
    out.central_mps2 = central_accel(r_ecef, mu_earth);
  }

  if (config_.use_sph) {
    if (!field_) {
      return GravitySphResult{.status = astroforces::core::Status::DataUnavailable};
    }

    const int nmax = std::min(config_.max_degree, field_->max_degree);
    if (nmax >= 2) {
      const bool needs_coeff_mutation = field_->has_time_variable || config_.use_solid_earth_tides || config_.use_pole_tide_solid
                                        || config_.use_pole_tide_ocean || config_.use_aod || config_.use_ocean_tide
                                        || config_.use_atmos_tide;
      const Eigen::MatrixXd* coeff_c = &field_->C;
      const Eigen::MatrixXd* coeff_s = &field_->S;
      Eigen::MatrixXd Ceff;
      Eigen::MatrixXd Seff;
      if (needs_coeff_mutation) {
        Ceff = field_->C;
        Seff = field_->S;
        coeff_c = &Ceff;
        coeff_s = &Seff;
        if (field_->has_time_variable) {
          apply_time_variable_coefficients(*field_, jd_utc, nmax, Ceff, Seff);
        }
      }

      Eigen::MatrixXd tmpC;
      Eigen::MatrixXd tmpS;
      if (needs_coeff_mutation) {
        tmpC = Eigen::MatrixXd::Zero(Ceff.rows(), Ceff.cols());
        tmpS = Eigen::MatrixXd::Zero(Seff.rows(), Seff.cols());
      }

      if (config_.use_solid_earth_tides && (config_.use_sun_tide || config_.use_moon_tide || config_.use_solid_earth_tide2)) {
        if ((config_.use_sun_tide || config_.use_moon_tide) && (!ephemeris_ || !workspace_)) {
          return GravitySphResult{.status = astroforces::core::Status::DataUnavailable};
        }

        Eigen::MatrixXd dC = Eigen::MatrixXd::Zero(Ceff.rows(), Ceff.cols());
        Eigen::MatrixXd dS = Eigen::MatrixXd::Zero(Seff.rows(), Seff.cols());

        if (config_.use_sun_tide) {
          const auto sun = ephemeris_->PlephSi(jd_utc, jpl::eph::Body::Sun, jpl::eph::Body::Earth, false, *workspace_);
          if (!sun.has_value()) {
            return GravitySphResult{.status = map_jpl_error(sun.error())};
          }
          const auto r_sun_ecef = rot_z(gmst, to_vec3(sun.value().pv));
          tides::add_solid_earth_tide1_delta(r_sun_ecef, config_.mu_sun_m3_s2, mu_earth, radius_m, nmax, dC, dS);
          out.solid_tide_sun_mps2 = accel_sph_noncentral(r_ecef, dC, dS, nmax, mu_earth, radius_m);
        }

        if (config_.use_moon_tide) {
          const auto moon = ephemeris_->PlephSi(jd_utc, jpl::eph::Body::Moon, jpl::eph::Body::Earth, false, *workspace_);
          if (!moon.has_value()) {
            return GravitySphResult{.status = map_jpl_error(moon.error())};
          }
          tmpC.setZero();
          tmpS.setZero();
          const auto r_moon_ecef = rot_z(gmst, to_vec3(moon.value().pv));
          tides::add_solid_earth_tide1_delta(r_moon_ecef, config_.mu_moon_m3_s2, mu_earth, radius_m, nmax, tmpC, tmpS);
          out.solid_tide_moon_mps2 = accel_sph_noncentral(r_ecef, tmpC, tmpS, nmax, mu_earth, radius_m);
          dC += tmpC;
          dS += tmpS;
        }

        if (config_.use_solid_earth_tide2) {
          tmpC.setZero();
          tmpS.setZero();
          tides::add_solid_earth_tide2_delta(jd_utc, gmst, tmpC, tmpS);
          out.solid_tide_freqdep_mps2 = accel_sph_noncentral(r_ecef, tmpC, tmpS, nmax, mu_earth, radius_m);
          dC += tmpC;
          dS += tmpS;
        }

        Ceff += dC;
        Seff += dS;
      }

      if (config_.use_pole_tide_solid || config_.use_pole_tide_ocean) {
        if (!eop_) {
          return GravitySphResult{.status = astroforces::core::Status::DataUnavailable};
        }
        const auto eop_now = eop_->at_utc_seconds(state.epoch.utc_seconds);
        if (!eop_now.has_value()) {
          return GravitySphResult{.status = astroforces::core::Status::DataUnavailable};
        }

        constexpr double kApproxTtMinusUtcSeconds = 69.184;
        const double jd_tt = jd_utc + kApproxTtMinusUtcSeconds / astroforces::core::constants::kSecondsPerDay;
        const double mjd_tt = jd_tt - 2400000.5;

        if (config_.use_pole_tide_solid) {
          tmpC.setZero();
          tmpS.setZero();
          tides::add_pole_solid_tide_delta(mjd_tt, eop_now->xp_rad, eop_now->yp_rad, tmpC, tmpS);
          out.pole_tide_solid_mps2 = accel_sph_noncentral(r_ecef, tmpC, tmpS, nmax, mu_earth, radius_m);
          Ceff += tmpC;
          Seff += tmpS;
        }

        if (config_.use_pole_tide_ocean) {
          tmpC.setZero();
          tmpS.setZero();
          if (ocean_pole_tide_) {
            double xpv_mas = 0.0;
            double ypv_mas = 0.0;
            const double t_years = (mjd_tt - 51544.5) / 365.25;
            xpv_mas = 55.0 + 1.677 * t_years;
            ypv_mas = 320.5 + 3.460 * t_years;
            const double m1 = +(eop_now->xp_rad / astroforces::core::constants::kArcsecToRad - xpv_mas / 1000.0);
            const double m2 = -(eop_now->yp_rad / astroforces::core::constants::kArcsecToRad - ypv_mas / 1000.0);
            ocean_pole_tide_->add_delta_coefficients(m1, m2, tmpC, tmpS);
          } else {
            tides::add_pole_ocean_tide_delta(mjd_tt, eop_now->xp_rad, eop_now->yp_rad, tmpC, tmpS);
          }
          out.pole_tide_ocean_mps2 = accel_sph_noncentral(r_ecef, tmpC, tmpS, nmax, mu_earth, radius_m);
          Ceff += tmpC;
          Seff += tmpS;
        }
      }

      if (config_.use_aod) {
        if (!aod_tide_) {
          return GravitySphResult{.status = astroforces::core::Status::DataUnavailable};
        }
        tmpC.setZero();
        tmpS.setZero();
        if (aod_tide_->interpolate_delta_coefficients(jd_utc, tmpC, tmpS)) {
          out.aod_mps2 = accel_sph_noncentral(r_ecef, tmpC, tmpS, nmax, mu_earth, radius_m);
          Ceff += tmpC;
          Seff += tmpS;
        }
      }

      if (config_.use_ocean_tide) {
        if (!ocean_tide_) {
          return GravitySphResult{.status = astroforces::core::Status::DataUnavailable};
        }
        tmpC.setZero();
        tmpS.setZero();
        ocean_tide_->add_delta_coefficients(jd_utc, gmst, tmpC, tmpS);
        out.ocean_tide_mps2 = accel_sph_noncentral(r_ecef, tmpC, tmpS, nmax, mu_earth, radius_m);
        Ceff += tmpC;
        Seff += tmpS;
      }

      if (config_.use_atmos_tide) {
        if (!atmos_tide_) {
          return GravitySphResult{.status = astroforces::core::Status::DataUnavailable};
        }
        tmpC.setZero();
        tmpS.setZero();
        atmos_tide_->add_delta_coefficients(jd_utc, gmst, tmpC, tmpS);
        out.atmos_tide_mps2 = accel_sph_noncentral(r_ecef, tmpC, tmpS, nmax, mu_earth, radius_m);
        Ceff += tmpC;
        Seff += tmpS;
      }

      out.sph_mps2 = accel_sph_noncentral(r_ecef, *coeff_c, *coeff_s, nmax, mu_earth, radius_m);
    }
  }

  out.acceleration_mps2 = out.central_mps2 + out.sph_mps2;

  if (state.frame == astroforces::core::Frame::ECI) {
    out.acceleration_mps2 = rot_z(-gmst, out.acceleration_mps2);
    out.central_mps2 = rot_z(-gmst, out.central_mps2);
    out.sph_mps2 = rot_z(-gmst, out.sph_mps2);
    out.solid_tide_sun_mps2 = rot_z(-gmst, out.solid_tide_sun_mps2);
    out.solid_tide_moon_mps2 = rot_z(-gmst, out.solid_tide_moon_mps2);
    out.solid_tide_freqdep_mps2 = rot_z(-gmst, out.solid_tide_freqdep_mps2);
    out.pole_tide_solid_mps2 = rot_z(-gmst, out.pole_tide_solid_mps2);
    out.pole_tide_ocean_mps2 = rot_z(-gmst, out.pole_tide_ocean_mps2);
    out.aod_mps2 = rot_z(-gmst, out.aod_mps2);
    out.ocean_tide_mps2 = rot_z(-gmst, out.ocean_tide_mps2);
    out.atmos_tide_mps2 = rot_z(-gmst, out.atmos_tide_mps2);
  }

  if (!finite_vec(out.acceleration_mps2)) {
    return GravitySphResult{.status = astroforces::core::Status::NumericalError};
  }

  out.status = astroforces::core::Status::Ok;
  return out;
}

}  // namespace astroforces::forces
