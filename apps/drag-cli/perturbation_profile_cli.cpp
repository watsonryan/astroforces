/**
 * @file perturbation_profile_cli.cpp
 * @brief Altitude sweep utility for perturbation-acceleration profiling.
 * @author Watosn
 */

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/adapters/dtm2020_adapter.hpp"
#include "astroforces/atmo/conversions.hpp"
#include "astroforces/atmo/constants.hpp"
#include "astroforces/forces/surface/drag/drag_perturbation.hpp"
#include "astroforces/forces/surface/erp/erp_perturbation.hpp"
#include "astroforces/forces/gravity/gravity_sph_model.hpp"
#include "astroforces/forces/core/perturbation.hpp"
#include "astroforces/forces/gravity/relativity_perturbation.hpp"
#include "astroforces/forces/gravity/relativity_model.hpp"
#include "astroforces/forces/gravity/third_body.hpp"
#include "astroforces/models/exponential_atmosphere.hpp"
#include "astroforces/sc/spacecraft.hpp"
#include "astroforces/forces/surface/srp/srp_perturbation.hpp"
#include "astroforces/weather/celestrak_csv_provider.hpp"

namespace {

constexpr double kDtmOperationalMinAltKm = 120.0;
constexpr double kDtmOperationalMaxAltKm = 1500.0;
constexpr const char* kRequiredDtmCoeffFilename = "DTM_2020_F107_Kp.dat";

bool is_required_dtm_coeff_file(const std::filesystem::path& path) { return path.filename() == kRequiredDtmCoeffFilename; }

double circular_speed_mps(double radius_m) {
  return std::sqrt(astroforces::core::constants::kEarthMuM3S2 / radius_m);
}

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

astroforces::core::Vec3 ecef_to_eci(const astroforces::core::Vec3& r_ecef_m, double utc_seconds) {
  constexpr double kPi = 3.1415926535897932384626433832795;
  const double jd = astroforces::core::utc_seconds_to_julian_date_utc(utc_seconds);
  const double t = (jd - 2451545.0) / 36525.0;
  double gmst_deg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * t * t - (t * t * t) / 38710000.0;
  gmst_deg = std::fmod(gmst_deg, 360.0);
  if (gmst_deg < 0.0) {
    gmst_deg += 360.0;
  }
  const double th = gmst_deg * kPi / 180.0;
  const double c = std::cos(th);
  const double s = std::sin(th);
  return astroforces::core::Vec3{
      c * r_ecef_m.x - s * r_ecef_m.y,
      s * r_ecef_m.x + c * r_ecef_m.y,
      r_ecef_m.z,
  };
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 7 || argc > 11) {
    spdlog::error(
        "usage: perturbation_profile_cli <output_csv> <alt_min_km> <alt_max_km> <samples> <dtm_coeff_file> <space_weather_csv> [jpl_ephemeris_file] [epoch_utc_s] [gravity_gfc_file] [gravity_max_degree]");
    return 1;
  }

  const std::filesystem::path out_csv = argv[1];
  const double alt_min_km = std::atof(argv[2]);
  const double alt_max_km = std::atof(argv[3]);
  const int samples = std::atoi(argv[4]);
  const std::filesystem::path dtm_coeff_file = argv[5];
  const std::filesystem::path space_weather_csv = argv[6];
  const std::string eph_file = (argc >= 8) ? argv[7] : "";
  const double epoch_utc_s = (argc >= 9) ? std::atof(argv[8]) : 1.0e9;
  const std::filesystem::path gravity_gfc_file = (argc >= 10) ? std::filesystem::path(argv[9]) : std::filesystem::path{};
  const int gravity_max_degree = (argc >= 11) ? std::atoi(argv[10]) : 120;

  if (!(alt_min_km >= 0.0) || !(alt_max_km > alt_min_km) || samples < 2) {
    spdlog::error("invalid sweep parameters: require alt_min>=0, alt_max>alt_min, samples>=2");
    return 2;
  }
  if (!std::filesystem::exists(dtm_coeff_file)) {
    spdlog::error("dtm coeff file not found: {}", dtm_coeff_file.string());
    return 4;
  }
  if (!is_required_dtm_coeff_file(dtm_coeff_file)) {
    spdlog::error("dtm coeff file must be {} (received: {})", kRequiredDtmCoeffFilename, dtm_coeff_file.filename().string());
    return 8;
  }
  if (!std::filesystem::exists(space_weather_csv)) {
    spdlog::error("space weather csv not found: {}", space_weather_csv.string());
    return 5;
  }
  if (!gravity_gfc_file.empty() && !std::filesystem::exists(gravity_gfc_file)) {
    spdlog::error("gravity gfc file not found: {}", gravity_gfc_file.string());
    return 9;
  }

  std::ofstream out(out_csv);
  if (!out) {
    spdlog::error("failed to open output csv: {}", out_csv.string());
    return 3;
  }

  const auto weather = astroforces::weather::CelesTrakCsvSpaceWeatherProvider::Create(
      astroforces::weather::CelesTrakCsvSpaceWeatherProvider::Config{.csv_file = space_weather_csv});
  const auto atmosphere = astroforces::adapters::Dtm2020AtmosphereAdapter::Create(
      astroforces::adapters::Dtm2020AtmosphereAdapter::Config{.coeff_file = dtm_coeff_file});
  astroforces::models::ZeroWindModel wind;
  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.25, .cr = 1.3, .use_surface_model = false, .surfaces = {}};
  if (!weather) {
    spdlog::error("failed to initialize CelesTrak weather provider");
    return 6;
  }
  if (!atmosphere) {
    spdlog::error("failed to initialize DTM2020 adapter");
    return 7;
  }
  std::unique_ptr<astroforces::forces::GravitySphAccelerationModel> gravity{};
  if (!gravity_gfc_file.empty()) {
    gravity = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_gfc_file,
         .ephemeris_file = std::filesystem::path(eph_file),
         .max_degree = gravity_max_degree,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = !eph_file.empty()});
  }

  struct ComponentModel {
    std::string label{};
    astroforces::core::Frame frame{astroforces::core::Frame::ECI};
    std::unique_ptr<astroforces::forces::IPerturbationModel> model{};
    bool has_valid_altitude_band{};
    double min_alt_km{};
    double max_alt_km{};
  };

  std::vector<ComponentModel> components{};
  components.push_back(ComponentModel{
      .label = "drag",
      .frame = astroforces::core::Frame::ECEF,
      .model = std::make_unique<astroforces::drag::DragPerturbationModel>(*weather, *atmosphere, wind, &sc, "drag"),
      .has_valid_altitude_band = true,
      .min_alt_km = kDtmOperationalMinAltKm,
      .max_alt_km = kDtmOperationalMaxAltKm,
  });
  components.push_back(ComponentModel{
      .label = "erp",
      .frame = astroforces::core::Frame::ECI,
      .model = std::make_unique<astroforces::erp::ErpPerturbationModel>(astroforces::erp::ErpAccelerationModel{}, &sc, "erp"),
  });
  components.push_back(ComponentModel{
      .label = "relativity",
      .frame = astroforces::core::Frame::ECI,
      .model = std::make_unique<astroforces::forces::RelativityPerturbationModel>(
          astroforces::forces::RelativityAccelerationModel::Create(
              {.ephemeris_file = std::filesystem::path(eph_file), .use_geodesic_precession = !eph_file.empty()}),
          "relativity"),
  });

  if (!eph_file.empty()) {
    components.push_back(ComponentModel{
        .label = "srp",
        .frame = astroforces::core::Frame::ECI,
        .model = std::make_unique<astroforces::srp::SrpPerturbationModel>(
            astroforces::srp::SrpAccelerationModel::Create({.ephemeris_file = std::filesystem::path(eph_file), .use_eclipse = false}),
            &sc,
            "srp"),
    });
    components.push_back(ComponentModel{
        .label = "third_body_sun",
        .frame = astroforces::core::Frame::ECI,
        .model = astroforces::forces::ThirdBodyPerturbationModel::Create(
            {.ephemeris_file = std::filesystem::path(eph_file), .use_sun = true, .use_moon = false, .name = "third_body_sun"}),
    });
    components.push_back(ComponentModel{
        .label = "third_body_moon",
        .frame = astroforces::core::Frame::ECI,
        .model = astroforces::forces::ThirdBodyPerturbationModel::Create(
            {.ephemeris_file = std::filesystem::path(eph_file), .use_sun = false, .use_moon = true, .name = "third_body_moon"}),
    });
  } else {
    spdlog::warn("no ephemeris path provided; third-body component curves will be omitted");
  }
  spdlog::warn("drag component is limited to DTM operational validity band: [{}, {}] km",
               kDtmOperationalMinAltKm,
               kDtmOperationalMaxAltKm);

  out << "altitude_km";
  if (gravity) {
    out << ",gravity_central_mps2,gravity_sph_mps2,gravity_tide_sun_mps2,gravity_tide_moon_mps2";
  }
  for (const auto& comp : components) {
    out << "," << comp.label << "_mps2";
  }
  out << ",total_mps2,status\n";

  for (int i = 0; i < samples; ++i) {
    const double u = static_cast<double>(i) / static_cast<double>(samples - 1);
    const double alt_km = alt_min_km + (alt_max_km - alt_min_km) * u;
    const double r_m = astroforces::core::constants::kEarthRadiusWgs84M + alt_km * 1000.0;

    const astroforces::core::Vec3 r_ecef_m{r_m, 0.0, 0.0};
    const astroforces::core::Vec3 v_ecef_mps{0.0, circular_speed_mps(r_m), 0.0};

    astroforces::core::StateVector drag_state{};
    drag_state.epoch.utc_seconds = epoch_utc_s;
    drag_state.frame = astroforces::core::Frame::ECEF;
    drag_state.position_m = r_ecef_m;
    drag_state.velocity_mps = v_ecef_mps;

    astroforces::core::StateVector third_state{};
    third_state.epoch.utc_seconds = epoch_utc_s;
    third_state.frame = astroforces::core::Frame::ECI;
    third_state.position_m = ecef_to_eci(r_ecef_m, epoch_utc_s);
    third_state.velocity_mps = astroforces::core::Vec3{};

    std::vector<double> gravity_mags{};
    if (gravity) {
      gravity_mags.reserve(4);
    }
    std::vector<double> comp_mags{};
    comp_mags.reserve(components.size());
    astroforces::core::Status status = astroforces::core::Status::Ok;
    double rss_sum = 0.0;
    if (gravity) {
      const auto g = gravity->evaluate(drag_state);
      const double g_c = magnitude(g.central_mps2);
      const double g_sph = magnitude(g.sph_mps2);
      const double g_ts = magnitude(g.solid_tide_sun_mps2);
      const double g_tm = magnitude(g.solid_tide_moon_mps2);
      gravity_mags.push_back(g_c);
      gravity_mags.push_back(g_sph);
      gravity_mags.push_back(g_ts);
      gravity_mags.push_back(g_tm);
      rss_sum += g_c * g_c + g_sph * g_sph + g_ts * g_ts + g_tm * g_tm;
      if (status == astroforces::core::Status::Ok && g.status != astroforces::core::Status::Ok) {
        status = g.status;
      }
    }

    for (const auto& comp : components) {
      if (comp.has_valid_altitude_band && (alt_km < comp.min_alt_km || alt_km > comp.max_alt_km)) {
        comp_mags.push_back(std::numeric_limits<double>::quiet_NaN());
        continue;
      }
      const auto& state = (comp.frame == astroforces::core::Frame::ECEF) ? drag_state : third_state;
      const auto c = comp.model->evaluate(astroforces::forces::PerturbationRequest{.state = state, .spacecraft = &sc});
      const double m = magnitude(c.acceleration_mps2);
      comp_mags.push_back(m);
      rss_sum += m * m;
      if (status == astroforces::core::Status::Ok && c.status != astroforces::core::Status::Ok) {
        status = c.status;
      }
    }
    out << fmt::format("{:.6f}", alt_km);
    for (const double m : gravity_mags) {
      out << fmt::format(",{:.12e}", m);
    }
    for (const double m : comp_mags) {
      out << fmt::format(",{:.12e}", m);
    }
    out << fmt::format(",{:.12e},{}\n", std::sqrt(rss_sum), static_cast<int>(status));
  }

  spdlog::info("wrote perturbation profile: {}", out_csv.string());
  return 0;
}
