/**
 * @file perturbation_profile_cli.cpp
 * @brief Altitude sweep utility for perturbation-acceleration profiling.
 * @author Watosn
 */

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/adapters/dtm2020_adapter.hpp"
#include "astroforces/atmo/conversions.hpp"
#include "astroforces/atmo/constants.hpp"
#include "astroforces/drag/drag_perturbation.hpp"
#include "astroforces/forces/erp_perturbation.hpp"
#include "astroforces/forces/perturbation.hpp"
#include "astroforces/forces/relativity_perturbation.hpp"
#include "astroforces/forces/relativity_model.hpp"
#include "astroforces/forces/third_body.hpp"
#include "astroforces/models/exponential_atmosphere.hpp"
#include "astroforces/sc/spacecraft.hpp"
#include "astroforces/srp/srp_perturbation.hpp"
#include "astroforces/weather/celestrak_csv_provider.hpp"

namespace {

constexpr double kDtmOperationalMinAltKm = 120.0;
constexpr double kDtmOperationalMaxAltKm = 1500.0;
constexpr const char* kRequiredDtmCoeffFilename = "DTM_2020_F107_Kp.dat";

bool is_required_dtm_coeff_file(const std::filesystem::path& path) { return path.filename() == kRequiredDtmCoeffFilename; }

double circular_speed_mps(double radius_m) {
  return std::sqrt(astroforces::atmo::constants::kEarthMuM3S2 / radius_m);
}

double magnitude(const astroforces::atmo::Vec3& v) { return astroforces::atmo::norm(v); }

astroforces::atmo::Vec3 ecef_to_eci(const astroforces::atmo::Vec3& r_ecef_m, double utc_seconds) {
  constexpr double kPi = 3.1415926535897932384626433832795;
  const double jd = astroforces::atmo::utc_seconds_to_julian_date_utc(utc_seconds);
  const double t = (jd - 2451545.0) / 36525.0;
  double gmst_deg = 280.46061837 + 360.98564736629 * (jd - 2451545.0) + 0.000387933 * t * t - (t * t * t) / 38710000.0;
  gmst_deg = std::fmod(gmst_deg, 360.0);
  if (gmst_deg < 0.0) {
    gmst_deg += 360.0;
  }
  const double th = gmst_deg * kPi / 180.0;
  const double c = std::cos(th);
  const double s = std::sin(th);
  return astroforces::atmo::Vec3{
      c * r_ecef_m.x - s * r_ecef_m.y,
      s * r_ecef_m.x + c * r_ecef_m.y,
      r_ecef_m.z,
  };
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 7 || argc > 9) {
    spdlog::error(
        "usage: perturbation_profile_cli <output_csv> <alt_min_km> <alt_max_km> <samples> <dtm_coeff_file> <space_weather_csv> [jpl_ephemeris_file] [epoch_utc_s]");
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

  struct ComponentModel {
    std::string label{};
    astroforces::atmo::Frame frame{astroforces::atmo::Frame::ECI};
    std::unique_ptr<astroforces::forces::IPerturbationModel> model{};
    bool has_valid_altitude_band{};
    double min_alt_km{};
    double max_alt_km{};
  };

  std::vector<ComponentModel> components{};
  components.push_back(ComponentModel{
      .label = "drag",
      .frame = astroforces::atmo::Frame::ECEF,
      .model = std::make_unique<astroforces::drag::DragPerturbationModel>(*weather, *atmosphere, wind, &sc, "drag"),
      .has_valid_altitude_band = true,
      .min_alt_km = kDtmOperationalMinAltKm,
      .max_alt_km = kDtmOperationalMaxAltKm,
  });
  components.push_back(ComponentModel{
      .label = "erp",
      .frame = astroforces::atmo::Frame::ECI,
      .model = std::make_unique<astroforces::erp::ErpPerturbationModel>(astroforces::erp::ErpAccelerationModel{}, &sc, "erp"),
  });
  components.push_back(ComponentModel{
      .label = "relativity",
      .frame = astroforces::atmo::Frame::ECI,
      .model = std::make_unique<astroforces::forces::RelativityPerturbationModel>(
          astroforces::forces::RelativityAccelerationModel::Create(
              {.ephemeris_file = std::filesystem::path(eph_file), .use_geodesic_precession = !eph_file.empty()}),
          "relativity"),
  });

  if (!eph_file.empty()) {
    components.push_back(ComponentModel{
        .label = "srp",
        .frame = astroforces::atmo::Frame::ECI,
        .model = std::make_unique<astroforces::srp::SrpPerturbationModel>(
            astroforces::srp::SrpAccelerationModel::Create({.ephemeris_file = std::filesystem::path(eph_file), .use_eclipse = false}),
            &sc,
            "srp"),
    });
    components.push_back(ComponentModel{
        .label = "third_body_sun",
        .frame = astroforces::atmo::Frame::ECI,
        .model = astroforces::forces::ThirdBodyPerturbationModel::Create(
            {.ephemeris_file = std::filesystem::path(eph_file), .use_sun = true, .use_moon = false, .name = "third_body_sun"}),
    });
    components.push_back(ComponentModel{
        .label = "third_body_moon",
        .frame = astroforces::atmo::Frame::ECI,
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
  for (const auto& comp : components) {
    out << "," << comp.label << "_mps2";
  }
  out << ",total_mps2,status\n";

  for (int i = 0; i < samples; ++i) {
    const double u = static_cast<double>(i) / static_cast<double>(samples - 1);
    const double alt_km = alt_min_km + (alt_max_km - alt_min_km) * u;
    const double r_m = astroforces::atmo::constants::kEarthRadiusWgs84M + alt_km * 1000.0;

    const astroforces::atmo::Vec3 r_ecef_m{r_m, 0.0, 0.0};
    const astroforces::atmo::Vec3 v_ecef_mps{0.0, circular_speed_mps(r_m), 0.0};

    astroforces::atmo::StateVector drag_state{};
    drag_state.epoch.utc_seconds = epoch_utc_s;
    drag_state.frame = astroforces::atmo::Frame::ECEF;
    drag_state.position_m = r_ecef_m;
    drag_state.velocity_mps = v_ecef_mps;

    astroforces::atmo::StateVector third_state{};
    third_state.epoch.utc_seconds = epoch_utc_s;
    third_state.frame = astroforces::atmo::Frame::ECI;
    third_state.position_m = ecef_to_eci(r_ecef_m, epoch_utc_s);
    third_state.velocity_mps = astroforces::atmo::Vec3{};

    std::vector<double> comp_mags{};
    comp_mags.reserve(components.size());
    astroforces::atmo::Status status = astroforces::atmo::Status::Ok;
    double rss_sum = 0.0;

    for (const auto& comp : components) {
      if (comp.has_valid_altitude_band && (alt_km < comp.min_alt_km || alt_km > comp.max_alt_km)) {
        comp_mags.push_back(0.0);
        continue;
      }
      const auto& state = (comp.frame == astroforces::atmo::Frame::ECEF) ? drag_state : third_state;
      const auto c = comp.model->evaluate(astroforces::forces::PerturbationRequest{.state = state, .spacecraft = &sc});
      const double m = magnitude(c.acceleration_mps2);
      comp_mags.push_back(m);
      rss_sum += m * m;
      if (status == astroforces::atmo::Status::Ok && c.status != astroforces::atmo::Status::Ok) {
        status = c.status;
      }
    }
    out << fmt::format("{:.6f}", alt_km);
    for (const double m : comp_mags) {
      out << fmt::format(",{:.12e}", m);
    }
    out << fmt::format(",{:.12e},{}\n", std::sqrt(rss_sum), static_cast<int>(status));
  }

  spdlog::info("wrote perturbation profile: {}", out_csv.string());
  return 0;
}
