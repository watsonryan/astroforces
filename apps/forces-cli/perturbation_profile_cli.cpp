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
#include "astroforces/core/transforms.hpp"
#include "astroforces/core/constants.hpp"
#include "astroforces/forces/surface/antenna_thrust/antenna_thrust_perturbation.hpp"
#include "astroforces/forces/surface/drag/drag_perturbation.hpp"
#include "astroforces/forces/surface/earth_radiation/earth_radiation_perturbation.hpp"
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

}  // namespace

int main(int argc, char** argv) {
  if (argc < 7 || argc > 17) {
    spdlog::error(
        "usage: perturbation_profile_cli <output_csv> <alt_min_km> <alt_max_km> <samples> <dtm_coeff_file> <space_weather_csv> [jpl_ephemeris_file] [epoch_utc_s] [gravity_gfc_file] [gravity_max_degree] [eop_finals_file] [ocean_tide_file] [atmos_tide_file] [ocean_pole_tide_file] [aod_file] [cip_xys_file]");
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
  const std::filesystem::path eop_finals_file = (argc >= 12) ? std::filesystem::path(argv[11]) : std::filesystem::path{};
  const std::filesystem::path ocean_tide_file = (argc >= 13) ? std::filesystem::path(argv[12]) : std::filesystem::path{};
  const std::filesystem::path atmos_tide_file = (argc >= 14) ? std::filesystem::path(argv[13]) : std::filesystem::path{};
  const std::filesystem::path ocean_pole_tide_file = (argc >= 15) ? std::filesystem::path(argv[14]) : std::filesystem::path{};
  const std::filesystem::path aod_file = (argc >= 16) ? std::filesystem::path(argv[15]) : std::filesystem::path{};
  const std::filesystem::path cip_xys_file = (argc >= 17) ? std::filesystem::path(argv[16]) : std::filesystem::path{};

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
  if (!eop_finals_file.empty() && !std::filesystem::exists(eop_finals_file)) {
    spdlog::error("EOP finals file not found: {}", eop_finals_file.string());
    return 10;
  }
  if (!ocean_tide_file.empty() && !std::filesystem::exists(ocean_tide_file)) {
    spdlog::error("ocean tide file not found: {}", ocean_tide_file.string());
    return 11;
  }
  if (!atmos_tide_file.empty() && !std::filesystem::exists(atmos_tide_file)) {
    spdlog::error("atmos tide file not found: {}", atmos_tide_file.string());
    return 12;
  }
  if (!ocean_pole_tide_file.empty() && !std::filesystem::exists(ocean_pole_tide_file)) {
    spdlog::error("ocean pole tide file not found: {}", ocean_pole_tide_file.string());
    return 13;
  }
  if (!aod_file.empty() && !std::filesystem::exists(aod_file)) {
    spdlog::error("aod file not found: {}", aod_file.string());
    return 14;
  }
  if (!cip_xys_file.empty() && !std::filesystem::exists(cip_xys_file)) {
    spdlog::error("CIP X/Y/s file not found: {}", cip_xys_file.string());
    return 15;
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
    const bool has_ephem = !eph_file.empty();
    const bool has_eop = !eop_finals_file.empty();
    const bool has_cip = !cip_xys_file.empty();
    const bool has_eop_cip = has_eop && has_cip;
    const bool enable_lunisolar_tides = has_ephem && has_eop_cip;
    const bool enable_tide2 = has_eop_cip;
    const bool enable_pole_tides = has_eop_cip;
    const bool enable_ocean_tide = has_eop_cip && !ocean_tide_file.empty();
    const bool enable_atmos_tide = has_eop_cip && !atmos_tide_file.empty();
    const bool enable_aod = has_eop_cip && !aod_file.empty();
    if (has_ephem && !has_eop_cip) {
      spdlog::warn(
          "gravity profile: ephemeris provided but EOP/CIP missing; disabling Sun/Moon "
          "and tide terms that require EOP/CIP");
    }
    gravity = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_gfc_file,
         .ephemeris_file = std::filesystem::path(eph_file),
         .eop_finals_file = eop_finals_file,
         .cip_xys_file = cip_xys_file,
         .ocean_pole_tide_file = ocean_pole_tide_file,
         .aod_file = aod_file,
         .ocean_tide_file = ocean_tide_file,
         .atmos_tide_file = atmos_tide_file,
         .max_degree = gravity_max_degree,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = enable_lunisolar_tides || enable_tide2,
         .use_sun_tide = enable_lunisolar_tides,
         .use_moon_tide = enable_lunisolar_tides,
         .use_solid_earth_tide2 = enable_tide2,
         .use_pole_tide_solid = enable_pole_tides,
         .use_pole_tide_ocean = enable_pole_tides,
         .use_ocean_tide = enable_ocean_tide,
         .use_atmos_tide = enable_atmos_tide,
         .use_aod = enable_aod});
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
      .model = std::make_unique<astroforces::forces::DragPerturbationModel>(*weather, *atmosphere, wind, &sc, "drag"),
      .has_valid_altitude_band = true,
      .min_alt_km = kDtmOperationalMinAltKm,
      .max_alt_km = kDtmOperationalMaxAltKm,
  });
  components.push_back(ComponentModel{
      .label = "antenna_thrust",
      .frame = astroforces::core::Frame::ECI,
      .model = std::make_unique<astroforces::forces::AntennaThrustPerturbationModel>(
          astroforces::forces::AntennaThrustAccelerationModel({
              .transmit_power_w = 20.0,
              .efficiency = 1.0,
              .direction_mode = astroforces::forces::AntennaThrustDirectionMode::Velocity,
          }),
          &sc,
          "antenna_thrust"),
  });
  components.push_back(ComponentModel{
      .label = "earth_radiation",
      .frame = astroforces::core::Frame::ECI,
      .model = std::make_unique<astroforces::forces::EarthRadiationPerturbationModel>(
          astroforces::forces::EarthRadiationAccelerationModel({.ephemeris_file = std::filesystem::path(eph_file)}), &sc, "earth_radiation"),
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
        .model = std::make_unique<astroforces::forces::SrpPerturbationModel>(
            astroforces::forces::SrpAccelerationModel::Create({.ephemeris_file = std::filesystem::path(eph_file), .use_eclipse = false}),
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
    out << ",gravity_central_mps2,gravity_sph_tides_mps2,gravity_tide_solid_sun_mps2,gravity_tide_solid_moon_mps2,gravity_tide_solid_freqdep_mps2,gravity_tide_pole_solid_mps2,gravity_tide_pole_ocean_mps2,gravity_tide_aod_mps2,gravity_tide_ocean_mps2,gravity_tide_atmos_mps2";
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
    third_state.position_m = astroforces::core::Vec3{r_m, 0.0, 0.0};
    third_state.velocity_mps = astroforces::core::Vec3{0.0, circular_speed_mps(r_m), 0.0};

    std::vector<double> gravity_mags{};
    if (gravity) {
      gravity_mags.assign(10, std::numeric_limits<double>::quiet_NaN());
    }
    std::vector<double> comp_mags{};
    comp_mags.reserve(components.size());
    astroforces::core::Status status = astroforces::core::Status::Ok;
    double rss_sum = 0.0;
    if (gravity) {
      const auto g = gravity->evaluate(drag_state);
      const double g_c = magnitude(g.central_mps2);
      const double g_sph = magnitude(g.sph_mps2);
      const double g_tss = magnitude(g.solid_tide_sun_mps2);
      const double g_tsm = magnitude(g.solid_tide_moon_mps2);
      const double g_tsf = magnitude(g.solid_tide_freqdep_mps2);
      const double g_tps = magnitude(g.pole_tide_solid_mps2);
      const double g_tpo = magnitude(g.pole_tide_ocean_mps2);
      const double g_aod = magnitude(g.aod_mps2);
      const double g_to = magnitude(g.ocean_tide_mps2);
      const double g_ta = magnitude(g.atmos_tide_mps2);
      gravity_mags[0] = g_c;
      gravity_mags[1] = g_sph;
      gravity_mags[2] = g_tss;
      gravity_mags[3] = g_tsm;
      gravity_mags[4] = g_tsf;
      gravity_mags[5] = g_tps;
      gravity_mags[6] = g_tpo;
      gravity_mags[7] = g_aod;
      gravity_mags[8] = g_to;
      gravity_mags[9] = g_ta;
      rss_sum += g_c * g_c + g_sph * g_sph + g_tss * g_tss + g_tsm * g_tsm + g_tsf * g_tsf + g_tps * g_tps + g_tpo * g_tpo
                 + g_aod * g_aod + g_to * g_to + g_ta * g_ta;
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
