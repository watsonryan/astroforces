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

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "dragcpp/atmo/constants.hpp"
#include "dragcpp/drag/drag_perturbation.hpp"
#include "dragcpp/forces/perturbation.hpp"
#include "dragcpp/forces/third_body.hpp"
#include "dragcpp/models/exponential_atmosphere.hpp"
#include "dragcpp/sc/spacecraft.hpp"
#include "dragcpp/weather/static_provider.hpp"

namespace {

double circular_speed_mps(double radius_m) {
  return std::sqrt(astroforces::atmo::constants::kEarthMuM3S2 / radius_m);
}

double magnitude(const astroforces::atmo::Vec3& v) { return astroforces::atmo::norm(v); }

}  // namespace

int main(int argc, char** argv) {
  if (argc < 5 || argc > 7) {
    spdlog::error("usage: perturbation_profile_cli <output_csv> <alt_min_km> <alt_max_km> <samples> [jpl_ephemeris_file] [epoch_utc_s]");
    return 1;
  }

  const std::filesystem::path out_csv = argv[1];
  const double alt_min_km = std::atof(argv[2]);
  const double alt_max_km = std::atof(argv[3]);
  const int samples = std::atoi(argv[4]);
  const std::string eph_file = (argc >= 6) ? argv[5] : "";
  const double epoch_utc_s = (argc >= 7) ? std::atof(argv[6]) : 1.0e9;

  if (!(alt_min_km >= 0.0) || !(alt_max_km > alt_min_km) || samples < 2) {
    spdlog::error("invalid sweep parameters: require alt_min>=0, alt_max>alt_min, samples>=2");
    return 2;
  }

  std::ofstream out(out_csv);
  if (!out) {
    spdlog::error("failed to open output csv: {}", out_csv.string());
    return 3;
  }

  const astroforces::atmo::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0, .status = astroforces::atmo::Status::Ok};
  astroforces::weather::StaticSpaceWeatherProvider weather(wx);
  astroforces::models::ExponentialAtmosphereModel atmosphere(3.0e-11, 400e3, 65e3, 900.0);
  astroforces::models::ZeroWindModel wind;
  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.25, .use_surface_model = false, .surfaces = {}};

  astroforces::forces::PerturbationStack stack;
  stack.add(std::make_unique<astroforces::drag::DragPerturbationModel>(weather, atmosphere, wind, &sc, "drag"));

  bool third_body_enabled = false;
  if (!eph_file.empty()) {
    auto third_body = astroforces::forces::ThirdBodyPerturbationModel::Create(
        {.ephemeris_file = std::filesystem::path(eph_file), .use_sun = true, .use_moon = true, .name = "third_body"});
    stack.add(std::move(third_body));
    third_body_enabled = true;
  } else {
    spdlog::warn("no ephemeris path provided; third-body curve will be zero");
  }

  out << "altitude_km,drag_mps2,third_body_mps2,total_mps2,status\n";
  for (int i = 0; i < samples; ++i) {
    const double u = static_cast<double>(i) / static_cast<double>(samples - 1);
    const double alt_km = alt_min_km + (alt_max_km - alt_min_km) * u;
    const double r_m = astroforces::atmo::constants::kEarthRadiusWgs84M + alt_km * 1000.0;

    astroforces::atmo::StateVector state{};
    state.epoch.utc_seconds = epoch_utc_s;
    state.frame = astroforces::atmo::Frame::ECI;
    state.position_m = astroforces::atmo::Vec3{r_m, 0.0, 0.0};
    state.velocity_mps = astroforces::atmo::Vec3{0.0, circular_speed_mps(r_m), 0.0};

    const auto result = stack.evaluate(astroforces::forces::PerturbationRequest{.state = state, .spacecraft = &sc});

    double drag_mag = 0.0;
    double third_body_mag = 0.0;
    for (const auto& c : result.contributions) {
      if (c.type == astroforces::forces::PerturbationType::Drag) {
        drag_mag += magnitude(c.acceleration_mps2);
      } else if (c.type == astroforces::forces::PerturbationType::ThirdBody) {
        third_body_mag += magnitude(c.acceleration_mps2);
      }
    }
    if (!third_body_enabled) {
      third_body_mag = 0.0;
    }

    out << fmt::format("{:.6f},{:.12e},{:.12e},{:.12e},{}\n",
                       alt_km,
                       drag_mag,
                       third_body_mag,
                       magnitude(result.total_acceleration_mps2),
                       static_cast<int>(result.status));
  }

  spdlog::info("wrote perturbation profile: {}", out_csv.string());
  return 0;
}

