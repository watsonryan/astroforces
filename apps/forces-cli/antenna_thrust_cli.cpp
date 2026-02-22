/**
 * @file antenna_thrust_cli.cpp
 * @brief Single-state antenna thrust perturbation CLI.
 * @author Watosn
 */

#include <cstdlib>
#include <string>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/antenna_thrust/antenna_thrust_model.hpp"

namespace {

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

astroforces::forces::AntennaThrustDirectionMode parse_mode(const std::string& s, bool* ok) {
  if (s == "velocity") {
    *ok = true;
    return astroforces::forces::AntennaThrustDirectionMode::Velocity;
  }
  if (s == "nadir") {
    *ok = true;
    return astroforces::forces::AntennaThrustDirectionMode::Nadir;
  }
  if (s == "custom_eci") {
    *ok = true;
    return astroforces::forces::AntennaThrustDirectionMode::CustomEci;
  }
  if (s == "body_fixed") {
    *ok = true;
    return astroforces::forces::AntennaThrustDirectionMode::BodyFixed;
  }
  *ok = false;
  return astroforces::forces::AntennaThrustDirectionMode::Velocity;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 8 || argc > 15) {
    spdlog::error(
        "usage: antenna_thrust_cli <x_eci_m> <y_eci_m> <z_eci_m> <vx_eci_mps> <vy_eci_mps> <vz_eci_mps> <epoch_utc_s> [mass_kg] [transmit_power_w] [efficiency] [mode:velocity|nadir|custom_eci|body_fixed] [dir_x] [dir_y] [dir_z]");
    return 1;
  }

  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
  state.velocity_mps = astroforces::core::Vec3{std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  state.epoch.utc_seconds = std::atof(argv[7]);
  state.frame = astroforces::core::Frame::ECI;

  const double mass_kg = (argc >= 9) ? std::atof(argv[8]) : 600.0;
  const double power_w = (argc >= 10) ? std::atof(argv[9]) : 20.0;
  const double efficiency = (argc >= 11) ? std::atof(argv[10]) : 1.0;
  const std::string mode_s = (argc >= 12) ? argv[11] : "velocity";

  bool mode_ok = false;
  const auto mode = parse_mode(mode_s, &mode_ok);
  if (!mode_ok) {
    spdlog::error("invalid mode: {}", mode_s);
    return 2;
  }

  const auto dir = astroforces::core::Vec3{
      (argc >= 13) ? std::atof(argv[12]) : 1.0,
      (argc >= 14) ? std::atof(argv[13]) : 0.0,
      (argc >= 15) ? std::atof(argv[14]) : 0.0,
  };

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = mass_kg, .reference_area_m2 = 4.0, .cd = 2.2, .cr = 1.3, .use_surface_model = false, .surfaces = {}};

  const astroforces::forces::AntennaThrustAccelerationModel model({
      .transmit_power_w = power_w,
      .efficiency = efficiency,
      .direction_mode = mode,
      .custom_direction_eci = dir,
      .body_axis = dir,
  });
  const auto out = model.evaluate(state, sc);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("antenna_thrust evaluation failed: status={}", static_cast<int>(out.status));
    return 3;
  }

  fmt::print("ax={} ay={} az={} amag={}\n",
             out.acceleration_mps2.x,
             out.acceleration_mps2.y,
             out.acceleration_mps2.z,
             magnitude(out.acceleration_mps2));
  fmt::print("thrust_n={} effective_power_w={} mass_kg={} dir_x={} dir_y={} dir_z={}\n",
             out.thrust_n,
             out.effective_power_w,
             out.mass_kg,
             out.direction_eci.x,
             out.direction_eci.y,
             out.direction_eci.z);
  return 0;
}
