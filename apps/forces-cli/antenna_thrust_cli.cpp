/**
 * @file antenna_thrust_cli.cpp
 * @brief Single-state antenna thrust perturbation CLI.
 * @author Watson
 */

#include <cstdlib>
#include <string>

#include <CLI/CLI.hpp>
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
  double x_eci_m{};
  double y_eci_m{};
  double z_eci_m{};
  double vx_eci_mps{};
  double vy_eci_mps{};
  double vz_eci_mps{};
  double epoch_utc_s{};
  double mass_kg{600.0};
  double power_w{20.0};
  double efficiency{1.0};
  std::string mode_s{"velocity"};
  double dir_x{1.0};
  double dir_y{0.0};
  double dir_z{0.0};

  CLI::App app{"Single-state antenna thrust perturbation CLI"};
  app.add_option("x_eci_m", x_eci_m)->required();
  app.add_option("y_eci_m", y_eci_m)->required();
  app.add_option("z_eci_m", z_eci_m)->required();
  app.add_option("vx_eci_mps", vx_eci_mps)->required();
  app.add_option("vy_eci_mps", vy_eci_mps)->required();
  app.add_option("vz_eci_mps", vz_eci_mps)->required();
  app.add_option("epoch_utc_s", epoch_utc_s)->required();
  app.add_option("mass_kg", mass_kg)->capture_default_str();
  app.add_option("transmit_power_w", power_w)->capture_default_str();
  app.add_option("efficiency", efficiency)->capture_default_str();
  app.add_option("mode", mode_s)
      ->check(CLI::IsMember({"velocity", "nadir", "custom_eci", "body_fixed"}))
      ->capture_default_str();
  app.add_option("dir_x", dir_x)->capture_default_str();
  app.add_option("dir_y", dir_y)->capture_default_str();
  app.add_option("dir_z", dir_z)->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{x_eci_m, y_eci_m, z_eci_m};
  state.velocity_mps = astroforces::core::Vec3{vx_eci_mps, vy_eci_mps, vz_eci_mps};
  state.epoch.utc_seconds = epoch_utc_s;
  state.frame = astroforces::core::Frame::ECI;

  bool mode_ok = false;
  const auto mode = parse_mode(mode_s, &mode_ok);
  if (!mode_ok) {
    spdlog::error("invalid mode: {}", mode_s);
    return 2;
  }

  const auto dir = astroforces::core::Vec3{dir_x, dir_y, dir_z};

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
