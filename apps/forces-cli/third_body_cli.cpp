/**
 * @file third_body_cli.cpp
 * @brief Single-state third-body perturbation CLI.
 * @author Watson
 */

#include <cstdlib>
#include <filesystem>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/gravity/third_body.hpp"

namespace {

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

}  // namespace

int main(int argc, char** argv) {
  double x_eci_m{};
  double y_eci_m{};
  double z_eci_m{};
  double vx_eci_mps{};
  double vy_eci_mps{};
  double vz_eci_mps{};
  double epoch_utc_s{};
  std::filesystem::path eph_file{};
  int use_sun_flag{1};
  int use_moon_flag{1};

  CLI::App app{"Single-state third-body perturbation CLI"};
  app.add_option("x_eci_m", x_eci_m)->required();
  app.add_option("y_eci_m", y_eci_m)->required();
  app.add_option("z_eci_m", z_eci_m)->required();
  app.add_option("vx_eci_mps", vx_eci_mps)->required();
  app.add_option("vy_eci_mps", vy_eci_mps)->required();
  app.add_option("vz_eci_mps", vz_eci_mps)->required();
  app.add_option("epoch_utc_s", epoch_utc_s)->required();
  app.add_option("jpl_ephemeris_file", eph_file)->required();
  app.add_option("use_sun", use_sun_flag)->check(CLI::Range(0, 1))->capture_default_str();
  app.add_option("use_moon", use_moon_flag)->check(CLI::Range(0, 1))->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  const bool use_sun = (use_sun_flag != 0);
  const bool use_moon = (use_moon_flag != 0);
  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{x_eci_m, y_eci_m, z_eci_m};
  state.velocity_mps = astroforces::core::Vec3{vx_eci_mps, vy_eci_mps, vz_eci_mps};
  state.epoch.utc_seconds = epoch_utc_s;
  state.frame = astroforces::core::Frame::ECI;

  if (!std::filesystem::exists(eph_file)) {
    spdlog::error("ephemeris file not found: {}", eph_file.string());
    return 2;
  }

  const auto sun = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_file, .use_sun = use_sun, .use_moon = false, .name = "third_body_sun"});
  const auto moon = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_file, .use_sun = false, .use_moon = use_moon, .name = "third_body_moon"});
  const auto total = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_file, .use_sun = use_sun, .use_moon = use_moon, .name = "third_body_total"});

  const auto req = astroforces::forces::PerturbationRequest{.state = state, .spacecraft = nullptr};
  const auto sun_c = sun->evaluate(req);
  const auto moon_c = moon->evaluate(req);
  const auto total_c = total->evaluate(req);

  if (sun_c.status != astroforces::core::Status::Ok || moon_c.status != astroforces::core::Status::Ok ||
      total_c.status != astroforces::core::Status::Ok) {
    spdlog::error("third-body evaluation failed: sun={}, moon={}, total={}", static_cast<int>(sun_c.status),
                  static_cast<int>(moon_c.status), static_cast<int>(total_c.status));
    return 3;
  }

  fmt::print("sun_ax={} sun_ay={} sun_az={} sun_mag={}\n", sun_c.acceleration_mps2.x, sun_c.acceleration_mps2.y,
             sun_c.acceleration_mps2.z, magnitude(sun_c.acceleration_mps2));
  fmt::print("moon_ax={} moon_ay={} moon_az={} moon_mag={}\n", moon_c.acceleration_mps2.x, moon_c.acceleration_mps2.y,
             moon_c.acceleration_mps2.z, magnitude(moon_c.acceleration_mps2));
  fmt::print("total_ax={} total_ay={} total_az={} total_mag={}\n", total_c.acceleration_mps2.x, total_c.acceleration_mps2.y,
             total_c.acceleration_mps2.z, magnitude(total_c.acceleration_mps2));

  return 0;
}
