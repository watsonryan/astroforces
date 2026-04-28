/**
 * @file gravity_sph_cli.cpp
 * @brief Single-state full SPH gravity CLI.
 * @author Watson
 */

#include <cstdlib>
#include <filesystem>
#include <string>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/gravity/gravity_sph_model.hpp"

namespace {

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

}  // namespace

int main(int argc, char** argv) {
  double x_m{};
  double y_m{};
  double z_m{};
  double vx_mps{};
  double vy_mps{};
  double vz_mps{};
  double epoch_utc_s{};
  std::filesystem::path gravity_file{};
  int max_degree{360};
  std::string frame{"eci"};
  std::filesystem::path eph_file{};
  int use_tides_flag{0};

  CLI::App app{"Single-state full SPH gravity CLI"};
  app.add_option("x_m", x_m)->required();
  app.add_option("y_m", y_m)->required();
  app.add_option("z_m", z_m)->required();
  app.add_option("vx_mps", vx_mps)->required();
  app.add_option("vy_mps", vy_mps)->required();
  app.add_option("vz_mps", vz_mps)->required();
  app.add_option("epoch_utc_s", epoch_utc_s)->required();
  app.add_option("gravity_gfc_file", gravity_file)->required();
  app.add_option("max_degree", max_degree)->capture_default_str();
  app.add_option("frame", frame)->check(CLI::IsMember({"eci", "ecef"}))->capture_default_str();
  app.add_option("jpl_ephemeris_file", eph_file)->capture_default_str();
  app.add_option("use_tides", use_tides_flag)->check(CLI::Range(0, 1))->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  const bool use_tides = (use_tides_flag != 0);
  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{x_m, y_m, z_m};
  state.velocity_mps = astroforces::core::Vec3{vx_mps, vy_mps, vz_mps};
  state.epoch.utc_seconds = epoch_utc_s;
  if (!std::filesystem::exists(gravity_file)) {
    spdlog::error("gravity coefficient file not found: {}", gravity_file.string());
    return 2;
  }

  if (frame == "ecef") {
    state.frame = astroforces::core::Frame::ECEF;
  } else if (frame == "eci") {
    state.frame = astroforces::core::Frame::ECI;
  } else {
    spdlog::error("invalid frame: {} (expected eci|ecef)", frame);
    return 3;
  }

  if (use_tides && (eph_file.empty() || !std::filesystem::exists(eph_file))) {
    spdlog::error("tides requested but ephemeris file not found: {}", eph_file.string());
    return 4;
  }

  const auto grav = astroforces::forces::GravitySphAccelerationModel::Create(
      {.gravity_model_file = gravity_file,
       .ephemeris_file = eph_file,
       .max_degree = max_degree,
       .use_central = true,
       .use_sph = true,
       .use_solid_earth_tides = use_tides});

  const auto out = grav->evaluate(state);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("gravity evaluation failed: status={}", static_cast<int>(out.status));
    return 5;
  }

  fmt::print("ax={} ay={} az={} amag={}\n", out.acceleration_mps2.x, out.acceleration_mps2.y, out.acceleration_mps2.z,
             magnitude(out.acceleration_mps2));
  fmt::print("central={} sph={} tide_sun={} tide_moon={}\n", magnitude(out.central_mps2), magnitude(out.sph_mps2),
             magnitude(out.solid_tide_sun_mps2), magnitude(out.solid_tide_moon_mps2));
  return 0;
}
