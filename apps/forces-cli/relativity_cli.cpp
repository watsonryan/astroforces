/**
 * @file relativity_cli.cpp
 * @brief Single-state relativistic perturbation CLI.
 * @author Watson
 */

#include <cstdlib>
#include <filesystem>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/gravity/relativity_model.hpp"

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
  int use_geodesic_flag{1};

  CLI::App app{"Single-state relativistic perturbation CLI"};
  app.add_option("x_eci_m", x_eci_m)->required();
  app.add_option("y_eci_m", y_eci_m)->required();
  app.add_option("z_eci_m", z_eci_m)->required();
  app.add_option("vx_eci_mps", vx_eci_mps)->required();
  app.add_option("vy_eci_mps", vy_eci_mps)->required();
  app.add_option("vz_eci_mps", vz_eci_mps)->required();
  app.add_option("epoch_utc_s", epoch_utc_s)->required();
  app.add_option("jpl_ephemeris_file", eph_file)->required();
  app.add_option("use_geodesic", use_geodesic_flag)->check(CLI::Range(0, 1))->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  const bool use_geodesic = (use_geodesic_flag != 0);
  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{x_eci_m, y_eci_m, z_eci_m};
  state.velocity_mps = astroforces::core::Vec3{vx_eci_mps, vy_eci_mps, vz_eci_mps};
  state.epoch.utc_seconds = epoch_utc_s;
  state.frame = astroforces::core::Frame::ECI;

  if (use_geodesic && !std::filesystem::exists(eph_file)) {
    spdlog::error("ephemeris file not found: {}", eph_file.string());
    return 2;
  }

  const auto rel = astroforces::forces::RelativityAccelerationModel::Create(
      {.ephemeris_file = eph_file, .use_geodesic_precession = use_geodesic});
  const auto out = rel->evaluate(state);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("relativity evaluation failed: status={}", static_cast<int>(out.status));
    return 3;
  }

  fmt::print("ax={} ay={} az={} amag={}\n", out.acceleration_mps2.x, out.acceleration_mps2.y, out.acceleration_mps2.z,
             magnitude(out.acceleration_mps2));
  fmt::print("schw={} geo={} lt={} j2={} rot={}\n", magnitude(out.spherical_central_body_mps2), magnitude(out.geodesic_precession_mps2),
             magnitude(out.lense_thirring_mps2), magnitude(out.oblateness_j2_mps2), magnitude(out.rotational_energy_mps2));
  return 0;
}
