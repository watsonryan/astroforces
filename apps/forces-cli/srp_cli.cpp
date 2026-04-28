/**
 * @file srp_cli.cpp
 * @brief Single-state SRP perturbation CLI.
 * @author Watson
 */

#include <cstdlib>
#include <filesystem>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/srp/srp_model.hpp"

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
  double mass_kg{600.0};
  double area_m2{4.0};
  double cr{1.3};
  int use_eclipse_flag{0};

  CLI::App app{"Single-state SRP perturbation CLI"};
  app.add_option("x_eci_m", x_eci_m)->required();
  app.add_option("y_eci_m", y_eci_m)->required();
  app.add_option("z_eci_m", z_eci_m)->required();
  app.add_option("vx_eci_mps", vx_eci_mps)->required();
  app.add_option("vy_eci_mps", vy_eci_mps)->required();
  app.add_option("vz_eci_mps", vz_eci_mps)->required();
  app.add_option("epoch_utc_s", epoch_utc_s)->required();
  app.add_option("jpl_ephemeris_file", eph_file)->required();
  app.add_option("mass_kg", mass_kg)->capture_default_str();
  app.add_option("area_m2", area_m2)->capture_default_str();
  app.add_option("cr", cr)->capture_default_str();
  app.add_option("use_eclipse", use_eclipse_flag)->check(CLI::Range(0, 1))->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  const bool use_eclipse = (use_eclipse_flag != 0);
  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{x_eci_m, y_eci_m, z_eci_m};
  state.velocity_mps = astroforces::core::Vec3{vx_eci_mps, vy_eci_mps, vz_eci_mps};
  state.epoch.utc_seconds = epoch_utc_s;
  state.frame = astroforces::core::Frame::ECI;

  if (!std::filesystem::exists(eph_file)) {
    spdlog::error("ephemeris file not found: {}", eph_file.string());
    return 2;
  }

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = mass_kg, .reference_area_m2 = area_m2, .cd = 2.2, .cr = cr, .use_surface_model = false, .surfaces = {}};

  auto srp = astroforces::forces::SrpAccelerationModel::Create(
      {.ephemeris_file = eph_file, .use_eclipse = use_eclipse});
  const auto out = srp->evaluate(state, sc);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("srp evaluation failed: status={}", static_cast<int>(out.status));
    return 3;
  }

  fmt::print("ax={} ay={} az={} amag={}\n", out.acceleration_mps2.x, out.acceleration_mps2.y, out.acceleration_mps2.z,
             magnitude(out.acceleration_mps2));
  fmt::print("p_pa={} eclipse_factor={} r_sun_m={} area={} cr={} eclipsed={}\n", out.solar_pressure_pa, out.eclipse_factor,
             out.sun_distance_m, out.area_m2, out.cr, out.eclipsed ? 1 : 0);
  return 0;
}
