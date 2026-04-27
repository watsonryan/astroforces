/**
 * @file earth_radiation_cli.cpp
 * @brief Single-state Earth radiation pressure perturbation CLI.
 * @author Watosn
 */

#include <cstdlib>
#include <string>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/earth_radiation/earth_radiation_model.hpp"

namespace {

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

}  // namespace

int main(int argc, char** argv) {
  double x_frame_m{};
  double y_frame_m{};
  double z_frame_m{};
  double vx_frame_mps{};
  double vy_frame_mps{};
  double vz_frame_mps{};
  double epoch_utc_s{};
  double mass_kg{600.0};
  double area_m2{4.0};
  double cr{1.3};
  std::string eph_file{};

  CLI::App app{"Single-state Earth radiation pressure perturbation CLI"};
  app.add_option("x_frame_m", x_frame_m)->required();
  app.add_option("y_frame_m", y_frame_m)->required();
  app.add_option("z_frame_m", z_frame_m)->required();
  app.add_option("vx_frame_mps", vx_frame_mps)->required();
  app.add_option("vy_frame_mps", vy_frame_mps)->required();
  app.add_option("vz_frame_mps", vz_frame_mps)->required();
  app.add_option("epoch_utc_s", epoch_utc_s)->required();
  app.add_option("mass_kg", mass_kg)->capture_default_str();
  app.add_option("area_m2", area_m2)->capture_default_str();
  app.add_option("cr", cr)->capture_default_str();
  app.add_option("jpl_ephemeris_file", eph_file)->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{x_frame_m, y_frame_m, z_frame_m};
  state.velocity_mps = astroforces::core::Vec3{vx_frame_mps, vy_frame_mps, vz_frame_mps};
  state.epoch.utc_seconds = epoch_utc_s;
  state.frame = astroforces::core::Frame::ECI;

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = mass_kg, .reference_area_m2 = area_m2, .cd = 2.2, .cr = cr, .use_surface_model = false, .surfaces = {}};

  auto earth_radiation = astroforces::forces::EarthRadiationAccelerationModel::Create({
      .ephemeris_file = eph_file,
  });
  if (!earth_radiation) {
    spdlog::error("failed to create Earth Radiation model");
    return 2;
  }
  const auto out = earth_radiation->evaluate(state, sc);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("earth_radiation evaluation failed: status={}", static_cast<int>(out.status));
    return 3;
  }

  fmt::print("ax={} ay={} az={} amag={}\n", out.acceleration_mps2.x, out.acceleration_mps2.y, out.acceleration_mps2.z,
             magnitude(out.acceleration_mps2));
  fmt::print("p_total_pa={} p_albedo_pa={} p_ir_pa={} albedo_phase={} albedo_eclipse_factor={} r_earth_m={} area={} cr={}\n",
             out.earth_radiation_pressure_pa, out.albedo_pressure_pa, out.ir_pressure_pa, out.albedo_phase_function,
             out.albedo_eclipse_factor, out.earth_distance_m, out.area_m2, out.cr);
  return 0;
}
