/**
 * @file test_drag_core.cpp
 * @brief Core drag equation tests.
 * @author Watosn
 */

#include <cmath>

#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/drag/drag_model.hpp"
#include "astroforces/models/exponential_atmosphere.hpp"
#include "astroforces/sc/spacecraft.hpp"
#include "astroforces/weather/static_provider.hpp"

namespace {

bool approx(double a, double b, double rel) {
  const double d = std::abs(a - b);
  const double n = std::max(std::abs(b), 1e-30);
  return d / n <= rel;
}

}  // namespace

int main() {
  using namespace astroforces;
  const core::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0, .status = core::Status::Ok};
  weather::StaticSpaceWeatherProvider weather(wx);
  models::ExponentialAtmosphereModel atmosphere(1.225, 0.0, 7000.0, 1000.0);
  models::ZeroWindModel wind;
  drag::DragAccelerationModel model(weather, atmosphere, wind);

  core::StateVector state{};
  state.frame = core::Frame::ECI;
  state.position_m = core::Vec3{6378137.0, 0.0, 0.0};
  state.velocity_mps = core::Vec3{7500.0, 0.0, 0.0};

  sc::SpacecraftProperties sc{.mass_kg = 1000.0, .reference_area_m2 = 10.0, .cd = 2.0, .use_surface_model = false};

  const auto out = model.evaluate(state, sc);
  if (out.status != core::Status::Ok) {
    spdlog::error("status failed");
    return 1;
  }

  const double expected_ax = -0.5 * 1.225 * 2.0 * 10.0 / 1000.0 * 7500.0 * 7500.0;
  if (!approx(out.acceleration_mps2.x, expected_ax, 1e-12)) {
    spdlog::error("ax mismatch");
    return 2;
  }
  if (!approx(out.acceleration_mps2.y, 0.0, 1e-12) || !approx(out.acceleration_mps2.z, 0.0, 1e-12)) {
    spdlog::error("vector mismatch");
    return 3;
  }
  if (!approx(out.relative_speed_mps, 7500.0, 1e-12) || !approx(out.dynamic_pressure_pa, 0.5 * 1.225 * 7500.0 * 7500.0, 1e-12)) {
    spdlog::error("derived drag scalars mismatch");
    return 7;
  }

  sc::SpacecraftProperties macro_sc{
      .mass_kg = 1000.0,
      .reference_area_m2 = 99.0,
      .cd = 2.0,
      .use_surface_model = true,
      .surfaces = {
          sc::Surface{.normal_body = core::Vec3{-1.0, 0.0, 0.0}, .area_m2 = 2.0, .cd = 2.5},
          sc::Surface{.normal_body = core::Vec3{0.0, -1.0, 0.0}, .area_m2 = 4.0, .cd = 3.5},
      },
  };

  core::StateVector state_identity = state;
  state_identity.body_from_frame_dcm = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  const auto macro_identity = model.evaluate(state_identity, macro_sc);
  if (macro_identity.status != core::Status::Ok || !approx(macro_identity.area_m2, 2.0, 1e-12) || !approx(macro_identity.cd, 2.5, 1e-12)) {
    spdlog::error("macro identity area mismatch");
    return 4;
  }

  core::StateVector state_rot = state;
  state_rot.body_from_frame_dcm = {0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};  // +90deg yaw
  const auto macro_rot = model.evaluate(state_rot, macro_sc);
  if (macro_rot.status != core::Status::Ok || !approx(macro_rot.area_m2, 4.0, 1e-12) || !approx(macro_rot.cd, 3.5, 1e-12)) {
    spdlog::error("macro rotated area mismatch");
    return 5;
  }

  if (!approx(out.area_m2, 10.0, 1e-12)) {
    spdlog::error("cannonball area mismatch");
    return 6;
  }

  sc::SpacecraftProperties aero_sc{
      .mass_kg = 1000.0,
      .reference_area_m2 = 99.0,
      .cd = 2.0,
      .use_surface_model = true,
      .surfaces = {
          sc::Surface{.normal_body = core::Vec3{-1.0, 0.0, 0.0}, .area_m2 = 2.0, .cd = 2.0, .specularity = 0.0, .accommodation = 1.0},
      },
  };
  const auto aero_out = model.evaluate(state_identity, aero_sc);
  if (aero_out.status != core::Status::Ok || !(aero_out.cd > 2.0)) {
    spdlog::error("surface aero modifier mismatch");
    return 8;
  }

  return 0;
}
