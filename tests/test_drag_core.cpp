/**
 * @file test_drag_core.cpp
 * @brief Core drag equation tests.
 * @author Watosn
 */

#include <cmath>
#include <iostream>

#include "dragcpp/drag/drag_model.hpp"
#include "dragcpp/models/exponential_atmosphere.hpp"
#include "dragcpp/sc/spacecraft.hpp"
#include "dragcpp/weather/static_provider.hpp"

namespace {

bool approx(double a, double b, double rel) {
  const double d = std::abs(a - b);
  const double n = std::max(std::abs(b), 1e-30);
  return d / n <= rel;
}

}  // namespace

int main() {
  using namespace dragcpp;
  const atmo::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0, .status = atmo::Status::Ok};
  weather::StaticSpaceWeatherProvider weather(wx);
  models::ExponentialAtmosphereModel atmosphere(1.225, 0.0, 7000.0, 1000.0);
  models::ZeroWindModel wind;
  drag::DragAccelerationModel model(weather, atmosphere, wind);

  atmo::StateVector state{};
  state.frame = atmo::Frame::ECI;
  state.position_m = atmo::Vec3{6378137.0, 0.0, 0.0};
  state.velocity_mps = atmo::Vec3{7500.0, 0.0, 0.0};

  sc::SpacecraftProperties sc{.mass_kg = 1000.0, .reference_area_m2 = 10.0, .cd = 2.0, .use_surface_model = false};

  const auto out = model.evaluate(state, sc);
  if (out.status != atmo::Status::Ok) {
    std::cerr << "status failed\n";
    return 1;
  }

  const double expected_ax = -0.5 * 1.225 * 2.0 * 10.0 / 1000.0 * 7500.0 * 7500.0;
  if (!approx(out.acceleration_mps2.x, expected_ax, 1e-12)) {
    std::cerr << "ax mismatch\n";
    return 2;
  }
  if (!approx(out.acceleration_mps2.y, 0.0, 1e-12) || !approx(out.acceleration_mps2.z, 0.0, 1e-12)) {
    std::cerr << "vector mismatch\n";
    return 3;
  }

  sc::SpacecraftProperties macro_sc{
      .mass_kg = 1000.0,
      .reference_area_m2 = 99.0,
      .cd = 2.0,
      .use_surface_model = true,
      .surfaces = {
          sc::Surface{.normal_body = atmo::Vec3{-1.0, 0.0, 0.0}, .area_m2 = 2.0, .cd = 2.0},
          sc::Surface{.normal_body = atmo::Vec3{0.0, -1.0, 0.0}, .area_m2 = 4.0, .cd = 2.0},
      },
  };

  atmo::StateVector state_identity = state;
  state_identity.body_from_frame_dcm = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  const auto macro_identity = model.evaluate(state_identity, macro_sc);
  if (macro_identity.status != atmo::Status::Ok || !approx(macro_identity.area_m2, 2.0, 1e-12)) {
    std::cerr << "macro identity area mismatch\n";
    return 4;
  }

  atmo::StateVector state_rot = state;
  state_rot.body_from_frame_dcm = {0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};  // +90deg yaw
  const auto macro_rot = model.evaluate(state_rot, macro_sc);
  if (macro_rot.status != atmo::Status::Ok || !approx(macro_rot.area_m2, 4.0, 1e-12)) {
    std::cerr << "macro rotated area mismatch\n";
    return 5;
  }

  if (!approx(out.area_m2, 10.0, 1e-12)) {
    std::cerr << "cannonball area mismatch\n";
    return 6;
  }

  return 0;
}
