/**
 * @file test_end_to_end.cpp
 * @brief End-to-end drag pipeline smoke test.
 * @author Watosn
 */

#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/drag/drag_model.hpp"
#include "astroforces/models/exponential_atmosphere.hpp"
#include "astroforces/sc/spacecraft.hpp"
#include "astroforces/weather/static_provider.hpp"

int main() {
  using namespace astroforces;

  const core::WeatherIndices wx{.f107 = 120.0, .f107a = 130.0, .ap = 8.0, .kp = 3.0, .status = core::Status::Ok};
  weather::StaticSpaceWeatherProvider weather(wx);
  models::ExponentialAtmosphereModel atmosphere(3.0e-11, 400e3, 65e3, 900.0);
  models::ZeroWindModel wind;
  drag::DragAccelerationModel model(weather, atmosphere, wind);

  core::StateVector state{};
  state.epoch.utc_seconds = 1.0e9;
  state.frame = core::Frame::ECI;
  state.position_m = core::Vec3{6778137.0, 0.0, 0.0};
  state.velocity_mps = core::Vec3{0.0, 7670.0, 0.0};

  sc::SpacecraftProperties sc{.mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.25, .use_surface_model = false};

  const auto out = model.evaluate(state, sc);
  if (out.status != core::Status::Ok) {
    spdlog::error("pipeline failed");
    return 1;
  }
  if (!(out.density_kg_m3 > 0.0)) {
    spdlog::error("density invalid");
    return 2;
  }
  return 0;
}
