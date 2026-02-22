/**
 * @file test_perturbation_stack.cpp
 * @brief Generic perturbation interface tests.
 * @author Watosn
 */

#include <cmath>
#include <memory>
#include <string>

#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/drag/drag_model.hpp"
#include "astroforces/forces/surface/drag/drag_perturbation.hpp"
#include "astroforces/forces/core/perturbation.hpp"
#include "astroforces/models/exponential_atmosphere.hpp"
#include "astroforces/sc/spacecraft.hpp"
#include "astroforces/weather/static_provider.hpp"

namespace {

bool approx(double a, double b, double rel = 1e-12) {
  const double d = std::abs(a - b);
  const double n = std::max(std::abs(b), 1e-30);
  return d / n <= rel;
}

class ConstantPerturbation final : public astroforces::forces::IPerturbationModel {
 public:
  ConstantPerturbation(astroforces::core::Vec3 a, const char* name) : a_(a), name_(name) {}

  [[nodiscard]] astroforces::forces::PerturbationContribution evaluate(
      const astroforces::forces::PerturbationRequest& /*request*/) const override {
    return astroforces::forces::PerturbationContribution{
        .name = name_, .type = astroforces::forces::PerturbationType::Unknown, .acceleration_mps2 = a_, .status = astroforces::core::Status::Ok};
  }

 private:
  astroforces::core::Vec3 a_{};
  std::string name_{};
};

}  // namespace

int main() {
  using namespace astroforces;

  const core::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0, .status = core::Status::Ok};
  weather::StaticSpaceWeatherProvider weather(wx);
  models::ExponentialAtmosphereModel atmosphere(3.0e-11, 400e3, 65e3, 900.0);
  models::ZeroWindModel wind;

  core::StateVector state{};
  state.epoch.utc_seconds = 1.0e9;
  state.frame = core::Frame::ECI;
  state.position_m = core::Vec3{6778137.0, 0.0, 0.0};
  state.velocity_mps = core::Vec3{0.0, 7670.0, 0.0};

  sc::SpacecraftProperties sc{.mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.25, .use_surface_model = false};

  drag::DragAccelerationModel drag_direct(weather, atmosphere, wind);
  const auto direct = drag_direct.evaluate(state, sc);
  if (direct.status != core::Status::Ok) {
    spdlog::error("direct drag failed");
    return 1;
  }

  forces::PerturbationStack stack;
  stack.add(std::make_unique<drag::DragPerturbationModel>(weather, atmosphere, wind, &sc, "drag_main"));
  stack.add(std::make_unique<ConstantPerturbation>(core::Vec3{1.0e-9, -2.0e-9, 3.0e-9}, "bias"));

  const auto result = stack.evaluate(forces::PerturbationRequest{.state = state, .spacecraft = &sc});
  if (result.status != core::Status::Ok || result.contributions.size() != 2U) {
    spdlog::error("stack evaluation failed");
    return 2;
  }
  if (!approx(result.contributions[0].acceleration_mps2.x, direct.acceleration_mps2.x) ||
      !approx(result.contributions[0].acceleration_mps2.y, direct.acceleration_mps2.y) ||
      !approx(result.contributions[0].acceleration_mps2.z, direct.acceleration_mps2.z)) {
    spdlog::error("drag contribution mismatch");
    return 3;
  }

  const core::Vec3 expected_total =
      direct.acceleration_mps2 + core::Vec3{1.0e-9, -2.0e-9, 3.0e-9};
  if (!approx(result.total_acceleration_mps2.x, expected_total.x) ||
      !approx(result.total_acceleration_mps2.y, expected_total.y) ||
      !approx(result.total_acceleration_mps2.z, expected_total.z)) {
    spdlog::error("total acceleration mismatch");
    return 4;
  }

  const drag::DragPerturbationModel drag_no_default(weather, atmosphere, wind, nullptr, "drag_nodefault");
  const auto missing_sc = drag_no_default.evaluate(forces::PerturbationRequest{.state = state, .spacecraft = nullptr});
  if (missing_sc.status != core::Status::InvalidInput) {
    spdlog::error("missing spacecraft should fail");
    return 5;
  }

  return 0;
}
