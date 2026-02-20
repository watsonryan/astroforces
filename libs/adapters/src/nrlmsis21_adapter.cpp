/**
 * @file nrlmsis21_adapter.cpp
 * @brief NRLMSIS 2.1 adapter implementation.
 * @author Watosn
 */

#include "dragcpp/adapters/nrlmsis21_adapter.hpp"

#include <utility>

#include "dragcpp/atmo/conversions.hpp"
#include "msis21/msis21.hpp"

namespace dragcpp::adapters {

class Nrlmsis21AtmosphereAdapter::Impl {
 public:
  explicit Impl(msis21::Model model) : model_(std::move(model)) {}

  msis21::Model model_;
};

std::unique_ptr<Nrlmsis21AtmosphereAdapter> Nrlmsis21AtmosphereAdapter::Create(const Config& config) {
  auto ptr = std::unique_ptr<Nrlmsis21AtmosphereAdapter>(new Nrlmsis21AtmosphereAdapter(config));
  msis21::Options options{};
  auto model = msis21::Model::load_from_file(config.parm_file, options);
  ptr->impl_ = std::make_shared<Impl>(std::move(model));
  return ptr;
}

dragcpp::atmo::AtmosphereSample Nrlmsis21AtmosphereAdapter::evaluate(const dragcpp::atmo::StateVector& state,
                                                                      const dragcpp::atmo::WeatherIndices& weather) const {
  if (!impl_) {
    return dragcpp::atmo::AtmosphereSample{.status = dragcpp::atmo::Status::DataUnavailable};
  }
  if (weather.status != dragcpp::atmo::Status::Ok) {
    return dragcpp::atmo::AtmosphereSample{.status = weather.status};
  }

  const auto geo = dragcpp::atmo::spherical_geodetic_from_ecef(state.position_m);
  const auto iyd_sec = dragcpp::atmo::utc_seconds_to_iyd_sec(state.epoch.utc_seconds);

  msis21::Input in{};
  in.iyd = iyd_sec.first;
  in.sec = iyd_sec.second;
  in.alt_km = geo.alt_m * 1e-3;
  in.glat_deg = geo.lat_deg;
  in.glon_deg = geo.lon_deg;
  in.stl_hr = dragcpp::atmo::local_solar_time_hours(state.epoch.utc_seconds, geo.lon_deg);
  in.f107a = weather.f107a;
  in.f107 = weather.f107;
  in.ap = weather.ap_3h_current;
  in.ap_history = weather.ap_msis_history;
  in.has_ap_history = weather.has_ap_msis_history;

  const auto out = impl_->model_.evaluate(in);
  if (out.status != msis21::Status::Ok) {
    return dragcpp::atmo::AtmosphereSample{.status = dragcpp::atmo::Status::NumericalError};
  }

  return dragcpp::atmo::AtmosphereSample{
      .density_kg_m3 = out.out.rho * 1000.0,
      .temperature_k = out.out.t,
      .status = dragcpp::atmo::Status::Ok};
}

}  // namespace dragcpp::adapters
