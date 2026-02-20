/**
 * @file dtm2020_adapter.cpp
 * @brief DTM2020 adapter implementation.
 * @author Watosn
 */

#include "dragcpp/adapters/dtm2020_adapter.hpp"

#include <utility>

#include "dragcpp/atmo/conversions.hpp"
#include "dtm2020/dtm2020_operational.hpp"

namespace dragcpp::adapters {

class Dtm2020AtmosphereAdapter::Impl {
 public:
  explicit Impl(dtm2020::Dtm2020Operational model) : model_(std::move(model)) {}

  dtm2020::Dtm2020Operational model_;
};

std::unique_ptr<Dtm2020AtmosphereAdapter> Dtm2020AtmosphereAdapter::Create(const Config& config) {
  auto ptr = std::unique_ptr<Dtm2020AtmosphereAdapter>(new Dtm2020AtmosphereAdapter(config));
  auto loaded = dtm2020::Dtm2020Operational::LoadFromFile(config.coeff_file);
  if (!loaded.has_value()) {
    return ptr;
  }
  ptr->impl_ = std::make_shared<Impl>(loaded.value());
  return ptr;
}

dragcpp::atmo::AtmosphereSample Dtm2020AtmosphereAdapter::evaluate(const dragcpp::atmo::StateVector& state,
                                                                    const dragcpp::atmo::WeatherIndices& weather) const {
  if (!impl_) {
    return dragcpp::atmo::AtmosphereSample{.status = dragcpp::atmo::Status::DataUnavailable};
  }
  if (weather.status != dragcpp::atmo::Status::Ok) {
    return dragcpp::atmo::AtmosphereSample{.status = weather.status};
  }

  const auto geo = dragcpp::atmo::spherical_geodetic_from_ecef(state.position_m);
  const auto iyd_sec = dragcpp::atmo::utc_seconds_to_iyd_sec(state.epoch.utc_seconds);
  const int doy = iyd_sec.first % 1000;

  dtm2020::OperationalInputs in{};
  in.altitude_km = geo.alt_m * 1e-3;
  in.latitude_deg = geo.lat_deg;
  in.longitude_deg = geo.lon_deg;
  in.local_time_h = dragcpp::atmo::local_solar_time_hours(state.epoch.utc_seconds, geo.lon_deg);
  in.day_of_year = static_cast<double>(doy);
  in.f107 = weather.f107;
  in.f107m = weather.f107a;
  in.kp_delayed_3h = weather.kp_3h_current;
  in.kp_mean_24h = weather.kp;

  const auto out = impl_->model_.Evaluate(in);
  if (!out.has_value()) {
    return dragcpp::atmo::AtmosphereSample{.status = dragcpp::atmo::Status::NumericalError};
  }

  return dragcpp::atmo::AtmosphereSample{
      .density_kg_m3 = out.value().density_g_cm3 * 1000.0,
      .temperature_k = out.value().temperature_k,
      .status = dragcpp::atmo::Status::Ok};
}

}  // namespace dragcpp::adapters
