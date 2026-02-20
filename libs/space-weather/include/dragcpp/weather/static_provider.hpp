/**
 * @file static_provider.hpp
 * @brief Constant-value weather provider.
 * @author Watosn
 */
#pragma once

#include "dragcpp/atmo/interfaces.hpp"

namespace dragcpp::weather {

class StaticSpaceWeatherProvider final : public dragcpp::atmo::ISpaceWeatherProvider {
 public:
  explicit StaticSpaceWeatherProvider(dragcpp::atmo::WeatherIndices indices) : indices_(indices) {
    indices_.source = dragcpp::atmo::WeatherSource::StaticProvider;
    indices_.ap_3h_current = indices_.ap;
    indices_.kp_3h_current = indices_.kp;
    indices_.ap_3h_utc.fill(indices_.ap);
    indices_.kp_3h_utc.fill(indices_.kp);
    indices_.ap_msis_history.fill(indices_.ap);
    indices_.has_ap_msis_history = true;
  }
  [[nodiscard]] dragcpp::atmo::WeatherIndices at(const dragcpp::atmo::Epoch& /*epoch*/) const override {
    return indices_;
  }

 private:
  dragcpp::atmo::WeatherIndices indices_{};
};

}  // namespace dragcpp::weather
