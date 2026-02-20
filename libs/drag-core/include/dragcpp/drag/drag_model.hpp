/**
 * @file drag_model.hpp
 * @brief Drag acceleration pipeline.
 * @author Watosn
 */
#pragma once

#include "dragcpp/atmo/interfaces.hpp"
#include "dragcpp/sc/spacecraft.hpp"

namespace dragcpp::drag {

struct DragResult {
  dragcpp::atmo::Vec3 acceleration_mps2{};
  dragcpp::atmo::Vec3 relative_velocity_mps{};
  double density_kg_m3{};
  double temperature_k{};
  double relative_speed_mps{};
  double dynamic_pressure_pa{};
  double area_m2{};
  double cd{};
  dragcpp::atmo::WeatherIndices weather{};
  dragcpp::atmo::Status status{dragcpp::atmo::Status::Ok};
};

class DragAccelerationModel {
 public:
  DragAccelerationModel(const dragcpp::atmo::ISpaceWeatherProvider& weather,
                        const dragcpp::atmo::IAtmosphereModel& atmosphere,
                        const dragcpp::atmo::IWindModel& wind)
      : weather_(weather), atmosphere_(atmosphere), wind_(wind) {}

  [[nodiscard]] DragResult evaluate(const dragcpp::atmo::StateVector& state,
                                    const dragcpp::sc::SpacecraftProperties& sc) const;

 private:
  const dragcpp::atmo::ISpaceWeatherProvider& weather_;
  const dragcpp::atmo::IAtmosphereModel& atmosphere_;
  const dragcpp::atmo::IWindModel& wind_;
};

}  // namespace dragcpp::drag
