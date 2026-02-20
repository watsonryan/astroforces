/**
 * @file drag_model.hpp
 * @brief Drag acceleration pipeline.
 * @author Watosn
 */
#pragma once

#include "astroforces/atmo/interfaces.hpp"
#include "astroforces/sc/spacecraft.hpp"

namespace astroforces::drag {

struct DragResult {
  astroforces::atmo::Vec3 acceleration_mps2{};
  astroforces::atmo::Vec3 relative_velocity_mps{};
  double density_kg_m3{};
  double temperature_k{};
  double relative_speed_mps{};
  double dynamic_pressure_pa{};
  double area_m2{};
  double cd{};
  astroforces::atmo::WeatherIndices weather{};
  astroforces::atmo::Status status{astroforces::atmo::Status::Ok};
};

class DragAccelerationModel {
 public:
  DragAccelerationModel(const astroforces::atmo::ISpaceWeatherProvider& weather,
                        const astroforces::atmo::IAtmosphereModel& atmosphere,
                        const astroforces::atmo::IWindModel& wind)
      : weather_(weather), atmosphere_(atmosphere), wind_(wind) {}

  [[nodiscard]] DragResult evaluate(const astroforces::atmo::StateVector& state,
                                    const astroforces::sc::SpacecraftProperties& sc) const;

 private:
  const astroforces::atmo::ISpaceWeatherProvider& weather_;
  const astroforces::atmo::IAtmosphereModel& atmosphere_;
  const astroforces::atmo::IWindModel& wind_;
};

}  // namespace astroforces::drag
