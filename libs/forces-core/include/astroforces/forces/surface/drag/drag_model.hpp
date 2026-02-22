/**
 * @file drag_model.hpp
 * @brief Drag acceleration pipeline.
 * @author Watosn
 */
#pragma once

#include <memory>

#include "astroforces/core/interfaces.hpp"
#include "astroforces/core/eop.hpp"
#include "astroforces/core/cip.hpp"
#include "astroforces/sc/spacecraft.hpp"

namespace astroforces::forces {

/**
 * @brief Drag frame transform strategy for ECI/ECEF conversions.
 */
enum class DragFrameTransformMode : unsigned char {
  AutoPreferStrict,
  ApproxGmst,
  StrictGcrfItrf,
};

/**
 * @brief Drag model output bundle.
 */
struct DragResult {
  astroforces::core::Vec3 acceleration_mps2{};
  astroforces::core::Vec3 relative_velocity_mps{};
  double density_kg_m3{};
  double temperature_k{};
  double relative_speed_mps{};
  double dynamic_pressure_pa{};
  double area_m2{};
  double cd{};
  astroforces::core::WeatherIndices weather{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

/**
 * @brief Atmospheric drag acceleration model.
 */
class DragAccelerationModel {
 public:
  /**
   * @brief Construct drag model from weather/atmosphere/wind providers.
   * @param weather Space weather provider.
   * @param atmosphere Atmosphere model provider.
   * @param wind Wind model provider.
   * @param transform_mode ECI/ECEF transform strategy.
   * @param eop_series Optional EOP series for strict transforms.
   * @param cip_series Optional CIP series for strict transforms.
   */
  DragAccelerationModel(const astroforces::core::ISpaceWeatherProvider& weather,
                        const astroforces::core::IAtmosphereModel& atmosphere,
                        const astroforces::core::IWindModel& wind,
                        DragFrameTransformMode transform_mode = DragFrameTransformMode::AutoPreferStrict,
                        std::shared_ptr<const astroforces::core::eop::Series> eop_series = nullptr,
                        std::shared_ptr<const astroforces::core::cip::Series> cip_series = nullptr)
      : weather_(weather),
        atmosphere_(atmosphere),
        wind_(wind),
        transform_mode_(transform_mode),
        eop_series_(std::move(eop_series)),
        cip_series_(std::move(cip_series)) {}

  /**
   * @brief Evaluate drag acceleration for a spacecraft state.
   * @param state State vector.
   * @param sc Spacecraft aerodynamic properties.
   * @return Drag output bundle with status and intermediate values.
   */
  [[nodiscard]] DragResult evaluate(const astroforces::core::StateVector& state,
                                    const astroforces::sc::SpacecraftProperties& sc) const;

 private:
  const astroforces::core::ISpaceWeatherProvider& weather_;
  const astroforces::core::IAtmosphereModel& atmosphere_;
  const astroforces::core::IWindModel& wind_;
  DragFrameTransformMode transform_mode_{};
  std::shared_ptr<const astroforces::core::eop::Series> eop_series_{};
  std::shared_ptr<const astroforces::core::cip::Series> cip_series_{};
};

}  // namespace astroforces::forces
