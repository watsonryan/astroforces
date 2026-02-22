/**
 * @file drag_perturbation.hpp
 * @brief Drag perturbation model adapter for the generic force interface.
 * @author Watosn
 */
#pragma once

#include <string>

#include "astroforces/forces/surface/drag/drag_model.hpp"
#include "astroforces/forces/core/perturbation.hpp"

namespace astroforces::forces {

/**
 * @brief Perturbation-stack wrapper around DragAccelerationModel.
 */
class DragPerturbationModel final : public astroforces::forces::IPerturbationModel {
 public:
  /**
   * @brief Construct drag wrapper with provider references.
   */
  DragPerturbationModel(const astroforces::core::ISpaceWeatherProvider& weather,
                        const astroforces::core::IAtmosphereModel& atmosphere,
                        const astroforces::core::IWindModel& wind,
                        const astroforces::sc::SpacecraftProperties* default_spacecraft = nullptr,
                        std::string name = "drag")
      : drag_(weather, atmosphere, wind), default_spacecraft_(default_spacecraft), name_(std::move(name)) {}

  /**
   * @brief Evaluate drag contribution.
   */
  [[nodiscard]] astroforces::forces::PerturbationContribution evaluate(
      const astroforces::forces::PerturbationRequest& request) const override;

 private:
  DragAccelerationModel drag_;
  const astroforces::sc::SpacecraftProperties* default_spacecraft_{nullptr};
  std::string name_{};
};

}  // namespace astroforces::forces
