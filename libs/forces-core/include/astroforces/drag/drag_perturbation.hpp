/**
 * @file drag_perturbation.hpp
 * @brief Drag perturbation model adapter for the generic force interface.
 * @author Watosn
 */
#pragma once

#include <string>

#include "astroforces/drag/drag_model.hpp"
#include "astroforces/forces/perturbation.hpp"

namespace astroforces::drag {

class DragPerturbationModel final : public astroforces::forces::IPerturbationModel {
 public:
  DragPerturbationModel(const astroforces::atmo::ISpaceWeatherProvider& weather,
                        const astroforces::atmo::IAtmosphereModel& atmosphere,
                        const astroforces::atmo::IWindModel& wind,
                        const astroforces::sc::SpacecraftProperties* default_spacecraft = nullptr,
                        std::string name = "drag")
      : drag_(weather, atmosphere, wind), default_spacecraft_(default_spacecraft), name_(std::move(name)) {}

  [[nodiscard]] astroforces::forces::PerturbationContribution evaluate(
      const astroforces::forces::PerturbationRequest& request) const override;

 private:
  DragAccelerationModel drag_;
  const astroforces::sc::SpacecraftProperties* default_spacecraft_{nullptr};
  std::string name_{};
};

}  // namespace astroforces::drag
