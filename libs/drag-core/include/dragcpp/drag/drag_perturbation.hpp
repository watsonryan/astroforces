/**
 * @file drag_perturbation.hpp
 * @brief Drag perturbation model adapter for the generic force interface.
 * @author Watosn
 */
#pragma once

#include <string>

#include "dragcpp/drag/drag_model.hpp"
#include "dragcpp/forces/perturbation.hpp"

namespace dragcpp::drag {

class DragPerturbationModel final : public dragcpp::forces::IPerturbationModel {
 public:
  DragPerturbationModel(const dragcpp::atmo::ISpaceWeatherProvider& weather,
                        const dragcpp::atmo::IAtmosphereModel& atmosphere,
                        const dragcpp::atmo::IWindModel& wind,
                        const dragcpp::sc::SpacecraftProperties* default_spacecraft = nullptr,
                        std::string name = "drag")
      : drag_(weather, atmosphere, wind), default_spacecraft_(default_spacecraft), name_(std::move(name)) {}

  [[nodiscard]] dragcpp::forces::PerturbationContribution evaluate(
      const dragcpp::forces::PerturbationRequest& request) const override;

 private:
  DragAccelerationModel drag_;
  const dragcpp::sc::SpacecraftProperties* default_spacecraft_{nullptr};
  std::string name_{};
};

}  // namespace dragcpp::drag
