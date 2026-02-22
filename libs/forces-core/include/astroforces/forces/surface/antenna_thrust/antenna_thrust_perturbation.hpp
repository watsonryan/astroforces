/**
 * @file antenna_thrust_perturbation.hpp
 * @brief Antenna thrust perturbation wrapper for generic force stack.
 * @author Watosn
 */
#pragma once

#include <string>

#include "astroforces/forces/core/perturbation.hpp"
#include "astroforces/forces/surface/antenna_thrust/antenna_thrust_model.hpp"

namespace astroforces::forces {

/**
 * @brief Perturbation-stack wrapper around AntennaThrustAccelerationModel.
 */
class AntennaThrustPerturbationModel final : public astroforces::forces::IPerturbationModel {
 public:
  /**
   * @brief Construct wrapper with antenna thrust model.
   */
  AntennaThrustPerturbationModel(AntennaThrustAccelerationModel antenna_thrust = AntennaThrustAccelerationModel{},
                                 const astroforces::sc::SpacecraftProperties* default_spacecraft = nullptr,
                                 std::string name = "antenna_thrust")
      : antenna_thrust_(antenna_thrust), default_spacecraft_(default_spacecraft), name_(std::move(name)) {}

  /**
   * @brief Evaluate antenna-thrust contribution.
   */
  [[nodiscard]] astroforces::forces::PerturbationContribution evaluate(
      const astroforces::forces::PerturbationRequest& request) const override;

 private:
  AntennaThrustAccelerationModel antenna_thrust_{};
  const astroforces::sc::SpacecraftProperties* default_spacecraft_{nullptr};
  std::string name_{};
};

}  // namespace astroforces::forces
