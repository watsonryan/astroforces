/**
 * @file earth_radiation_perturbation.hpp
 * @brief Earth Radiation perturbation wrapper for generic force stack.
 * @author Watosn
 */
#pragma once

#include <string>

#include "astroforces/forces/surface/earth_radiation/earth_radiation_model.hpp"
#include "astroforces/forces/core/perturbation.hpp"

namespace astroforces::forces {

/**
 * @brief Perturbation-stack wrapper around EarthRadiationAccelerationModel.
 */
class EarthRadiationPerturbationModel final : public astroforces::forces::IPerturbationModel {
 public:
  /**
   * @brief Construct wrapper with Earth radiation model.
   */
  EarthRadiationPerturbationModel(EarthRadiationAccelerationModel earth_radiation = EarthRadiationAccelerationModel{},
                       const astroforces::sc::SpacecraftProperties* default_spacecraft = nullptr,
                       std::string name = "earth_radiation")
      : earth_radiation_(earth_radiation), default_spacecraft_(default_spacecraft), name_(std::move(name)) {}

  /**
   * @brief Evaluate Earth radiation contribution.
   */
  [[nodiscard]] astroforces::forces::PerturbationContribution evaluate(
      const astroforces::forces::PerturbationRequest& request) const override;

 private:
  EarthRadiationAccelerationModel earth_radiation_{};
  const astroforces::sc::SpacecraftProperties* default_spacecraft_{nullptr};
  std::string name_{};
};

}  // namespace astroforces::forces
