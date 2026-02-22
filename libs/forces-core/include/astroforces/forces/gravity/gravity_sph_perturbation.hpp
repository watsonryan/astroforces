/**
 * @file gravity_sph_perturbation.hpp
 * @brief Gravity perturbation wrapper for the force stack.
 * @author Watosn
 */
#pragma once

#include <memory>
#include <string>

#include "astroforces/forces/gravity/gravity_sph_model.hpp"
#include "astroforces/forces/core/perturbation.hpp"

namespace astroforces::forces {

/**
 * @brief Perturbation-stack wrapper around GravitySphAccelerationModel.
 */
class GravitySphPerturbationModel final : public IPerturbationModel {
 public:
  /**
   * @brief Construct wrapper with owned gravity model.
   */
  explicit GravitySphPerturbationModel(std::unique_ptr<GravitySphAccelerationModel> gravity, std::string name = "gravity_sph")
      : gravity_(std::move(gravity)), name_(std::move(name)) {}

  /**
   * @brief Evaluate gravity contribution.
   */
  [[nodiscard]] PerturbationContribution evaluate(const PerturbationRequest& request) const override;

 private:
  std::unique_ptr<GravitySphAccelerationModel> gravity_{};
  std::string name_{};
};

}  // namespace astroforces::forces
