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

class GravitySphPerturbationModel final : public IPerturbationModel {
 public:
  explicit GravitySphPerturbationModel(std::unique_ptr<GravitySphAccelerationModel> gravity, std::string name = "gravity_sph")
      : gravity_(std::move(gravity)), name_(std::move(name)) {}

  [[nodiscard]] PerturbationContribution evaluate(const PerturbationRequest& request) const override;

 private:
  std::unique_ptr<GravitySphAccelerationModel> gravity_{};
  std::string name_{};
};

}  // namespace astroforces::forces
