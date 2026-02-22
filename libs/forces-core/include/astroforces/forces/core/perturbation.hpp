/**
 * @file perturbation.hpp
 * @brief Generic perturbation model interfaces and aggregation utilities.
 * @author Watosn
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "astroforces/core/types.hpp"
#include "astroforces/sc/spacecraft.hpp"

namespace astroforces::forces {

/**
 * @brief High-level perturbation category identifier.
 */
enum class PerturbationType : unsigned char { Unknown, Drag, Gravity, SRP, EarthRadiation, Relativity, ThirdBody, AntennaThrust };

/**
 * @brief Evaluation request passed to a perturbation model.
 */
struct PerturbationRequest {
  astroforces::core::StateVector state{};
  const astroforces::sc::SpacecraftProperties* spacecraft{nullptr};
};

/**
 * @brief One model contribution to total acceleration.
 */
struct PerturbationContribution {
  std::string name{};
  PerturbationType type{PerturbationType::Unknown};
  astroforces::core::Vec3 acceleration_mps2{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

/**
 * @brief Aggregated perturbation output across all active models.
 */
struct PerturbationResult {
  astroforces::core::Vec3 total_acceleration_mps2{};
  std::vector<PerturbationContribution> contributions{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

/**
 * @brief Interface for one force-model contribution.
 */
class IPerturbationModel {
 public:
  virtual ~IPerturbationModel() = default;
  /**
   * @brief Evaluate one perturbation contribution.
   * @param request Input state and spacecraft properties.
   * @return Contribution with acceleration and status.
   */
  [[nodiscard]] virtual PerturbationContribution evaluate(const PerturbationRequest& request) const = 0;
};

/**
 * @brief Composes multiple perturbation models into one aggregate result.
 */
class PerturbationStack final {
 public:
  /**
   * @brief Add one model to the stack.
   */
  void add(std::unique_ptr<IPerturbationModel> model);
  /**
   * @brief Evaluate all added models and return total acceleration.
   */
  [[nodiscard]] PerturbationResult evaluate(const PerturbationRequest& request) const;

 private:
  std::vector<std::unique_ptr<IPerturbationModel>> models_{};
};

}  // namespace astroforces::forces
