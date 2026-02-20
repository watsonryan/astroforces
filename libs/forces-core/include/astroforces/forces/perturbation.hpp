/**
 * @file perturbation.hpp
 * @brief Generic perturbation model interfaces and aggregation utilities.
 * @author Watosn
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "astroforces/atmo/types.hpp"
#include "astroforces/sc/spacecraft.hpp"

namespace astroforces::forces {

enum class PerturbationType : unsigned char { Unknown, Drag, Gravity, SRP, ThirdBody };

struct PerturbationRequest {
  astroforces::atmo::StateVector state{};
  const astroforces::sc::SpacecraftProperties* spacecraft{nullptr};
};

struct PerturbationContribution {
  std::string name{};
  PerturbationType type{PerturbationType::Unknown};
  astroforces::atmo::Vec3 acceleration_mps2{};
  astroforces::atmo::Status status{astroforces::atmo::Status::Ok};
};

struct PerturbationResult {
  astroforces::atmo::Vec3 total_acceleration_mps2{};
  std::vector<PerturbationContribution> contributions{};
  astroforces::atmo::Status status{astroforces::atmo::Status::Ok};
};

class IPerturbationModel {
 public:
  virtual ~IPerturbationModel() = default;
  [[nodiscard]] virtual PerturbationContribution evaluate(const PerturbationRequest& request) const = 0;
};

class PerturbationStack final {
 public:
  void add(std::unique_ptr<IPerturbationModel> model);
  [[nodiscard]] PerturbationResult evaluate(const PerturbationRequest& request) const;

 private:
  std::vector<std::unique_ptr<IPerturbationModel>> models_{};
};

}  // namespace astroforces::forces
