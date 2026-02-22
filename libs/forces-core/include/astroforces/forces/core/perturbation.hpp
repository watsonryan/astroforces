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

enum class PerturbationType : unsigned char { Unknown, Drag, Gravity, SRP, EarthRadiation, Relativity, ThirdBody, AntennaThrust };

struct PerturbationRequest {
  astroforces::core::StateVector state{};
  const astroforces::sc::SpacecraftProperties* spacecraft{nullptr};
};

struct PerturbationContribution {
  std::string name{};
  PerturbationType type{PerturbationType::Unknown};
  astroforces::core::Vec3 acceleration_mps2{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

struct PerturbationResult {
  astroforces::core::Vec3 total_acceleration_mps2{};
  std::vector<PerturbationContribution> contributions{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
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
