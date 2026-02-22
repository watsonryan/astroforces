/**
 * @file aod1b_tide.hpp
 * @brief AOD1B gravity de-aliasing SH correction model.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>
#include <vector>

#include <Eigen/Dense>

namespace astroforces::forces::tides {

class Aod1bTideModel final {
 public:
  static std::unique_ptr<Aod1bTideModel> load_from_file(const std::filesystem::path& path, int max_degree);
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] bool interpolate_delta_coefficients(double jd_utc, Eigen::MatrixXd& dC, Eigen::MatrixXd& dS) const;

 private:
  struct Snapshot {
    double jd_utc{0.0};
    Eigen::MatrixXd Cnm{};
    Eigen::MatrixXd Snm{};
  };

  explicit Aod1bTideModel(int max_degree) : max_degree_(max_degree) {}

  int max_degree_{0};
  std::vector<Snapshot> snapshots_{};
};

}  // namespace astroforces::forces::tides
