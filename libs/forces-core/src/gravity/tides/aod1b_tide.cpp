/**
 * @file aod1b_tide.cpp
 * @brief AOD1B gravity de-aliasing SH correction model.
 * @author Watosn
 */

#include "astroforces/forces/gravity/tides/aod1b_tide.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

namespace astroforces::forces::tides {
namespace {

double julian_date_utc_from_ymdhms(int year, int month, int day, int hour, int minute, int second) {
  int y = year;
  int m = month;
  if (m <= 2) {
    y -= 1;
    m += 12;
  }
  const int A = y / 100;
  const int B = 2 - A + A / 4;
  const double day_frac =
      (static_cast<double>(hour) + static_cast<double>(minute) / 60.0 + static_cast<double>(second) / 3600.0) / 24.0;
  return std::floor(365.25 * (y + 4716)) + std::floor(30.6001 * (m + 1)) + static_cast<double>(day) + day_frac + B -
         1524.5;
}

}  // namespace

std::unique_ptr<Aod1bTideModel> Aod1bTideModel::load_from_file(const std::filesystem::path& path, int max_degree) {
  std::ifstream in(path);
  if (!in) {
    return {};
  }

  auto out = std::unique_ptr<Aod1bTideModel>(new Aod1bTideModel(std::max(0, max_degree)));
  std::string line;
  while (std::getline(in, line)) {
    if (line.find("DATA SET") == std::string::npos) {
      continue;
    }

    int data_set_number = 0;
    int total_coefficients = 0;
    int epoch[6] = {0, 0, 0, 0, 0, 0};
    char type[4] = {0, 0, 0, 0};
    const int found = std::sscanf(line.c_str(),
                                  "DATA SET %d: %d COEFFICIENTS FOR %d-%d-%d %d:%d:%d OF TYPE %3s",
                                  &data_set_number,
                                  &total_coefficients,
                                  &epoch[0],
                                  &epoch[1],
                                  &epoch[2],
                                  &epoch[3],
                                  &epoch[4],
                                  &epoch[5],
                                  type);
    if (found != 9 || std::string(type) != "glo") {
      continue;
    }

    (void)data_set_number;
    Snapshot snapshot{};
    snapshot.jd_utc =
        julian_date_utc_from_ymdhms(epoch[0], epoch[1], epoch[2], epoch[3], epoch[4], epoch[5]);
    snapshot.Cnm = Eigen::MatrixXd::Zero(out->max_degree_ + 1, out->max_degree_ + 1);
    snapshot.Snm = Eigen::MatrixXd::Zero(out->max_degree_ + 1, out->max_degree_ + 1);

    for (int i = 0; i < total_coefficients; ++i) {
      if (!std::getline(in, line)) {
        break;
      }
      int n = 0;
      int m = 0;
      double C = 0.0;
      double S = 0.0;
      const int got = std::sscanf(line.c_str(), "%d %d %lf %lf", &n, &m, &C, &S);
      if (got != 4) {
        continue;
      }
      if (n < 0 || m < 0 || m > n || n > out->max_degree_) {
        continue;
      }
      snapshot.Cnm(n, m) = C;
      snapshot.Snm(n, m) = S;
    }
    out->snapshots_.push_back(std::move(snapshot));
  }

  std::sort(out->snapshots_.begin(), out->snapshots_.end(), [](const Snapshot& a, const Snapshot& b) {
    return a.jd_utc < b.jd_utc;
  });
  return out;
}

bool Aod1bTideModel::empty() const noexcept { return snapshots_.empty(); }

bool Aod1bTideModel::interpolate_delta_coefficients(const double jd_utc, Eigen::MatrixXd& dC, Eigen::MatrixXd& dS) const {
  if (empty() || dC.rows() == 0 || dS.rows() == 0) {
    return false;
  }

  const Snapshot* s0 = nullptr;
  const Snapshot* s1 = nullptr;

  if (snapshots_.size() == 1) {
    s0 = &snapshots_.front();
    s1 = &snapshots_.front();
  } else if (jd_utc <= snapshots_.front().jd_utc) {
    s0 = &snapshots_.front();
    s1 = &snapshots_[1];
  } else if (jd_utc >= snapshots_.back().jd_utc) {
    s0 = &snapshots_[snapshots_.size() - 2];
    s1 = &snapshots_.back();
  } else {
    auto it = std::upper_bound(
        snapshots_.begin(), snapshots_.end(), jd_utc, [](double jd, const Snapshot& s) { return jd < s.jd_utc; });
    const std::size_t idx1 = static_cast<std::size_t>(std::distance(snapshots_.begin(), it));
    s0 = &snapshots_[idx1 - 1];
    s1 = &snapshots_[idx1];
  }

  const double dt = s1->jd_utc - s0->jd_utc;
  const double w1 = (std::abs(dt) > 0.0) ? (jd_utc - s0->jd_utc) / dt : 0.0;
  const double w0 = 1.0 - w1;

  const int nr = std::min({dC.rows(), dS.rows(), s0->Cnm.rows(), s1->Cnm.rows()});
  const int nc = std::min({dC.cols(), dS.cols(), s0->Cnm.cols(), s1->Cnm.cols()});
  dC.topLeftCorner(nr, nc) += w0 * s0->Cnm.topLeftCorner(nr, nc) + w1 * s1->Cnm.topLeftCorner(nr, nc);
  dS.topLeftCorner(nr, nc) += w0 * s0->Snm.topLeftCorner(nr, nc) + w1 * s1->Snm.topLeftCorner(nr, nc);
  return true;
}

}  // namespace astroforces::forces::tides
