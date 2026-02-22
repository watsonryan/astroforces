/**
 * @file celestrak_csv_provider.cpp
 * @brief CelesTrak SW-Last5Years CSV provider implementation.
 * @author Watosn
 */

#include "astroforces/weather/celestrak_csv_provider.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>

namespace astroforces::weather {
namespace {

constexpr std::size_t kMinCelesTrakColumns = 31;
constexpr std::size_t kDateCol = 0;
constexpr std::size_t kNdCol = 2;
constexpr std::size_t kKp1Col = 3;
constexpr std::size_t kKpSumCol = 11;
constexpr std::size_t kAp1Col = 12;
constexpr std::size_t kApAvgCol = 20;
constexpr std::size_t kF107ObsCol = 24;
constexpr std::size_t kF107DataTypeCol = 26;
constexpr std::size_t kF107ObsCenter81Col = 27;
constexpr double kSecondsPerDay = 86400.0;
constexpr double kSecondsPer3h = 10800.0;

int days_from_civil(int y, unsigned m, unsigned d) {
  y -= static_cast<int>(m <= 2);
  const int era = (y >= 0 ? y : y - 399) / 400;
  const unsigned yoe = static_cast<unsigned>(y - era * 400);
  const unsigned doy = (153U * (m + (m > 2 ? -3U : 9U)) + 2U) / 5U + d - 1U;
  const unsigned doe = yoe * 365U + yoe / 4U - yoe / 100U + doy;
  return era * 146097 + static_cast<int>(doe) - 719468;
}

double ymd_to_utc_seconds(int y, unsigned m, unsigned d) {
  return static_cast<double>(days_from_civil(y, m, d)) * 86400.0;
}

std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> fields;
  fields.reserve(40);
  std::string token;
  std::stringstream ss(line);
  while (std::getline(ss, token, ',')) {
    fields.push_back(token);
  }
  return fields;
}

bool parse_date(const std::string& text, int& year, unsigned& month, unsigned& day) {
  if (text.size() != 10 || text[4] != '-' || text[7] != '-') {
    return false;
  }
  try {
    year = std::stoi(text.substr(0, 4));
    month = static_cast<unsigned>(std::stoul(text.substr(5, 2)));
    day = static_cast<unsigned>(std::stoul(text.substr(8, 2)));
  } catch (...) {
    return false;
  }
  return month >= 1U && month <= 12U && day >= 1U && day <= 31U;
}

bool parse_double(const std::string& text, double& value) {
  try {
    value = std::stod(text);
  } catch (...) {
    return false;
  }
  return true;
}

bool parse_int(const std::string& text, int& value) {
  try {
    value = std::stoi(text);
  } catch (...) {
    return false;
  }
  return true;
}

bool is_observed_type(const std::string& type_token) {
  return type_token == "OBS";
}

std::pair<std::size_t, int> resolve_day_slot(const std::vector<CelesTrakCsvSpaceWeatherProvider::DailySample>& samples, double t) {
  if (samples.empty()) {
    return {0U, 0};
  }

  if (t <= samples.front().day_start_utc_s) {
    return {0U, 0};
  }
  const auto it = std::upper_bound(samples.begin(), samples.end(), t,
                                   [](double ts, const CelesTrakCsvSpaceWeatherProvider::DailySample& s) {
                                     return ts < s.day_start_utc_s;
                                   });
  std::size_t idx = static_cast<std::size_t>(std::distance(samples.begin(), it));
  if (idx == 0U) {
    idx = 0U;
  } else {
    idx -= 1U;
  }
  if (idx >= samples.size()) {
    idx = samples.size() - 1U;
  }

  const double dt = t - samples[idx].day_start_utc_s;
  int slot = static_cast<int>(dt / kSecondsPer3h);
  if (slot < 0) {
    slot = 0;
  }
  if (slot > 7) {
    slot = 7;
  }
  return {idx, slot};
}

double ap_at_time(const std::vector<CelesTrakCsvSpaceWeatherProvider::DailySample>& samples, double t) {
  const auto idx_slot = resolve_day_slot(samples, t);
  return samples[idx_slot.first].ap_3h_utc[static_cast<std::size_t>(idx_slot.second)];
}

double avg_ap_range_hours_back(const std::vector<CelesTrakCsvSpaceWeatherProvider::DailySample>& samples, double t, int first_hour,
                               int last_hour) {
  double sum = 0.0;
  int n = 0;
  for (int hour = first_hour; hour <= last_hour; hour += 3) {
    sum += ap_at_time(samples, t - static_cast<double>(hour) * 3600.0);
    ++n;
  }
  return (n > 0) ? (sum / static_cast<double>(n)) : 0.0;
}

astroforces::core::WeatherIndices make_indices(const CelesTrakCsvSpaceWeatherProvider::DailySample& s, bool interpolated,
                                           bool extrapolated, int slot, bool has_msis_history,
                                           const std::array<double, 7>& ap_msis_history) {
  astroforces::core::WeatherIndices out{
      .f107 = s.f107_obs,
      .f107a = s.f107_obs_center81,
      .ap = s.ap_avg,
      .kp = s.kp_avg,
      .ap_3h_current = s.ap_3h_utc[static_cast<std::size_t>(slot)],
      .kp_3h_current = s.kp_3h_utc[static_cast<std::size_t>(slot)],
      .ap_3h_utc = s.ap_3h_utc,
      .kp_3h_utc = s.kp_3h_utc,
      .ap_msis_history = ap_msis_history,
      .has_ap_msis_history = has_msis_history,
      .f107_observed = s.f107_observed,
      .geomagnetic_observed = s.geomagnetic_observed,
      .source = astroforces::core::WeatherSource::CelesTrakLast5YearsCsv,
      .interpolated = interpolated,
      .extrapolated = extrapolated,
      .status = astroforces::core::Status::Ok};
  return out;
}

}  // namespace

std::unique_ptr<CelesTrakCsvSpaceWeatherProvider> CelesTrakCsvSpaceWeatherProvider::Create(const Config& config) {
  std::ifstream in(config.csv_file);
  if (!in) {
    return std::unique_ptr<CelesTrakCsvSpaceWeatherProvider>(
        new CelesTrakCsvSpaceWeatherProvider(std::vector<DailySample>{}));
  }

  std::vector<DailySample> samples;
  std::string line;
  bool header_consumed = false;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    if (!header_consumed) {
      header_consumed = true;
      if (line.find("DATE") != std::string::npos) {
        continue;
      }
    }

    const auto fields = split_csv_line(line);
    if (fields.size() < kMinCelesTrakColumns) {
      continue;
    }

    int year = 0;
    unsigned month = 0;
    unsigned day = 0;
    if (!parse_date(fields[kDateCol], year, month, day)) {
      continue;
    }

    double kp_sum = 0.0;
    double ap_avg = 0.0;
    double f107 = 0.0;
    double f107a = 0.0;
    int nd = 0;
    std::array<double, 8> kp_slots{};
    std::array<double, 8> ap_slots{};
    bool ok_slots = true;
    for (std::size_t i = 0; i < 8; ++i) {
      double kp_raw = 0.0;
      double ap_raw = 0.0;
      if (!parse_double(fields[kKp1Col + i], kp_raw) || !parse_double(fields[kAp1Col + i], ap_raw)) {
        ok_slots = false;
        break;
      }
      kp_slots[i] = kp_raw / 10.0;  // CelesTrak KP columns are tenths.
      ap_slots[i] = ap_raw;
    }
    if (!parse_double(fields[kKpSumCol], kp_sum) || !parse_double(fields[kApAvgCol], ap_avg) ||
        !parse_double(fields[kF107ObsCol], f107) || !parse_double(fields[kF107ObsCenter81Col], f107a) ||
        !parse_int(fields[kNdCol], nd) || !ok_slots) {
      continue;
    }
    const bool f107_observed = is_observed_type(fields[kF107DataTypeCol]);
    const bool geomagnetic_observed = nd == 0;

    samples.push_back(DailySample{
        .day_start_utc_s = ymd_to_utc_seconds(year, month, day),
        .f107_obs = f107,
        .f107_obs_center81 = f107a,
        .ap_avg = ap_avg,
        .kp_avg = (kp_sum / 8.0) / 10.0,
        .ap_3h_utc = ap_slots,
        .kp_3h_utc = kp_slots,
        .f107_observed = f107_observed,
        .geomagnetic_observed = geomagnetic_observed,
    });
  }

  std::sort(samples.begin(), samples.end(),
            [](const DailySample& a, const DailySample& b) { return a.day_start_utc_s < b.day_start_utc_s; });
  return std::unique_ptr<CelesTrakCsvSpaceWeatherProvider>(new CelesTrakCsvSpaceWeatherProvider(std::move(samples)));
}

astroforces::core::WeatherIndices CelesTrakCsvSpaceWeatherProvider::at(const astroforces::core::Epoch& epoch) const {
  if (samples_.empty()) {
    return astroforces::core::WeatherIndices{
        .source = astroforces::core::WeatherSource::CelesTrakLast5YearsCsv,
        .status = astroforces::core::Status::DataUnavailable};
  }

  const double t = epoch.utc_seconds;
  const auto idx_slot = resolve_day_slot(samples_, t);
  const std::size_t idx = idx_slot.first;
  const int slot = idx_slot.second;

  std::array<double, 7> ap_msis_history{};
  ap_msis_history[0] = ap_at_time(samples_, t);
  ap_msis_history[1] = ap_at_time(samples_, t - 3.0 * 3600.0);
  ap_msis_history[2] = ap_at_time(samples_, t - 6.0 * 3600.0);
  ap_msis_history[3] = ap_at_time(samples_, t - 9.0 * 3600.0);
  ap_msis_history[4] = avg_ap_range_hours_back(samples_, t, 12, 33);
  ap_msis_history[5] = avg_ap_range_hours_back(samples_, t, 36, 57);
  ap_msis_history[6] = avg_ap_range_hours_back(samples_, t, 60, 81);

  const bool has_msis_history = t >= (samples_.front().day_start_utc_s + 81.0 * 3600.0);

  if (t <= samples_.front().day_start_utc_s) {
    return make_indices(samples_.front(), false, true, 0, false, ap_msis_history);
  }
  if (t >= (samples_.back().day_start_utc_s + kSecondsPerDay)) {
    return make_indices(samples_.back(), false, true, 7, has_msis_history, ap_msis_history);
  }

  const auto it = std::upper_bound(samples_.begin(), samples_.end(), t,
                                   [](double ts, const DailySample& s) { return ts < s.day_start_utc_s; });
  if (it == samples_.begin()) {
    return make_indices(samples_.front(), false, true, slot, false, ap_msis_history);
  }

  const auto& hi = *it;
  const auto& lo = *(it - 1);
  const double dt = hi.day_start_utc_s - lo.day_start_utc_s;
  if (dt <= 0.0) {
    return make_indices(lo, false, false, slot, has_msis_history, ap_msis_history);
  }
  const double alpha = std::clamp((t - lo.day_start_utc_s) / dt, 0.0, 1.0);

  DailySample interp{};
  interp.f107_obs = lo.f107_obs + alpha * (hi.f107_obs - lo.f107_obs);
  interp.f107_obs_center81 = lo.f107_obs_center81 + alpha * (hi.f107_obs_center81 - lo.f107_obs_center81);
  interp.ap_avg = lo.ap_avg + alpha * (hi.ap_avg - lo.ap_avg);
  interp.kp_avg = lo.kp_avg + alpha * (hi.kp_avg - lo.kp_avg);
  interp.ap_3h_utc = samples_[idx].ap_3h_utc;
  interp.kp_3h_utc = samples_[idx].kp_3h_utc;
  interp.f107_observed = lo.f107_observed && hi.f107_observed;
  interp.geomagnetic_observed = samples_[idx].geomagnetic_observed;

  const bool extrapolated = (t < samples_.front().day_start_utc_s) || (t > (samples_.back().day_start_utc_s + kSecondsPerDay));
  return make_indices(interp, alpha > 0.0 && alpha < 1.0, extrapolated, slot, has_msis_history, ap_msis_history);
}

}  // namespace astroforces::weather
