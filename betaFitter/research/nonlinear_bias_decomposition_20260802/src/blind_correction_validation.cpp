#define main bias_decomposition_standalone_main
#include "bias_decomposition.cpp"
#undef main

#include <numeric>
#include <tuple>

namespace blind_validation {

constexpr int kLearnBins = 20;
constexpr int kIntervals = 3;
constexpr int kSegments = 11;
constexpr double kLearnBetaMin = 0.25;
constexpr double kLearnBetaMax = 0.95;

struct ModelEvent {
  Event event;
  std::array<int, 4> bar{};
  std::array<int, 4> tofPrimary{};
  std::array<int, 9> trackerPrimary{};
  bool primaryStable{false};
};

int learnBin(double beta) {
  if (!(beta >= kLearnBetaMin && beta < kLearnBetaMax))
    return -1;
  int bin = static_cast<int>((beta - kLearnBetaMin) /
                             (kLearnBetaMax - kLearnBetaMin) * kLearnBins);
  return std::max(0, std::min(kLearnBins - 1, bin));
}

double learnCenter(int bin) {
  return kLearnBetaMin +
         (bin + 0.5) * (kLearnBetaMax - kLearnBetaMin) / kLearnBins;
}

double median(const std::vector<double> &values) {
  return quantile(values, 0.5);
}

struct IntervalGeometry {
  double topOuter{};
  double upper{};
  double lower{};
  double bottomOuter{};
};

bool geometry(const Event &event, IntervalGeometry &output) {
  output.topOuter = std::abs(event.path[0] - event.path[1]);
  output.upper = std::abs(event.path[1]);
  output.lower = std::abs(event.path[2]);
  output.bottomOuter = std::abs(event.path[3] - event.path[2]);
  return output.topOuter > 0 && output.topOuter < 30 && output.upper > 0 &&
         output.upper < 150 && output.lower > 0 && output.lower < 150 &&
         output.bottomOuter > 0 && output.bottomOuter < 30;
}

struct SegmentGeometry {
  std::array<double, kSegments> length{};
  std::array<double, 7> upperZ{};
  std::array<double, 6> lowerZ{};
};

bool segmentGeometry(const Event &event, SegmentGeometry &output) {
  const double denominator = event.tkPos[4][2] - event.tkPos[5][2];
  if (!(denominator > 0) || !(event.tkPos[4][2] > 0) ||
      !(event.tkPos[5][2] < 0))
    return false;
  const double fraction = event.tkPos[4][2] / denominator;
  Point center{};
  for (int coordinate = 0; coordinate < 3; ++coordinate)
    center[coordinate] =
        event.tkPos[4][coordinate] +
        fraction * (event.tkPos[5][coordinate] - event.tkPos[4][coordinate]);
  const std::array<Point, 7> upper{ center,         event.tkPos[4],
                                    event.tkPos[3], event.tkPos[2],
                                    event.tkPos[1], event.tofPos[1],
                                    event.tofPos[0] };
  const std::array<Point, 6> lower{ center,          event.tkPos[5],
                                    event.tkPos[6],  event.tkPos[7],
                                    event.tofPos[2], event.tofPos[3] };
  double upperGeometry = 0;
  for (int segment = 0; segment < 5; ++segment)
    upperGeometry += distance(upper[segment], upper[segment + 1]);
  double lowerGeometry = 0;
  for (int segment = 0; segment < 4; ++segment)
    lowerGeometry += distance(lower[segment], lower[segment + 1]);
  if (!(upperGeometry > 0 && lowerGeometry > 0))
    return false;
  const double upperScale = std::abs(event.path[1]) / upperGeometry;
  const double lowerScale = std::abs(event.path[2]) / lowerGeometry;
  for (int segment = 0; segment < 5; ++segment)
    output.length[segment] =
        distance(upper[segment], upper[segment + 1]) * upperScale;
  output.length[5] = std::abs(event.path[0] - event.path[1]);
  for (int segment = 0; segment < 4; ++segment)
    output.length[6 + segment] =
        distance(lower[segment], lower[segment + 1]) * lowerScale;
  output.length[10] = std::abs(event.path[3] - event.path[2]);
  for (int index = 0; index < 7; ++index)
    output.upperZ[index] = upper[index][2];
  for (int index = 0; index < 6; ++index)
    output.lowerZ[index] = lower[index][2];
  return std::all_of(output.length.begin(), output.length.end(),
                     [](double value) {
    return value > 0 && value < 150 && std::isfinite(value);
  });
}

class MaterialModel {
public:
  void add(const Event &event) {
    const int bin = learnBin(event.betaCenter);
    if (bin < 0)
      return;
    SegmentGeometry g;
    if (!segmentGeometry(event, g))
      return;
    std::array<double, 7> upperEnergy{};
    std::array<double, 6> lowerEnergy{};
    for (int index = 0; index < 7; ++index)
      upperEnergy[index] = truthEnergyAt(event, g.upperZ[index]);
    for (int index = 0; index < 6; ++index)
      lowerEnergy[index] = truthEnergyAt(event, g.lowerZ[index]);
    for (int segment = 0; segment < 6; ++segment) {
      const double loss = upperEnergy[segment + 1] - upperEnergy[segment];
      if (std::isfinite(loss))
        samples_[segment][bin].add(loss / g.length[segment]);
    }
    for (int segment = 0; segment < 5; ++segment) {
      const double loss = lowerEnergy[segment] - lowerEnergy[segment + 1];
      if (std::isfinite(loss))
        samples_[6 + segment][bin].add(loss / g.length[6 + segment]);
    }
  }

  void finalize() {
    for (int segment = 0; segment < kSegments; ++segment) {
      for (int bin = 0; bin < kLearnBins; ++bin) {
        counts_[segment][bin] = samples_[segment][bin].v.size();
        if (counts_[segment][bin] >= 20)
          rate_[segment][bin] = std::max(0.0, median(samples_[segment][bin].v));
      }
      fillMissing(rate_[segment]);
      smooth(rate_[segment]);
    }
  }

  bool times(const Event &event, double invBeta,
             std::array<double, 4> &times) const {
    if (!(invBeta > 1.0 && invBeta < 5.0))
      return false;
    SegmentGeometry g;
    if (!segmentGeometry(event, g))
      return false;
    const double beta = 1.0 / invBeta;
    const double mass = event.fitMass;
    const double center = mass / std::sqrt(1.0 - beta * beta);
    if (!(center > mass))
      return false;

    auto inverseBeta = [&](double energy)->double {
      if (!(energy > mass))
        return std::numeric_limits<double>::quiet_NaN();
      const double argument = 1.0 - mass * mass / (energy * energy);
      return argument > 0 ? 1.0 / std::sqrt(argument)
                          : std::numeric_limits<double>::quiet_NaN();
    };
    auto integrate = [&](double length, double first, double second)->double {
      const double a = inverseBeta(first);
      const double m = inverseBeta(0.5 * (first + second));
      const double b = inverseBeta(second);
      if (!std::isfinite(a + m + b))
        return std::numeric_limits<double>::quiet_NaN();
      return length * (a + 4.0 * m + b) / (6.0 * kC);
    };

    times.fill(std::numeric_limits<double>::quiet_NaN());
    double energy = center;
    double accumulated = 0;
    for (int segment = 0; segment < 5; ++segment) {
      const double next = energy + lookup(segment, beta) * g.length[segment];
      const double dt = integrate(g.length[segment], energy, next);
      if (!std::isfinite(dt))
        return false;
      accumulated += dt;
      energy = next;
    }
    times[1] = -accumulated;
    {
      const double next = energy + lookup(5, beta) * g.length[5];
      const double dt = integrate(g.length[5], energy, next);
      if (!std::isfinite(dt))
        return false;
      times[0] = times[1] - dt;
    }

    energy = center;
    accumulated = 0;
    for (int segment = 6; segment < 10; ++segment) {
      const double next = energy - lookup(segment, beta) * g.length[segment];
      const double dt = integrate(g.length[segment], energy, next);
      if (!std::isfinite(dt))
        return false;
      accumulated += dt;
      energy = next;
    }
    times[2] = accumulated;
    {
      const double next = energy - lookup(10, beta) * g.length[10];
      const double dt = integrate(g.length[10], energy, next);
      if (!std::isfinite(dt))
        return false;
      times[3] = times[2] + dt;
    }
    return std::all_of(times.begin(), times.end(),
                       [](double value) { return std::isfinite(value); });
  }

  void write(const std::string &path) const {
    std::ofstream output(path);
    output << "segment,beta_center,n,median_dE_per_cm_GeV\n";
    output << std::setprecision(10);
    for (int segment = 0; segment < kSegments; ++segment)
      for (int bin = 0; bin < kLearnBins; ++bin)
        output << segment << ',' << learnCenter(bin) << ','
               << counts_[segment][bin] << ',' << rate_[segment][bin] << '\n';
  }

private:
  std::array<std::array<Stat, kLearnBins>, kSegments> samples_{};
  std::array<std::array<long long, kLearnBins>, kSegments> counts_{};
  std::array<std::array<double, kLearnBins>, kSegments> rate_{};

  static void fillMissing(std::array<double, kLearnBins> &table) {
    for (int bin = 0; bin < kLearnBins; ++bin) {
      if (table[bin] > 0)
        continue;
      for (int delta = 1; delta < kLearnBins; ++delta) {
        if (bin - delta >= 0 && table[bin - delta] > 0) {
          table[bin] = table[bin - delta];
          break;
        }
        if (bin + delta < kLearnBins && table[bin + delta] > 0) {
          table[bin] = table[bin + delta];
          break;
        }
      }
    }
  }

  static void smooth(std::array<double, kLearnBins> &table) {
    const auto original = table;
    for (int bin = 0; bin < kLearnBins; ++bin) {
      std::vector<double> nearby;
      for (int offset = -1; offset <= 1; ++offset) {
        const int index = bin + offset;
        if (index >= 0 && index < kLearnBins && original[index] > 0)
          nearby.push_back(original[index]);
      }
      if (!nearby.empty())
        table[bin] = median(nearby);
    }
  }

  double lookup(int segment, double beta) const {
    const double position =
        (std::clamp(beta, kLearnBetaMin,
                    std::nextafter(kLearnBetaMax, kLearnBetaMin)) -
         kLearnBetaMin) /
            (kLearnBetaMax - kLearnBetaMin) * kLearnBins -
        0.5;
    const int lower = std::max(
        0, std::min(kLearnBins - 1, static_cast<int>(std::floor(position))));
    const int upper = std::min(kLearnBins - 1, lower + 1);
    const double fraction = std::clamp(position - lower, 0.0, 1.0);
    return rate_[segment][lower] * (1.0 - fraction) +
           rate_[segment][upper] * fraction;
  }
};

#include "constrained_material_model.inc"
#include "kernel_material_model.inc"

struct PairKey {
  int interval{};
  int first{};
  int second{};
  bool operator<(const PairKey &other) const {
    return std::tie(interval, first, second) <
           std::tie(other.interval, other.first, other.second);
  }
};

class TimingModel {
public:
  void train(const std::vector<ModelEvent> &events) {
    for (const auto &item : events) {
      const auto &event = item.event;
      for (int interval = 0; interval < kIntervals; ++interval) {
        const double residual =
            (event.timeMeasured[interval + 1] - event.timeMeasured[interval]) -
            (event.timeTruth[interval + 1] - event.timeTruth[interval]);
        if (event.betaSeed >= 0.85 && event.betaSeed < 0.95) {
          const PairKey key{ interval, item.bar[interval],
                             item.bar[interval + 1] };
          pairSamples_[key].add(residual);
          globalSamples_[interval].add(residual);
        }
      }
    }
    for (int interval = 0; interval < kIntervals; ++interval)
      globalBaseline_[interval] = median(globalSamples_[interval].v);
    for (const auto &entry : pairSamples_) {
      if (entry.second.v.size() >= 10) {
        pairBaseline_[entry.first] = median(entry.second.v);
        pairCounts_[entry.first] = entry.second.v.size();
      }
    }

    for (const auto &item : events) {
      const auto &event = item.event;
      const int bin = learnBin(event.betaSeed);
      if (bin < 0)
        continue;
      for (int interval = 0; interval < kIntervals; ++interval) {
        const double residual =
            (event.timeMeasured[interval + 1] - event.timeMeasured[interval]) -
            (event.timeTruth[interval + 1] - event.timeTruth[interval]);
        responseSamples_[interval][bin]
            .add(residual - baseline(item, interval));
      }
    }
    for (int interval = 0; interval < kIntervals; ++interval) {
      for (int bin = 0; bin < kLearnBins; ++bin) {
        responseCounts_[interval][bin] =
            responseSamples_[interval][bin].v.size();
        if (responseCounts_[interval][bin] >= 20)
          response_[interval][bin] = median(responseSamples_[interval][bin].v);
      }
      fillMissing(response_[interval]);
      smooth(response_[interval]);
    }
  }

  std::array<double, 4> correct(const ModelEvent &item) const {
    std::array<double, 4> output{};
    output[0] = item.event.timeMeasured[0];
    for (int interval = 0; interval < kIntervals; ++interval) {
      const double observed = item.event.timeMeasured[interval + 1] -
                              item.event.timeMeasured[interval];
      output[interval + 1] =
          output[interval] + observed - predict(item, interval);
    }
    return output;
  }

  double predict(const ModelEvent &item, int interval) const {
    return baseline(item, interval) +
           interpolate(response_[interval], item.event.betaSeed);
  }

  void write(const std::string &prefix) const {
    {
      std::ofstream output(prefix + "_pair_baselines.csv");
      output << "interval,bar_first,bar_second,n,baseline_ns\n";
      output << std::setprecision(10);
      for (const auto &entry : pairBaseline_)
        output << entry.first.interval << ',' << entry.first.first << ','
               << entry.first.second << ',' << pairCounts_.at(entry.first)
               << ',' << entry.second << '\n';
    }
    {
      std::ofstream output(prefix + "_beta_response.csv");
      output << "interval,beta_center,n,response_ns\n";
      output << std::setprecision(10);
      for (int interval = 0; interval < kIntervals; ++interval)
        for (int bin = 0; bin < kLearnBins; ++bin)
          output << interval << ',' << learnCenter(bin) << ','
                 << responseCounts_[interval][bin] << ','
                 << response_[interval][bin] << '\n';
    }
  }

private:
  std::map<PairKey, Stat> pairSamples_;
  std::map<PairKey, double> pairBaseline_;
  std::map<PairKey, long long> pairCounts_;
  std::array<Stat, kIntervals> globalSamples_{};
  std::array<double, kIntervals> globalBaseline_{};
  std::array<std::array<Stat, kLearnBins>, kIntervals> responseSamples_{};
  std::array<std::array<long long, kLearnBins>, kIntervals> responseCounts_{};
  std::array<std::array<double, kLearnBins>, kIntervals> response_{};

  double baseline(const ModelEvent &item, int interval) const {
    const PairKey key{ interval, item.bar[interval], item.bar[interval + 1] };
    const auto found = pairBaseline_.find(key);
    return found == pairBaseline_.end() ? globalBaseline_[interval]
                                        : found->second;
  }

  static void fillMissing(std::array<double, kLearnBins> &table) {
    for (int bin = 0; bin < kLearnBins; ++bin) {
      if (std::isfinite(table[bin]) && table[bin] != 0)
        continue;
      for (int delta = 1; delta < kLearnBins; ++delta) {
        if (bin - delta >= 0 && std::isfinite(table[bin - delta]) &&
            table[bin - delta] != 0) {
          table[bin] = table[bin - delta];
          break;
        }
        if (bin + delta < kLearnBins && std::isfinite(table[bin + delta]) &&
            table[bin + delta] != 0) {
          table[bin] = table[bin + delta];
          break;
        }
      }
    }
  }

  static void smooth(std::array<double, kLearnBins> &table) {
    const auto original = table;
    for (int bin = 0; bin < kLearnBins; ++bin) {
      std::vector<double> nearby;
      for (int offset = -1; offset <= 1; ++offset) {
        const int index = bin + offset;
        if (index >= 0 && index < kLearnBins && std::isfinite(original[index]))
          nearby.push_back(original[index]);
      }
      if (!nearby.empty())
        table[bin] = median(nearby);
    }
  }

  static double interpolate(const std::array<double, kLearnBins> &table,
                            double beta) {
    const double position =
        (std::clamp(beta, kLearnBetaMin,
                    std::nextafter(kLearnBetaMax, kLearnBetaMin)) -
         kLearnBetaMin) /
            (kLearnBetaMax - kLearnBetaMin) * kLearnBins -
        0.5;
    const int lower = std::max(
        0, std::min(kLearnBins - 1, static_cast<int>(std::floor(position))));
    const int upper = std::min(kLearnBins - 1, lower + 1);
    const double fraction = std::clamp(position - lower, 0.0, 1.0);
    return table[lower] * (1.0 - fraction) + table[upper] * fraction;
  }
};

std::vector<ModelEvent> loadEvents(const std::vector<std::string> &files,
                                   int begin, int end, int charge) {
  std::vector<ModelEvent> output;
  for (int fileIndex = begin; fileIndex < end; ++fileIndex) {
    TFile file(files[fileIndex].c_str(), "READ");
    auto *tree = dynamic_cast<TTree *>(file.Get("amstreea"));
    if (!tree)
      continue;
    int mpar = 0;
    float mch = 0;
    float momentum[21]{}, coordinate[21][3]{}, time[4]{}, edep[4]{};
    float path[4]{}, tofPosition[4][3]{}, trackerEdep[9]{};
    float trackerPosition[9][3]{}, betaSeed = 0;
    int bar[4]{};
    int tofPrimary[4]{};
    int trackerPrimary[9]{};
    tree->SetBranchStatus("*", 0);
    for (
             const char *branch :
             { "mpar",    "mch",      "mevmom1",   "mevcoo1",
               "tof_tl",  "tof_edep", "tof_leng",  "tof_pos",
               "tk_edep", "tk_pos",   "tof_betah", "tof_barid",
               "mtofpri", "mtrpri" })
      tree->SetBranchStatus(branch, 1);
    tree->SetBranchAddress("mpar", &mpar);
    tree->SetBranchAddress("mch", &mch);
    tree->SetBranchAddress("mevmom1", momentum);
    tree->SetBranchAddress("mevcoo1", coordinate);
    tree->SetBranchAddress("tof_tl", time);
    tree->SetBranchAddress("tof_edep", edep);
    tree->SetBranchAddress("tof_leng", path);
    tree->SetBranchAddress("tof_pos", tofPosition);
    tree->SetBranchAddress("tk_edep", trackerEdep);
    tree->SetBranchAddress("tk_pos", trackerPosition);
    tree->SetBranchAddress("tof_betah", &betaSeed);
    tree->SetBranchAddress("tof_barid", bar);
    tree->SetBranchAddress("mtofpri", tofPrimary);
    tree->SetBranchAddress("mtrpri", trackerPrimary);

    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
      tree->GetEntry(entry);
      if (momentum[10] == -1000)
        continue;
      bool valid = betaSeed > 0 && std::isfinite(betaSeed);
      for (int station = 0; station < 4; ++station)
        valid = valid && std::isfinite(time[station]) && time[station] != -1 &&
                edep[station] > 0 && std::isfinite(path[station]);
      for (int station = 1; station < 4; ++station)
        valid = valid && path[station] - path[station - 1] > 0 &&
                path[station] - path[station - 1] < 300;
      for (int layer = 1; layer <= 7; ++layer)
        valid = valid && trackerEdep[layer] > 0 &&
                std::isfinite(trackerPosition[layer][2]);
      valid = valid && trackerPosition[4][2] > 0 && trackerPosition[5][2] < 0;
      if (!valid)
        continue;

      ModelEvent item;
      Event &event = item.event;
      event.mcMass = mcMassFor(mpar, mch);
      event.fitMass = 2.0 * charge * 0.9314941;
      event.betaSeed = betaSeed;
      for (int index = 0; index < 21; ++index)
        if (momentum[index] > 0 && std::isfinite(momentum[index]) &&
            std::isfinite(coordinate[index][2]))
          event.momentumPoints.push_back(
              MomentumPoint{ coordinate[index][2], momentum[index] });
      normalize(event.momentumPoints);
      const double centerMomentum = momentumAt(event.momentumPoints, 0.0);
      if (!(centerMomentum > 0))
        continue;
      event.centerTruthEnergy = std::sqrt(centerMomentum * centerMomentum +
                                          event.mcMass * event.mcMass);
      event.betaCenter = centerMomentum / event.centerTruthEnergy;
      if (learnBin(event.betaCenter) < 0)
        continue;
      for (int station = 0; station < 4; ++station) {
        event.timeMeasured[station] = time[station];
        event.tofEdepGeV[station] = edep[station] * 1e-3;
        event.path[station] = path[station];
        item.bar[station] = bar[station];
        item.tofPrimary[station] = tofPrimary[station];
        for (int coordinateIndex = 0; coordinateIndex < 3; ++coordinateIndex)
          event.tofPos[station][coordinateIndex] =
              tofPosition[station][coordinateIndex];
      }
      item.primaryStable = true;
      for (int layer = 0; layer < 9; ++layer) {
        event.tkEdepGeV[layer] = trackerEdep[layer] * 1e-3;
        item.trackerPrimary[layer] = trackerPrimary[layer];
        for (int coordinateIndex = 0; coordinateIndex < 3; ++coordinateIndex)
          event.tkPos[layer][coordinateIndex] =
              trackerPosition[layer][coordinateIndex];
      }
      for (int station = 0; station < 4; ++station)
        item.primaryStable = item.primaryStable && item.tofPrimary[station] == 1;
      for (int layer = 1; layer <= 7; ++layer)
        item.primaryStable =
            item.primaryStable && item.trackerPrimary[layer] == 1;
      IntervalGeometry g;
      if (!truthTimes(event) || !geometry(event, g))
        continue;
      output.push_back(std::move(item));
    }
    std::cerr << "loaded file " << fileIndex + 1 << '/' << end
              << ", selected=" << output.size() << '\n';
  }
  return output;
}

void
addResidual(std::map<std::string, std::array<Stat, kLearnBins> > &statistics,
            const std::string &name, const Event &event, const FitResult &fit) {
  const int bin = learnBin(event.betaCenter);
  if (bin >= 0 && fit.ok)
    statistics[name][bin].add(fit.invBeta - 1.0 / event.betaCenter);
}

void writeStatistics(
    const std::string &path,
    std::map<std::string, std::array<Stat, kLearnBins> > &statistics,
    bool timing) {
  std::ofstream output(path);
  output << (timing ? "method,beta_center,n,mean_ns,median_ns,sigma68_ns,rms_"
                      "ns,q05_ns,q95_ns\n"
                    : "method,beta_center,n,mean,median,sigma68,rms,q05,q95\n");
  output << std::setprecision(10);
  for (auto &entry : statistics) {
    for (int bin = 0; bin < kLearnBins; ++bin) {
      const auto result = summarize(entry.second[bin]);
      if (!result.n)
        continue;
      output << entry.first << ',' << learnCenter(bin) << ',' << result.n << ','
             << result.mean << ',' << result.median << ',' << result.sigma68
             << ',' << result.rms << ',' << result.q05 << ',' << result.q95
             << '\n';
    }
  }
}

} // namespace blind_validation

#ifndef BLIND_VALIDATION_NO_MAIN
int main(int argc, char **argv) {
  using namespace blind_validation;
  if (argc < 4) {
    std::cerr << "Usage: " << argv[0]
              << " input.list Z output_prefix [train_files=50] "
                 "[validation_files=50] [max_validation_events=0]\n";
    return 1;
  }
  const std::string listPath = argv[1];
  const int charge = std::stoi(argv[2]);
  const std::string prefix = argv[3];
  const int trainFiles = argc > 4 ? std::stoi(argv[4]) : 50;
  const int validationFiles = argc > 5 ? std::stoi(argv[5]) : 50;
  const int maxValidationEvents = argc > 6 ? std::stoi(argv[6]) : 0;

  const auto files = readList(listPath, 0);
  if (trainFiles <= 0 || validationFiles <= 0 ||
      trainFiles + validationFiles > static_cast<int>(files.size()))
    throw std::runtime_error("Invalid train/validation split");

  std::cerr << "Loading training events...\n";
  auto training = loadEvents(files, 0, trainFiles, charge);
  KernelMaterialModel material;
  for (const auto &item : training)
    material.add(item.event);
  material.finalize();
  TimingModel timing;
  timing.train(training);
  material.write(prefix + "_material_model.csv");
  timing.write(prefix + "_timing_model");

  std::cerr << "Loading validation events...\n";
  auto validation =
      loadEvents(files, trainFiles, trainFiles + validationFiles, charge);
  if (maxValidationEvents > 0 &&
      static_cast<int>(validation.size()) > maxValidationEvents)
    validation.resize(maxValidationEvents);

  std::map<std::string, std::array<Stat, kLearnBins> > residualStats;
  std::map<std::string, std::array<Stat, kLearnBins> > timingStats;
  long long processed = 0;
  for (const auto &item : validation) {
    const Event &event = item.event;
    const auto correctedTimes = timing.correct(item);
    const double seed = 1.0 / std::clamp(event.betaSeed, 0.251, 0.999999);

    Model current = [&](double invBeta, std::array<double, 4> &times) {
      return currentTimes(event, invBeta, 1.9, 1.0, false, times);
    };
    Model learned = [&](double invBeta, std::array<double, 4> &times) {
      return material.times(event, invBeta, times);
    };

    const auto linearMeasured = fitLinear(event, event.timeMeasured);
    const auto linearCorrected = fitLinear(event, correctedTimes);
    const auto currentMeasured = fitModel(event.timeMeasured, current, seed);
    const auto currentCorrected = fitModel(correctedTimes, current, seed);
    const auto materialMeasured = fitModel(event.timeMeasured, learned, seed);
    const auto materialCorrected = fitModel(correctedTimes, learned, seed);
    const auto currentTruth =
        fitModel(event.timeTruth, current, 1.0 / event.betaCenter);
    const auto materialTruth =
        fitModel(event.timeTruth, learned, 1.0 / event.betaCenter);

    addResidual(residualStats, "linear_measured", event, linearMeasured);
    addResidual(residualStats, "linear_timing_corrected", event,
                linearCorrected);
    addResidual(residualStats, "current_measured", event, currentMeasured);
    addResidual(residualStats, "current_timing_corrected", event,
                currentCorrected);
    addResidual(residualStats, "material_measured", event, materialMeasured);
    addResidual(residualStats, "material_timing_corrected", event,
                materialCorrected);
    addResidual(residualStats, "current_truth", event, currentTruth);
    addResidual(residualStats, "material_truth", event, materialTruth);

    std::array<double, 4> currentAtTruth{};
    std::array<double, 4> materialAtTruth{};
    const bool currentValid = currentTimes(event, 1.0 / event.betaCenter, 1.9,
                                           1.0, false, currentAtTruth);
    const bool materialValid =
        material.times(event, 1.0 / event.betaCenter, materialAtTruth);
    const int bin = learnBin(event.betaCenter);
    for (int interval = 0; interval < kIntervals; ++interval) {
      const std::string suffix = "_I" + std::to_string(interval);
      timingStats["measured" + suffix][bin].add(
          (event.timeMeasured[interval + 1] - event.timeMeasured[interval]) -
          (event.timeTruth[interval + 1] - event.timeTruth[interval]));
      timingStats["corrected" + suffix][bin]
          .add((correctedTimes[interval + 1] - correctedTimes[interval]) -
               (event.timeTruth[interval + 1] - event.timeTruth[interval]));
      if (currentValid)
        timingStats["current_model" + suffix][bin]
            .add((currentAtTruth[interval + 1] - currentAtTruth[interval]) -
                 (event.timeTruth[interval + 1] - event.timeTruth[interval]));
      if (materialValid)
        timingStats["material_model" + suffix][bin]
            .add((materialAtTruth[interval + 1] - materialAtTruth[interval]) -
                 (event.timeTruth[interval + 1] - event.timeTruth[interval]));
    }
    ++processed;
    if (processed % 5000 == 0)
      std::cerr << "validated " << processed << '/' << validation.size()
                << '\n';
  }

  writeStatistics(prefix + "_residuals.csv", residualStats, false);
  writeStatistics(prefix + "_timing_residuals.csv", timingStats, true);
  {
    std::ofstream output(prefix + "_summary.txt");
    output << "Z=" << charge << '\n';
    output << "train_files=" << trainFiles << '\n';
    output << "validation_files=" << validationFiles << '\n';
    output << "training_events=" << training.size() << '\n';
    output << "validation_events=" << validation.size() << '\n';
    output << "processed_events=" << processed << '\n';
  }
  std::cerr << "Done. Processed validation events: " << processed << '\n';
  return 0;
}
#endif
#endif
