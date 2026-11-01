#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr double kC = 29.9792458;      // cm/ns
constexpr double kSigmaT = 0.1544809;  // ns, current betaFitter convention
constexpr double kInvBetaMin = 1.000001;
constexpr double kInvBetaMax = 3.5;
constexpr std::array<int, 3> kCharges{{2, 6, 8}};

struct RawEvent {
    unsigned int run{};
    unsigned int event{};
    float tofBeta{};
    float time[4]{};
    float edep[4]{};
    float length[4]{};
    float tkQ[2]{};
    float tkQin[2][3]{};
    int tofQs{};
};

struct PreparedEvent {
    unsigned int run{};
    unsigned int event{};
    int fileIndex{-1};
    int Z{};
    double betaSeed{};
    double mass{};
    std::array<double, 4> t{};
    std::array<double, 4> edepMeV{};
    std::array<double, 4> s{};
    std::array<double, 4> seg{};
    double edepUpperPerZ2{};
    double linearFourLayerChi2{std::numeric_limits<double>::infinity()};
};

struct FourLayerFit {
    bool ok{false};
    double invBeta{};
    double sse{};
    std::array<double, 4> residual{};
};

struct PairFit {
    bool ok{false};
    double invBeta{};
    int boundary{0};  // -1: beta=1 boundary, +1: slow/energy-domain boundary
};

struct BranchReader {
    RawEvent r;
    explicit BranchReader(TTree *tree) {
        tree->SetBranchStatus("*", 0);
        for (const char *branch : {"run", "event", "tof_betah", "tof_tl", "tof_edep",
                                   "tof_leng", "tk_q", "tk_qin", "tof_qs"}) {
            tree->SetBranchStatus(branch, 1);
        }
        tree->SetBranchAddress("run", &r.run);
        tree->SetBranchAddress("event", &r.event);
        tree->SetBranchAddress("tof_betah", &r.tofBeta);
        tree->SetBranchAddress("tof_tl", r.time);
        tree->SetBranchAddress("tof_edep", r.edep);
        tree->SetBranchAddress("tof_leng", r.length);
        tree->SetBranchAddress("tk_q", r.tkQ);
        tree->SetBranchAddress("tk_qin", r.tkQin);
        tree->SetBranchAddress("tof_qs", &r.tofQs);
    }
};

bool isChargeOfInterest(int z) {
    return std::find(kCharges.begin(), kCharges.end(), z) != kCharges.end();
}

int reconstructedCharge(const RawEvent &r) {
    const double charge = (r.tkQin[0][2] < 2.5) ? r.tkQ[1] : r.tkQin[0][2];
    return static_cast<int>(charge + 0.5);
}

double massForCharge(int z) {
    return 2.0 * z * 0.9314941;  // match the current production mass assumption
}

bool prepareRange(const RawEvent &r, PreparedEvent &e, double betaMin, double betaMax) {
    if (r.tofQs != 1111 || !(r.tofBeta > betaMin && r.tofBeta < betaMax)) return false;
    const int z = reconstructedCharge(r);
    if (!isChargeOfInterest(z)) return false;

    for (int layer = 0; layer < 4; ++layer) {
        if (!std::isfinite(r.time[layer]) || r.time[layer] == -1 ||
            std::abs(r.time[layer]) > 1000) return false;
        if (!std::isfinite(r.edep[layer]) || r.edep[layer] <= 0) return false;
        if (!std::isfinite(r.length[layer])) return false;
    }

    e.run = r.run;
    e.event = r.event;
    e.Z = z;
    e.betaSeed = r.tofBeta;
    e.mass = massForCharge(z);
    for (int layer = 0; layer < 4; ++layer) {
        e.t[layer] = r.time[layer];
        e.edepMeV[layer] = r.edep[layer];
    }
    e.seg[0] = 0;
    e.s[0] = 0;
    for (int layer = 1; layer < 4; ++layer) {
        e.seg[layer] = static_cast<double>(r.length[layer] - r.length[layer - 1]);
        if (!(e.seg[layer] > 0 && e.seg[layer] < 250)) return false;
        e.s[layer] = e.s[layer - 1] + e.seg[layer];
    }
    if (!(e.s[3] > 100 && e.s[3] < 180)) return false;
    e.edepUpperPerZ2 = (e.edepMeV[0] + e.edepMeV[1]) / (z * z);
    return true;
}

FourLayerFit fitLinearFourLayer(const PreparedEvent &e,
                                double invBetaMin = kInvBetaMin,
                                double invBetaMax = kInvBetaMax) {
    FourLayerFit fit;
    double meanS = 0;
    double meanT = 0;
    for (int layer = 0; layer < 4; ++layer) {
        meanS += e.s[layer];
        meanT += e.t[layer];
    }
    meanS /= 4.0;
    meanT /= 4.0;

    double numerator = 0;
    double denominator = 0;
    for (int layer = 0; layer < 4; ++layer) {
        numerator += (e.s[layer] - meanS) * (e.t[layer] - meanT);
        denominator += (e.s[layer] - meanS) * (e.s[layer] - meanS);
    }
    if (!(denominator > 0)) return fit;
    fit.invBeta = std::clamp(kC * numerator / denominator, invBetaMin, invBetaMax);

    double offset = 0;
    for (int layer = 0; layer < 4; ++layer) {
        offset += e.t[layer] - e.s[layer] * fit.invBeta / kC;
    }
    offset /= 4.0;
    fit.sse = 0;
    for (int layer = 0; layer < 4; ++layer) {
        fit.residual[layer] = e.t[layer] - (offset + e.s[layer] * fit.invBeta / kC);
        fit.sse += fit.residual[layer] * fit.residual[layer];
    }
    fit.ok = std::isfinite(fit.sse);
    return fit;
}

double quantileSorted(const std::vector<float> &values, double probability) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double position = probability * (values.size() - 1);
    const auto index = static_cast<size_t>(position);
    const double fraction = position - index;
    return values[index] * (1.0 - fraction) +
           values[std::min(index + 1, values.size() - 1)] * fraction;
}

struct CalibrationAccumulator {
    std::map<int, std::array<std::vector<float>, 4>> residuals;
    std::map<int, long long> counts;

    void add(const PreparedEvent &event) {
        const auto linear = fitLinearFourLayer(event, 0.6, 1.4);
        if (!linear.ok || linear.sse / (kSigmaT * kSigmaT) > 50) return;
        ++counts[event.Z];
        for (int layer = 0; layer < 4; ++layer) {
            residuals[event.Z][layer].push_back(static_cast<float>(linear.residual[layer]));
        }
    }
};

struct CalibrationResult {
    std::map<int, std::array<double, 4>> residualMedian;
    std::map<int, long long> counts;
};

CalibrationResult finishCalibration(CalibrationAccumulator &accumulator) {
    CalibrationResult result;
    result.counts = accumulator.counts;
    for (int z : kCharges) {
        std::array<double, 4> medians{};
        for (int layer = 0; layer < 4; ++layer) {
            auto &values = accumulator.residuals[z][layer];
            if (values.empty()) continue;
            std::sort(values.begin(), values.end());
            medians[layer] = quantileSorted(values, 0.5);
        }
        result.residualMedian[z] = medians;
    }
    return result;
}

void applyCalibration(PreparedEvent &event, const CalibrationResult &calibration) {
    const auto found = calibration.residualMedian.find(event.Z);
    if (found == calibration.residualMedian.end()) return;
    for (int layer = 0; layer < 4; ++layer) event.t[layer] -= found->second[layer];
}

std::vector<std::string> readList(const std::string &path, int maxFiles) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open input list: " + path);
    std::vector<std::string> files;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') continue;
        files.push_back(line);
        if (maxFiles > 0 && static_cast<int>(files.size()) >= maxFiles) break;
    }
    return files;
}

CalibrationResult computeTrainingCalibration(const std::vector<std::string> &files,
                                              int trainFiles) {
    CalibrationAccumulator accumulator;
    for (int fileIndex = 0; fileIndex < trainFiles; ++fileIndex) {
        TFile input(files[fileIndex].c_str(), "READ");
        auto *tree = dynamic_cast<TTree *>(input.Get("amstreea"));
        if (!tree) {
            std::cerr << "Skip missing training tree: " << files[fileIndex] << '\n';
            continue;
        }
        BranchReader reader(tree);
        for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
            tree->GetEntry(entry);
            PreparedEvent event;
            if (prepareRange(reader.r, event, 0.97, 1.03)) accumulator.add(event);
        }
        if ((fileIndex + 1) % 10 == 0 || fileIndex + 1 == trainFiles) {
            std::cerr << "high-beta training file " << (fileIndex + 1) << '/' << trainFiles << '\n';
        }
    }
    return finishCalibration(accumulator);
}

void writeCalibration(const fs::path &path, const CalibrationResult &training,
                      const CalibrationResult &validationDiagnostic) {
    std::ofstream output(path);
    output << "split,usage,Z,n,residual_S1_ps,residual_S2_ps,residual_S3_ps,residual_S4_ps\n";
    output << std::setprecision(10);
    for (const auto &item : {std::make_pair(std::string("training"), &training),
                             std::make_pair(std::string("validation"), &validationDiagnostic)}) {
        const bool isTraining = item.first == "training";
        for (int z : kCharges) {
            const auto medians = item.second->residualMedian.count(z)
                                     ? item.second->residualMedian.at(z)
                                     : std::array<double, 4>{};
            const long long n = item.second->counts.count(z) ? item.second->counts.at(z) : 0;
            output << item.first << ','
                   << (isTraining ? "applied_to_blind_validation" : "diagnostic_only_not_applied")
                   << ',' << z << ',' << n;
            for (double value : medians) output << ',' << 1000.0 * value;
            output << '\n';
        }
    }
}

double inverseBetaAtShiftedEnergy(double energy, double mass) {
    if (!(energy > mass) || !std::isfinite(energy)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double momentum2 = energy * energy - mass * mass;
    if (!(momentum2 > 0)) return std::numeric_limits<double>::quiet_NaN();
    return energy / std::sqrt(momentum2);
}

double centerCoordinate(const PreparedEvent &event, int layer, double invBeta,
                        double zeta, double edepS2, double edepS3) {
    if (!(invBeta >= kInvBetaMin && invBeta <= kInvBetaMax)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (std::abs(zeta) < 1e-15) return event.s[layer] * invBeta / kC;

    if (layer == 1) return 0.0;
    if (layer == 2) return event.seg[2] * invBeta / kC;

    const double denominator = std::sqrt((invBeta - 1.0) * (invBeta + 1.0));
    if (!(denominator > 0)) return std::numeric_limits<double>::quiet_NaN();
    const double centerEnergy = event.mass * invBeta / denominator;
    if (layer == 0) {
        const double upperInvBeta = inverseBetaAtShiftedEnergy(
            centerEnergy + zeta * edepS2 * 1e-3, event.mass);
        if (!std::isfinite(upperInvBeta)) return std::numeric_limits<double>::quiet_NaN();
        return -event.seg[1] * upperInvBeta / kC;
    }
    if (layer == 3) {
        const double lowerInvBeta = inverseBetaAtShiftedEnergy(
            centerEnergy - zeta * edepS3 * 1e-3, event.mass);
        if (!std::isfinite(lowerInvBeta)) return std::numeric_limits<double>::quiet_NaN();
        return event.seg[2] * invBeta / kC + event.seg[3] * lowerInvBeta / kC;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double pairModelTime(const PreparedEvent &event, int firstLayer, int secondLayer,
                     double invBeta, double zeta, double edepS2, double edepS3) {
    const double first = centerCoordinate(event, firstLayer, invBeta, zeta, edepS2, edepS3);
    const double second = centerCoordinate(event, secondLayer, invBeta, zeta, edepS2, edepS3);
    if (!std::isfinite(first) || !std::isfinite(second)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return second - first;
}

PairFit fitPair(const PreparedEvent &event, int firstLayer, int secondLayer,
                double zeta, double edepS2, double edepS3) {
    PairFit fit;
    const double observed = event.t[secondLayer] - event.t[firstLayer];
    if (!std::isfinite(observed)) return fit;

    if (std::abs(zeta) < 1e-15) {
        const double path = event.s[secondLayer] - event.s[firstLayer];
        if (!(path > 0)) return fit;
        const double unconstrained = kC * observed / path;
        fit.invBeta = std::clamp(unconstrained, kInvBetaMin, kInvBetaMax);
        fit.boundary = unconstrained <= kInvBetaMin ? -1 :
                       unconstrained >= kInvBetaMax ? +1 : 0;
        fit.ok = std::isfinite(fit.invBeta);
        return fit;
    }

    double low = kInvBetaMin;
    double high = kInvBetaMax;
    const double modelLow = pairModelTime(event, firstLayer, secondLayer, low,
                                          zeta, edepS2, edepS3);
    if (!std::isfinite(modelLow)) return fit;

    double modelHigh = pairModelTime(event, firstLayer, secondLayer, high,
                                     zeta, edepS2, edepS3);
    if (!std::isfinite(modelHigh)) {
        double valid = low;
        double invalid = high;
        for (int iteration = 0; iteration < 48; ++iteration) {
            const double middle = 0.5 * (valid + invalid);
            const double value = pairModelTime(event, firstLayer, secondLayer, middle,
                                               zeta, edepS2, edepS3);
            if (std::isfinite(value)) valid = middle;
            else invalid = middle;
        }
        high = valid;
        modelHigh = pairModelTime(event, firstLayer, secondLayer, high,
                                  zeta, edepS2, edepS3);
    }
    if (!std::isfinite(modelHigh) || !(modelHigh >= modelLow)) return fit;

    if (observed <= modelLow) {
        fit.ok = true;
        fit.invBeta = low;
        fit.boundary = -1;
        return fit;
    }
    if (observed >= modelHigh) {
        fit.ok = true;
        fit.invBeta = high;
        fit.boundary = +1;
        return fit;
    }

    for (int iteration = 0; iteration < 40; ++iteration) {
        const double middle = 0.5 * (low + high);
        const double modelMiddle = pairModelTime(event, firstLayer, secondLayer, middle,
                                                 zeta, edepS2, edepS3);
        if (!std::isfinite(modelMiddle)) {
            high = middle;
        } else if (modelMiddle < observed) {
            low = middle;
        } else {
            high = middle;
        }
    }
    fit.ok = true;
    fit.invBeta = 0.5 * (low + high);
    fit.boundary = 0;
    return fit;
}

std::string betaBin(double beta) {
    if (beta < 0.55) return "b0.40_0.55";
    if (beta < 0.65) return "b0.55_0.65";
    if (beta < 0.75) return "b0.65_0.75";
    if (beta < 0.85) return "b0.75_0.85";
    return "b0.85_0.90";
}

std::string edepBin(double edepPerZ2) {
    if (edepPerZ2 < 3.5) return "e0_3.5";
    if (edepPerZ2 < 4.0) return "e3.5_4.0";
    if (edepPerZ2 < 4.5) return "e4.0_4.5";
    if (edepPerZ2 < 5.0) return "e4.5_5.0";
    if (edepPerZ2 < 5.5) return "e5.0_5.5";
    if (edepPerZ2 < 6.5) return "e5.5_6.5";
    return "e6.5_inf";
}

std::vector<std::string> groupKeys(const PreparedEvent &event) {
    const std::string charge = "Z" + std::to_string(event.Z);
    const std::string beta = betaBin(event.betaSeed);
    const std::string edep = edepBin(event.edepUpperPerZ2);
    return {"all", charge, charge + "_" + beta, charge + "_" + edep,
            charge + "_" + beta + "_" + edep};
}

struct RunningStats {
    long long n{0};
    long long betterAbsThanLinear{0};
    long long boundaryA{0};
    long long boundaryB{0};
    long long tail02{0};
    long long tail05{0};
    long double sumDelta{0};
    long double sumDelta2{0};
    long double sumAbsDelta{0};
    long double sumDeltaChange{0};
    long double sumA{0};
    long double sumB{0};
    long double sumA2{0};
    long double sumB2{0};
    long double sumAB{0};
    std::vector<float> deltaValues;

    void add(double estimateA, double estimateB, double linearDelta,
             int pairBoundaryA, int pairBoundaryB, bool isLinear) {
        const double delta = estimateA - estimateB;
        if (!std::isfinite(delta + estimateA + estimateB)) return;
        ++n;
        if (!isLinear && std::abs(delta) < std::abs(linearDelta)) ++betterAbsThanLinear;
        if (pairBoundaryA != 0) ++boundaryA;
        if (pairBoundaryB != 0) ++boundaryB;
        if (std::abs(delta) > 0.2) ++tail02;
        if (std::abs(delta) > 0.5) ++tail05;
        sumDelta += delta;
        sumDelta2 += delta * delta;
        sumAbsDelta += std::abs(delta);
        sumDeltaChange += delta - linearDelta;
        sumA += estimateA;
        sumB += estimateB;
        sumA2 += estimateA * estimateA;
        sumB2 += estimateB * estimateB;
        sumAB += estimateA * estimateB;
        deltaValues.push_back(static_cast<float>(delta));
    }
};

struct Summary {
    long long n{};
    double mean{}, meanErrorNaive{}, rms{}, meanAbs{}, median{}, q16{}, q84{}, sigma68{};
    double medianAbs{}, q95Abs{}, tail02{}, tail05{}, fracBetter{}, meanDeltaChange{};
    double boundaryA{}, boundaryB{};
    double meanA{}, meanB{}, varA{}, varB{}, covariance{}, correlation{};
    double sigmaAFromCov{}, sigmaBFromCov{};
};

Summary summarize(RunningStats &stats, bool isLinear) {
    Summary summary;
    summary.n = stats.n;
    if (!stats.n) return summary;
    const double n = static_cast<double>(stats.n);
    summary.mean = static_cast<double>(stats.sumDelta / n);
    summary.rms = std::sqrt(std::max(0.0,
        static_cast<double>(stats.sumDelta2 / n) - summary.mean * summary.mean));
    summary.meanErrorNaive = summary.rms / std::sqrt(n);
    summary.meanAbs = static_cast<double>(stats.sumAbsDelta / n);
    summary.tail02 = static_cast<double>(stats.tail02) / n;
    summary.tail05 = static_cast<double>(stats.tail05) / n;
    summary.fracBetter = isLinear ? std::numeric_limits<double>::quiet_NaN()
                                  : static_cast<double>(stats.betterAbsThanLinear) / n;
    summary.meanDeltaChange = static_cast<double>(stats.sumDeltaChange / n);
    summary.boundaryA = static_cast<double>(stats.boundaryA) / n;
    summary.boundaryB = static_cast<double>(stats.boundaryB) / n;

    std::sort(stats.deltaValues.begin(), stats.deltaValues.end());
    summary.median = quantileSorted(stats.deltaValues, 0.5);
    summary.q16 = quantileSorted(stats.deltaValues, 0.16);
    summary.q84 = quantileSorted(stats.deltaValues, 0.84);
    summary.sigma68 = 0.5 * (summary.q84 - summary.q16);
    std::vector<float> absolute;
    absolute.reserve(stats.deltaValues.size());
    for (float value : stats.deltaValues) absolute.push_back(std::abs(value));
    std::sort(absolute.begin(), absolute.end());
    summary.medianAbs = quantileSorted(absolute, 0.5);
    summary.q95Abs = quantileSorted(absolute, 0.95);

    summary.meanA = static_cast<double>(stats.sumA / n);
    summary.meanB = static_cast<double>(stats.sumB / n);
    summary.varA = std::max(0.0, static_cast<double>(stats.sumA2 / n) -
                                summary.meanA * summary.meanA);
    summary.varB = std::max(0.0, static_cast<double>(stats.sumB2 / n) -
                                summary.meanB * summary.meanB);
    summary.covariance = static_cast<double>(stats.sumAB / n) - summary.meanA * summary.meanB;
    if (summary.varA > 0 && summary.varB > 0) {
        summary.correlation = summary.covariance / std::sqrt(summary.varA * summary.varB);
    } else {
        summary.correlation = std::numeric_limits<double>::quiet_NaN();
    }
    const double resolutionA2 = summary.varA - summary.covariance;
    const double resolutionB2 = summary.varB - summary.covariance;
    summary.sigmaAFromCov = resolutionA2 >= 0 ? std::sqrt(resolutionA2)
                                              : std::numeric_limits<double>::quiet_NaN();
    summary.sigmaBFromCov = resolutionB2 >= 0 ? std::sqrt(resolutionB2)
                                              : std::numeric_limits<double>::quiet_NaN();
    return summary;
}

struct StatsKey {
    std::string pair;
    std::string quality;
    std::string group;
    std::string model;
    bool operator<(const StatsKey &other) const {
        return std::tie(pair, quality, group, model) <
               std::tie(other.pair, other.quality, other.group, other.model);
    }
};

struct FileStatsKey {
    int fileIndex{};
    std::string pair;
    std::string quality;
    std::string group;
    std::string model;
    bool operator<(const FileStatsKey &other) const {
        return std::tie(fileIndex, pair, quality, group, model) <
               std::tie(other.fileIndex, other.pair, other.quality, other.group, other.model);
    }
};

struct PairDefinition {
    std::string name;
    int aFirst{}, aSecond{}, bFirst{}, bSecond{};
};

struct ModelDefinition {
    std::string name;
};

struct ModelPairResult {
    PairFit a;
    PairFit b;
};

void writeSummaryCSV(const fs::path &path, std::map<StatsKey, RunningStats> &statistics) {
    std::ofstream output(path);
    output << "pair,quality,group,model,n,mean_delta_invbeta,mean_error_naive,rms_delta_invbeta,"
              "mean_abs_delta,median_delta_invbeta,q16,q84,sigma68,median_abs_delta,q95_abs_delta,"
              "tail_abs_gt_0p2,tail_abs_gt_0p5,frac_abs_better_than_linear,mean_change_from_linear,"
              "frac_boundary_A,frac_boundary_B,mean_invbeta_A,mean_invbeta_B,var_A,var_B,cov_AB,"
              "corr_AB,sigma_A_from_cov,sigma_B_from_cov\n";
    output << std::setprecision(10);
    for (auto &[key, stats] : statistics) {
        const bool isLinear = key.model == "linear";
        const auto s = summarize(stats, isLinear);
        output << key.pair << ',' << key.quality << ',' << key.group << ',' << key.model << ','
               << s.n << ',' << s.mean << ',' << s.meanErrorNaive << ',' << s.rms << ','
               << s.meanAbs << ',' << s.median << ',' << s.q16 << ',' << s.q84 << ',' << s.sigma68
               << ',' << s.medianAbs << ',' << s.q95Abs << ',' << s.tail02 << ',' << s.tail05 << ','
               << s.fracBetter << ',' << s.meanDeltaChange << ',' << s.boundaryA << ',' << s.boundaryB
               << ',' << s.meanA << ',' << s.meanB << ',' << s.varA << ',' << s.varB << ','
               << s.covariance << ',' << s.correlation << ',' << s.sigmaAFromCov << ','
               << s.sigmaBFromCov << '\n';
    }
}

void writeFileCSV(const fs::path &path, const std::vector<std::string> &files,
                  std::map<FileStatsKey, RunningStats> &statistics) {
    std::ofstream output(path);
    output << "file_index,file,pair,quality,group,model,n,mean_delta_invbeta,median_delta_invbeta,"
              "sigma68,median_abs_delta,q95_abs_delta,tail_abs_gt_0p2,frac_boundary_A,frac_boundary_B\n";
    output << std::setprecision(10);
    for (auto &[key, stats] : statistics) {
        const auto s = summarize(stats, key.model == "linear");
        output << key.fileIndex << ',' << files[key.fileIndex] << ',' << key.pair << ','
               << key.quality << ',' << key.group << ',' << key.model << ',' << s.n << ','
               << s.mean << ',' << s.median << ',' << s.sigma68 << ',' << s.medianAbs << ','
               << s.q95Abs << ',' << s.tail02 << ',' << s.boundaryA << ',' << s.boundaryB << '\n';
    }
}

}  // namespace

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " input.list output_dir [train_files=65] [max_files=0]"
                  << " [global_zeta_Z2=3.258719208] [global_zeta_Z6=2.838085203]"
                  << " [global_zeta_Z8=2.601233755]\n";
        return 1;
    }

    const std::string listPath = argv[1];
    const fs::path outputDir = argv[2];
    int trainFiles = argc > 3 ? std::stoi(argv[3]) : 65;
    const int maxFiles = argc > 4 ? std::stoi(argv[4]) : 0;
    const std::map<int, double> globalZeta{{2, argc > 5 ? std::stod(argv[5]) : 3.258719208},
                                           {6, argc > 6 ? std::stod(argv[6]) : 2.838085203},
                                           {8, argc > 7 ? std::stod(argv[7]) : 2.601233755}};
    fs::create_directories(outputDir);

    const auto files = readList(listPath, maxFiles);
    if (files.size() < 2) throw std::runtime_error("Need at least two input files");
    trainFiles = std::clamp(trainFiles, 1, static_cast<int>(files.size()) - 1);

    std::cerr << "Computing frozen high-beta calibration from training files only...\n";
    const auto trainingCalibration = computeTrainingCalibration(files, trainFiles);

    std::vector<PreparedEvent> events;
    events.reserve(650000);
    CalibrationAccumulator validationCalibrationAccumulator;
    long long validationEntries = 0;
    std::vector<long long> selectedPerFile(files.size(), 0);
    for (int fileIndex = trainFiles; fileIndex < static_cast<int>(files.size()); ++fileIndex) {
        TFile input(files[fileIndex].c_str(), "READ");
        auto *tree = dynamic_cast<TTree *>(input.Get("amstreea"));
        if (!tree) {
            std::cerr << "Skip missing validation tree: " << files[fileIndex] << '\n';
            continue;
        }
        BranchReader reader(tree);
        validationEntries += tree->GetEntries();
        for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
            tree->GetEntry(entry);

            PreparedEvent highBeta;
            if (prepareRange(reader.r, highBeta, 0.97, 1.03)) {
                validationCalibrationAccumulator.add(highBeta);
            }

            PreparedEvent event;
            if (!prepareRange(reader.r, event, 0.4, 0.9)) continue;
            event.fileIndex = fileIndex;
            applyCalibration(event, trainingCalibration);
            const auto fourLayer = fitLinearFourLayer(event);
            if (!fourLayer.ok) continue;
            event.linearFourLayerChi2 = fourLayer.sse / (kSigmaT * kSigmaT);
            events.push_back(event);
            ++selectedPerFile[fileIndex];
        }
        const int done = fileIndex - trainFiles + 1;
        if (done % 10 == 0 || fileIndex + 1 == static_cast<int>(files.size())) {
            std::cerr << "validation file " << done << '/' << (files.size() - trainFiles)
                      << ", selected=" << events.size() << '\n';
        }
    }
    const auto validationCalibrationDiagnostic =
        finishCalibration(validationCalibrationAccumulator);
    writeCalibration(outputDir / "high_beta_layer_calibration.csv", trainingCalibration,
                     validationCalibrationDiagnostic);

    // Deterministic event-level null: jointly move S2/S3 deposits to another event
    // with the same charge and a 0.01-wide built-in-beta stratum.  This preserves the
    // charge/beta deposit distribution while destroying event-specific association.
    std::map<std::pair<int, int>, std::vector<size_t>> shuffleStrata;
    for (size_t index = 0; index < events.size(); ++index) {
        const int fineBetaBin = std::clamp(
            static_cast<int>(std::floor((events[index].betaSeed - 0.4) / 0.01)), 0, 49);
        shuffleStrata[{events[index].Z, fineBetaBin}].push_back(index);
    }
    std::vector<size_t> donor(events.size());
    std::iota(donor.begin(), donor.end(), 0);
    std::mt19937_64 generator(0x13A24B20260728ULL);
    long long singletonShuffleEvents = 0;
    for (auto &[stratum, indices] : shuffleStrata) {
        (void)stratum;
        if (indices.size() < 2) {
            singletonShuffleEvents += indices.size();
            continue;
        }
        std::shuffle(indices.begin(), indices.end(), generator);
        for (size_t position = 0; position < indices.size(); ++position) {
            donor[indices[position]] = indices[(position + 1) % indices.size()];
        }
    }

    const std::array<PairDefinition, 2> pairs{{
        {"S1S3_vs_S2S4", 0, 2, 1, 3},
        {"S1S4_vs_S2S3", 0, 3, 1, 2},
    }};
    const std::array<ModelDefinition, 6> models{{
        {"linear"},
        {"center_global"},
        {"center_zeta_1p9"},
        {"center_global_shuffled_edep"},
        {"center_global_swapped_S2S3_edep"},
        {"center_global_reverse_loss"},
    }};

    std::map<StatsKey, RunningStats> statistics;
    std::map<FileStatsKey, RunningStats> fileStatistics;
    std::array<long long, 2> commonPairEvents{{0, 0}};
    std::array<long long, 2> failedPairEvents{{0, 0}};

    for (size_t eventIndex = 0; eventIndex < events.size(); ++eventIndex) {
        const auto &event = events[eventIndex];
        const auto &donorEvent = events[donor[eventIndex]];
        const auto groups = groupKeys(event);
        const std::string chargeBetaGroup = "Z" + std::to_string(event.Z) + "_" +
                                            betaBin(event.betaSeed);
        const double zeta = globalZeta.at(event.Z);

        for (size_t pairIndex = 0; pairIndex < pairs.size(); ++pairIndex) {
            const auto &pair = pairs[pairIndex];
            std::array<ModelPairResult, 6> fitResults;
            bool allModelsValid = true;
            for (size_t modelIndex = 0; modelIndex < models.size(); ++modelIndex) {
                double modelZeta = zeta;
                double edepS2 = event.edepMeV[1];
                double edepS3 = event.edepMeV[2];
                if (modelIndex == 0) {
                    modelZeta = 0;
                } else if (modelIndex == 2) {
                    modelZeta = 1.9;
                } else if (modelIndex == 3) {
                    edepS2 = donorEvent.edepMeV[1];
                    edepS3 = donorEvent.edepMeV[2];
                } else if (modelIndex == 4) {
                    std::swap(edepS2, edepS3);
                } else if (modelIndex == 5) {
                    modelZeta = -zeta;
                }
                fitResults[modelIndex].a = fitPair(event, pair.aFirst, pair.aSecond,
                                                   modelZeta, edepS2, edepS3);
                fitResults[modelIndex].b = fitPair(event, pair.bFirst, pair.bSecond,
                                                   modelZeta, edepS2, edepS3);
                if (!fitResults[modelIndex].a.ok || !fitResults[modelIndex].b.ok) {
                    allModelsValid = false;
                }
            }
            if (!allModelsValid) {
                ++failedPairEvents[pairIndex];
                continue;
            }
            ++commonPairEvents[pairIndex];
            const double linearDelta = fitResults[0].a.invBeta - fitResults[0].b.invBeta;
            std::array<std::string, 2> qualities{{"base", "linear_core10"}};
            const int qualityCount = event.linearFourLayerChi2 < 10 ? 2 : 1;
            for (int qualityIndex = 0; qualityIndex < qualityCount; ++qualityIndex) {
                const auto &quality = qualities[qualityIndex];
                for (const auto &group : groups) {
                    for (size_t modelIndex = 0; modelIndex < models.size(); ++modelIndex) {
                        const auto &result = fitResults[modelIndex];
                        statistics[{pair.name, quality, group, models[modelIndex].name}].add(
                            result.a.invBeta, result.b.invBeta, linearDelta,
                            result.a.boundary, result.b.boundary, modelIndex == 0);
                    }
                }

                // File-level stability is kept compact: the physically primary Z-beta
                // bins, both disjoint comparisons, and only the two tested models.
                for (size_t modelIndex = 0; modelIndex < 2; ++modelIndex) {
                    const auto &result = fitResults[modelIndex];
                    fileStatistics[{event.fileIndex, pair.name, quality,
                                    chargeBetaGroup, models[modelIndex].name}].add(
                        result.a.invBeta, result.b.invBeta, linearDelta,
                        result.a.boundary, result.b.boundary, modelIndex == 0);
                }
            }
        }

        if ((eventIndex + 1) % 100000 == 0 || eventIndex + 1 == events.size()) {
            std::cerr << "pair inversion " << (eventIndex + 1) << '/' << events.size() << '\n';
        }
    }

    writeSummaryCSV(outputDir / "pair_closure_summary.csv", statistics);
    writeFileCSV(outputDir / "pair_closure_file_stability.csv", files, fileStatistics);

    {
        std::ofstream output(outputDir / "event_counts.csv");
        output << "file_index,file,selected_low_beta\n";
        for (int fileIndex = trainFiles; fileIndex < static_cast<int>(files.size()); ++fileIndex) {
            output << fileIndex << ',' << files[fileIndex] << ',' << selectedPerFile[fileIndex] << '\n';
        }
    }
    {
        std::ofstream output(outputDir / "run_summary.txt");
        output << std::setprecision(12);
        output << "files_total=" << files.size() << '\n';
        output << "train_files=" << trainFiles << '\n';
        output << "validation_files=" << files.size() - trainFiles << '\n';
        output << "validation_entries_read=" << validationEntries << '\n';
        output << "validation_selected=" << events.size() << '\n';
        output << "calibration_source=training_high_beta_only\n";
        output << "validation_high_beta_calibration=diagnostic_only_not_applied\n";
        output << "shuffle_matching=joint_S2S3_edep_within_Z_and_0.01_tof_betah\n";
        output << "shuffle_seed=0x13A24B20260728\n";
        output << "shuffle_singleton_events=" << singletonShuffleEvents << '\n';
        output << "pair_S1S3_vs_S2S4_common=" << commonPairEvents[0] << '\n';
        output << "pair_S1S3_vs_S2S4_failed=" << failedPairEvents[0] << '\n';
        output << "pair_S1S4_vs_S2S3_common=" << commonPairEvents[1] << '\n';
        output << "pair_S1S4_vs_S2S3_failed=" << failedPairEvents[1] << '\n';
        output << "reference=AMSCenter\n";
        for (int z : kCharges) output << "mc_global_zeta_Z" << z << '=' << globalZeta.at(z) << '\n';
        output << "primary_selection=base_event_common_no_four_layer_SSE_cut\n";
        output << "linear_core10_selection=diagnostic_only_model_asymmetric\n";
    }

    std::cerr << "Done. Results in " << outputDir << '\n';
    return 0;
}
