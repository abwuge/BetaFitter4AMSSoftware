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
#include <set>
#include <sstream>
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
constexpr double kInf = 1e100;
constexpr std::array<int, 3> kCharges{{2, 6, 8}};
bool gCenterReference = false;

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
    int Z{};
    double betaSeed{};
    double mass{};
    std::array<double, 4> t{};
    std::array<double, 4> edepMeV{};
    std::array<double, 4> s{};
    std::array<double, 4> seg{};
    double edepUpperPerZ2{};
};

struct FitResult {
    bool ok{false};
    double invBeta{0};
    double beta{0};
    double offset{0};
    double sse{kInf};
    std::array<double, 4> q{};
    std::array<double, 4> residual{};
};

bool isChargeOfInterest(int z) {
    return std::find(kCharges.begin(), kCharges.end(), z) != kCharges.end();
}

double massForCharge(int z) {
    return 2.0 * z * 0.9314941;  // exactly match current production assumption
}

int reconstructedCharge(const RawEvent &r) {
    const double q = (r.tkQin[0][2] < 2.5) ? r.tkQ[1] : r.tkQin[0][2];
    return static_cast<int>(q + 0.5);
}

bool prepareRange(const RawEvent &r, PreparedEvent &e, double betaMin, double betaMax) {
    if (r.tofQs != 1111 || !(r.tofBeta > betaMin && r.tofBeta < betaMax)) return false;
    const int z = reconstructedCharge(r);
    if (!isChargeOfInterest(z)) return false;

    for (int i = 0; i < 4; ++i) {
        if (!std::isfinite(r.time[i]) || r.time[i] == -1 || std::abs(r.time[i]) > 1000) return false;
        if (!std::isfinite(r.edep[i]) || r.edep[i] <= 0) return false;
        if (!std::isfinite(r.length[i])) return false;
    }

    e.run = r.run;
    e.event = r.event;
    e.Z = z;
    e.betaSeed = r.tofBeta;
    e.mass = massForCharge(z);
    for (int i = 0; i < 4; ++i) {
        e.t[i] = r.time[i];
        e.edepMeV[i] = r.edep[i];
    }
    e.seg[0] = 0;
    e.s[0] = 0;
    for (int i = 1; i < 4; ++i) {
        e.seg[i] = static_cast<double>(r.length[i] - r.length[i - 1]);
        if (!(e.seg[i] > 0 && e.seg[i] < 250)) return false;
        e.s[i] = e.s[i - 1] + e.seg[i];
    }
    if (!(e.s[3] > 100 && e.s[3] < 180)) return false;
    e.edepUpperPerZ2 = (e.edepMeV[0] + e.edepMeV[1]) / (z * z);
    return true;
}

bool prepare(const RawEvent &r, PreparedEvent &e) {
    return prepareRange(r, e, 0.4, 0.9);
}

bool modelTimes(const PreparedEvent &e, double invBeta, double zeta,
                std::array<double, 4> &q) {
    if (!(invBeta >= kInvBetaMin && invBeta <= kInvBetaMax)) return false;
    const double beta = 1.0 / invBeta;
    if (!(beta > 0 && beta < 1)) return false;

    q.fill(0.0);
    if (std::abs(zeta) < 1e-14) {
        for (int i = 1; i < 4; ++i) q[i] = q[i - 1] + e.seg[i] * invBeta / kC;
        return true;
    }

    const double m2 = e.mass * e.mass;
    double energy = e.mass / std::sqrt(1.0 - beta * beta);
    if (gCenterReference) {
        // Match the current BetaNL::propagate(AMSCenter) definition exactly.
        // An arbitrary common shift is omitted because profileSSE profiles it out.
        const double energyUpper = energy + zeta * e.edepMeV[1] * 1e-3;
        const double energyLower = energy - zeta * e.edepMeV[2] * 1e-3;
        if (!(energyUpper > e.mass) || !(energyLower > e.mass)) return false;
        const double upperArg = 1.0 - m2 / (energyUpper * energyUpper);
        const double lowerArg = 1.0 - m2 / (energyLower * energyLower);
        if (!(upperArg > 0 && upperArg < 1) || !(lowerArg > 0 && lowerArg < 1)) return false;
        q[1] = 0.0;
        q[2] = e.seg[2] * invBeta / kC;
        q[0] = -e.seg[1] / (kC * std::sqrt(upperArg));
        q[3] = q[2] + e.seg[3] / (kC * std::sqrt(lowerArg));
        return true;
    }

    for (int i = 1; i < 4; ++i) {
        energy -= zeta * e.edepMeV[i - 1] * 1e-3;
        if (!(energy > e.mass) || !std::isfinite(energy)) return false;
        const double arg = 1.0 - m2 / (energy * energy);
        if (!(arg > 0 && arg < 1)) return false;
        const double nextInvBeta = 1.0 / std::sqrt(arg);
        q[i] = q[i - 1] + e.seg[i] * nextInvBeta / kC;
    }
    return true;
}

double profileSSE(const PreparedEvent &e, double invBeta, double zeta,
                  double *offset = nullptr,
                  std::array<double, 4> *qOut = nullptr,
                  std::array<double, 4> *rOut = nullptr) {
    std::array<double, 4> q{};
    if (!modelTimes(e, invBeta, zeta, q)) return kInf;
    double a = 0;
    for (int i = 0; i < 4; ++i) a += e.t[i] - q[i];
    a /= 4.0;
    double sse = 0;
    std::array<double, 4> residual{};
    for (int i = 0; i < 4; ++i) {
        residual[i] = e.t[i] - (a + q[i]);
        sse += residual[i] * residual[i];
    }
    if (offset) *offset = a;
    if (qOut) *qOut = q;
    if (rOut) *rOut = residual;
    return sse;
}

FitResult fitLinear(const PreparedEvent &e,
                    double invBetaMin = kInvBetaMin,
                    double invBetaMax = kInvBetaMax) {
    FitResult out;
    double meanS = 0, meanT = 0;
    for (int i = 0; i < 4; ++i) { meanS += e.s[i]; meanT += e.t[i]; }
    meanS /= 4.0;
    meanT /= 4.0;
    double num = 0, den = 0;
    for (int i = 0; i < 4; ++i) {
        num += (e.s[i] - meanS) * (e.t[i] - meanT);
        den += (e.s[i] - meanS) * (e.s[i] - meanS);
    }
    if (!(den > 0)) return out;
    double invBeta = kC * num / den;
    invBeta = std::clamp(invBeta, invBetaMin, invBetaMax);
    out.invBeta = invBeta;
    out.beta = 1.0 / invBeta;
    for (int i = 0; i < 4; ++i) out.q[i] = e.s[i] * invBeta / kC;
    out.offset = 0;
    for (int i = 0; i < 4; ++i) out.offset += e.t[i] - out.q[i];
    out.offset /= 4.0;
    out.sse = 0;
    for (int i = 0; i < 4; ++i) {
        out.residual[i] = e.t[i] - (out.offset + out.q[i]);
        out.sse += out.residual[i] * out.residual[i];
    }
    out.ok = std::isfinite(out.sse);
    return out;
}

FitResult fitAtZeta(const PreparedEvent &e, double zeta, const FitResult *linearHint = nullptr) {
    if (std::abs(zeta) < 1e-14) return linearHint ? *linearHint : fitLinear(e);
    FitResult out;

    // Global coarse scan avoids dependence on the production Minuit starting point.
    constexpr int nGrid = 18;
    double bestX = kInvBetaMin;
    double bestF = profileSSE(e, bestX, zeta);
    int bestJ = 0;
    for (int j = 1; j <= nGrid; ++j) {
        const double x = kInvBetaMin + (kInvBetaMax - kInvBetaMin) * j / nGrid;
        const double f = profileSSE(e, x, zeta);
        if (f < bestF) { bestF = f; bestX = x; bestJ = j; }
    }
    if (linearHint && linearHint->ok) {
        const double f = profileSSE(e, linearHint->invBeta, zeta);
        if (f < bestF) { bestF = f; bestX = linearHint->invBeta; }
    }
    if (!(bestF < kInf / 2)) return out;

    const double gridStep = (kInvBetaMax - kInvBetaMin) / nGrid;
    double lo = std::max(kInvBetaMin, bestX - gridStep);
    double hi = std::min(kInvBetaMax, bestX + gridStep);
    constexpr double gr = 0.6180339887498948482;
    double c = hi - gr * (hi - lo);
    double d = lo + gr * (hi - lo);
    double fc = profileSSE(e, c, zeta);
    double fd = profileSSE(e, d, zeta);
    for (int it = 0; it < 24; ++it) {
        if (fc < fd) {
            hi = d; d = c; fd = fc; c = hi - gr * (hi - lo); fc = profileSSE(e, c, zeta);
        } else {
            lo = c; c = d; fc = fd; d = lo + gr * (hi - lo); fd = profileSSE(e, d, zeta);
        }
    }
    double x = 0.5 * (lo + hi);
    double f = profileSSE(e, x, zeta);
    if (bestF < f) { x = bestX; f = bestF; }

    out.sse = profileSSE(e, x, zeta, &out.offset, &out.q, &out.residual);
    out.ok = std::isfinite(out.sse) && out.sse < kInf / 2;
    out.invBeta = x;
    out.beta = 1.0 / x;
    return out;
}

std::string betaBin(double b) {
    if (b < 0.55) return "b0.40_0.55";
    if (b < 0.65) return "b0.55_0.65";
    if (b < 0.75) return "b0.65_0.75";
    if (b < 0.85) return "b0.75_0.85";
    return "b0.85_0.90";
}

std::string edepBin(double e) {
    if (e < 3.5) return "e0_3.5";
    if (e < 4.0) return "e3.5_4.0";
    if (e < 4.5) return "e4.0_4.5";
    if (e < 5.0) return "e4.5_5.0";
    if (e < 5.5) return "e5.0_5.5";
    if (e < 6.5) return "e5.5_6.5";
    return "e6.5_inf";
}

std::vector<std::string> groupKeys(const PreparedEvent &e) {
    const std::string z = "Z" + std::to_string(e.Z);
    return {"all", z, z + "_" + betaBin(e.betaSeed), z + "_" + edepBin(e.edepUpperPerZ2)};
}


std::vector<std::string> groupKeysWithQuality(const PreparedEvent &e, const FitResult &lin) {
    auto keys = groupKeys(e);
    const double chi2 = lin.sse / (kSigmaT * kSigmaT);
    const std::string z = "Z" + std::to_string(e.Z);
    if (chi2 < 10) {
        keys.push_back("all_core10");
        keys.push_back(z + "_core10");
        keys.push_back(z + "_" + betaBin(e.betaSeed) + "_core10");
    }
    if (chi2 < 50) {
        keys.push_back("all_core50");
        keys.push_back(z + "_core50");
    }
    return keys;
}

double quantileSorted(const std::vector<float> &v, double p) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    const double x = p * (v.size() - 1);
    const size_t i = static_cast<size_t>(x);
    const double f = x - i;
    return v[i] * (1.0 - f) + v[std::min(i + 1, v.size() - 1)] * f;
}

struct DeltaStats {
    long long n{0};
    long long better{0};
    long double sum{0};
    long double sum2{0};
    long double sumClip{0};
    long double sumSSE0{0};
    long double sumSSE1{0};
    std::vector<float> values;

    void add(double deltaChi2, double sse0, double sse1, bool keepValues = true) {
        ++n;
        if (deltaChi2 < 0) ++better;
        sum += deltaChi2;
        sum2 += deltaChi2 * deltaChi2;
        sumClip += std::clamp(deltaChi2, -25.0, 25.0);
        sumSSE0 += sse0;
        sumSSE1 += sse1;
        if (keepValues) values.push_back(static_cast<float>(std::clamp(deltaChi2, -1e4, 1e4)));
    }
};

struct DeltaSummary {
    long long n{};
    double mean{}, rms{}, clippedMean{}, median{}, q16{}, q84{}, q05{}, q95{}, fracBetter{};
    double meanSSE0{}, meanSSE1{};
};

DeltaSummary summarize(DeltaStats &s) {
    DeltaSummary x;
    x.n = s.n;
    if (!s.n) return x;
    x.mean = static_cast<double>(s.sum / s.n);
    x.rms = std::sqrt(std::max(0.0, static_cast<double>(s.sum2 / s.n) - x.mean * x.mean));
    x.clippedMean = static_cast<double>(s.sumClip / s.n);
    x.fracBetter = static_cast<double>(s.better) / s.n;
    x.meanSSE0 = static_cast<double>(s.sumSSE0 / s.n);
    x.meanSSE1 = static_cast<double>(s.sumSSE1 / s.n);
    if (!s.values.empty()) {
        std::sort(s.values.begin(), s.values.end());
        x.q05 = quantileSorted(s.values, .05);
        x.q16 = quantileSorted(s.values, .16);
        x.median = quantileSorted(s.values, .50);
        x.q84 = quantileSorted(s.values, .84);
        x.q95 = quantileSorted(s.values, .95);
    }
    return x;
}

struct ShiftStats {
    long long n{0};
    long double sumInv{0}, sumInv2{0}, sumBeta{0}, sumBeta2{0};
    std::vector<float> invValues;
    std::vector<float> betaValues;

    void add(double dInv, double dBeta) {
        ++n;
        sumInv += dInv;
        sumInv2 += dInv * dInv;
        sumBeta += dBeta;
        sumBeta2 += dBeta * dBeta;
        invValues.push_back(static_cast<float>(dInv));
        betaValues.push_back(static_cast<float>(dBeta));
    }
};

struct CurvatureStats {
    long long n{0};
    long double sumRD{0};
    long double sumDD{0};
    long double sumRR{0};
    long double sumE{0};
    long double sumBeta{0};

    void add(double rd, double dd, double rr, double e, double beta) {
        if (!(dd > 0) || !std::isfinite(rd + dd + rr)) return;
        ++n; sumRD += rd; sumDD += dd; sumRR += rr; sumE += e; sumBeta += beta;
    }
};

bool curvatureTerms(const PreparedEvent &e, const FitResult &lin,
                    double &rd, double &dd, double &rr) {
    if (!lin.ok) return false;
    std::array<double, 4> q1{};
    if (!modelTimes(e, lin.invBeta, 1.0, q1)) return false;
    std::array<double, 4> d{};
    for (int i = 0; i < 4; ++i) d[i] = q1[i] - lin.q[i];

    // Remove the two nuisance directions: a common time shift and a linear slope in path.
    double meanS = 0, meanD = 0;
    for (int i = 0; i < 4; ++i) { meanS += e.s[i]; meanD += d[i]; }
    meanS /= 4.0; meanD /= 4.0;
    double num = 0, den = 0;
    for (int i = 0; i < 4; ++i) {
        num += (e.s[i] - meanS) * (d[i] - meanD);
        den += (e.s[i] - meanS) * (e.s[i] - meanS);
    }
    if (!(den > 0)) return false;
    const double slope = num / den;
    rd = dd = rr = 0;
    for (int i = 0; i < 4; ++i) {
        const double dp = d[i] - (meanD + slope * (e.s[i] - meanS));
        rd += lin.residual[i] * dp;
        dd += dp * dp;
        rr += lin.residual[i] * lin.residual[i];
    }
    return dd > 1e-16;
}

struct ScanCell {
    DeltaStats stats;
};

struct BranchReader {
    RawEvent r;
    explicit BranchReader(TTree *t) {
        t->SetBranchStatus("*", 0);
        for (const char *b : {"run", "event", "tof_betah", "tof_tl", "tof_edep", "tof_leng", "tk_q", "tk_qin", "tof_qs"})
            t->SetBranchStatus(b, 1);
        t->SetBranchAddress("run", &r.run);
        t->SetBranchAddress("event", &r.event);
        t->SetBranchAddress("tof_betah", &r.tofBeta);
        t->SetBranchAddress("tof_tl", r.time);
        t->SetBranchAddress("tof_edep", r.edep);
        t->SetBranchAddress("tof_leng", r.length);
        t->SetBranchAddress("tk_q", r.tkQ);
        t->SetBranchAddress("tk_qin", r.tkQin);
        t->SetBranchAddress("tof_qs", &r.tofQs);
    }
};

bool passPrescale(const PreparedEvent &e, int prescale) {
    if (prescale <= 1) return true;
    uint64_t x = (static_cast<uint64_t>(e.run) << 32) ^ e.event;
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return static_cast<int>(x % prescale) == 0;
}

std::vector<std::string> readList(const std::string &path, int maxFiles) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open input list: " + path);
    std::vector<std::string> files;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        files.push_back(line);
        if (maxFiles > 0 && static_cast<int>(files.size()) >= maxFiles) break;
    }
    return files;
}


struct CalibrationResult {
    std::map<int, std::array<double, 4>> residualMedian;
    std::map<int, long long> counts;
};

CalibrationResult computeHighBetaCalibration(const std::vector<std::string> &files,
                                             int begin, int end) {
    std::map<int, std::array<std::vector<float>, 4>> values;
    CalibrationResult result;
    for (int fi = begin; fi < end; ++fi) {
        TFile f(files[fi].c_str(), "READ");
        auto *t = dynamic_cast<TTree *>(f.Get("amstreea"));
        if (!t) continue;
        BranchReader br(t);
        for (Long64_t i = 0; i < t->GetEntries(); ++i) {
            t->GetEntry(i);
            PreparedEvent e;
            if (!prepareRange(br.r, e, 0.97, 1.03)) continue;
            const auto lin = fitLinear(e, 0.6, 1.4);
            if (!lin.ok || lin.sse / (kSigmaT * kSigmaT) > 50) continue;
            ++result.counts[e.Z];
            for (int l = 0; l < 4; ++l) values[e.Z][l].push_back(static_cast<float>(lin.residual[l]));
        }
    }
    for (int z : kCharges) {
        std::array<double, 4> c{};
        for (int l = 0; l < 4; ++l) {
            auto &v = values[z][l];
            if (v.empty()) continue;
            std::sort(v.begin(), v.end());
            c[l] = quantileSorted(v, .5);
        }
        result.residualMedian[z] = c;
    }
    return result;
}

void applyCalibration(PreparedEvent &e, const CalibrationResult &cal) {
    auto it = cal.residualMedian.find(e.Z);
    if (it == cal.residualMedian.end()) return;
    for (int i = 0; i < 4; ++i) e.t[i] -= it->second[i];
}

void writeCalibrationCSV(const fs::path &path,
                         const CalibrationResult &train,
                         const CalibrationResult &validation) {
    std::ofstream out(path);
    out << "split,Z,n,residual_S1_ps,residual_S2_ps,residual_S3_ps,residual_S4_ps\n";
    out << std::setprecision(10);
    for (const auto *item : {&train, &validation}) {
        const std::string split = item == &train ? "train" : "validation";
        for (int z : kCharges) {
            auto c = item->residualMedian.count(z) ? item->residualMedian.at(z) : std::array<double,4>{};
            const long long n = item->counts.count(z) ? item->counts.at(z) : 0;
            out << split << ',' << z << ',' << n;
            for (double x : c) out << ',' << 1000.0 * x;
            out << '\n';
        }
    }
}

void writeScanCSV(const fs::path &path,
                  const std::vector<double> &zetas,
                  std::map<std::string, std::vector<ScanCell>> &scan) {
    std::ofstream out(path);
    out << "group,zeta,n,mean_delta_chi2,rms_delta_chi2,clipped_mean_delta_chi2,median_delta_chi2,q16,q84,frac_better,mean_sse_linear_ns2,mean_sse_model_ns2\n";
    out << std::setprecision(10);
    for (auto &[group, cells] : scan) {
        for (size_t iz = 0; iz < zetas.size(); ++iz) {
            auto s = summarize(cells[iz].stats);
            out << group << ',' << zetas[iz] << ',' << s.n << ',' << s.mean << ',' << s.rms << ',' << s.clippedMean << ','
                << s.median << ',' << s.q16 << ',' << s.q84 << ',' << s.fracBetter << ','
                << s.meanSSE0 << ',' << s.meanSSE1 << '\n';
        }
    }
}

struct BestZeta { double raw{0}, clipped{0}, median{0}; };

std::map<std::string, BestZeta> chooseBest(const std::vector<double> &zetas,
                                           std::map<std::string, std::vector<ScanCell>> &scan) {
    std::map<std::string, BestZeta> best;
    for (auto &[group, cells] : scan) {
        double br = kInf, bc = kInf, bm = kInf;
        BestZeta bz;
        for (size_t iz = 0; iz < zetas.size(); ++iz) {
            auto s = summarize(cells[iz].stats);
            if (!s.n) continue;
            if (s.mean < br) { br = s.mean; bz.raw = zetas[iz]; }
            if (s.clippedMean < bc) { bc = s.clippedMean; bz.clipped = zetas[iz]; }
            if (s.median < bm) { bm = s.median; bz.median = zetas[iz]; }
        }
        best[group] = bz;
    }
    return best;
}

void writeCurvatureCSV(const fs::path &path,
                       const std::map<std::string, CurvatureStats> &curv) {
    std::ofstream out(path);
    out << "group,n,effective_zeta,stat_error_assuming_sigma,aggregate_projection_sigma,corr_residual_template,rms_template_ps,rms_linear_residual_ps,mean_beta,mean_upper_edep_per_z2_mev\n";
    out << std::setprecision(10);
    for (const auto &[group, s] : curv) {
        if (!s.n || !(s.sumDD > 0)) continue;
        const double alpha = static_cast<double>(s.sumRD / s.sumDD);
        const double err = kSigmaT / std::sqrt(static_cast<double>(s.sumDD));
        const double projectionSigma = static_cast<double>(s.sumRD / (kSigmaT * std::sqrt(s.sumDD)));
        const double corr = static_cast<double>(s.sumRD / std::sqrt(s.sumDD * s.sumRR));
        const double rmsTemplatePs = 1000.0 * std::sqrt(static_cast<double>(s.sumDD / s.n / 4.0));
        const double rmsResidualPs = 1000.0 * std::sqrt(static_cast<double>(s.sumRR / s.n / 4.0));
        out << group << ',' << s.n << ',' << alpha << ',' << err << ',' << projectionSigma << ',' << corr << ',' << rmsTemplatePs << ',' << rmsResidualPs << ','
            << static_cast<double>(s.sumBeta / s.n) << ',' << static_cast<double>(s.sumE / s.n) << '\n';
    }
}

void addCurvature(std::map<std::string, CurvatureStats> &curv,
                  const PreparedEvent &e, const FitResult &lin) {
    double rd, dd, rr;
    if (!curvatureTerms(e, lin, rd, dd, rr)) return;
    for (const auto &g : groupKeysWithQuality(e, lin)) curv[g].add(rd, dd, rr, e.edepUpperPerZ2, e.betaSeed);
}

struct ModelSpec {
    std::string name;
    double commonZeta{0};
    bool perZ{false};
    std::map<int, double> zetaByZ;
    double zeta(int z) const {
        if (!perZ) return commonZeta;
        auto it = zetaByZ.find(z);
        return it == zetaByZ.end() ? commonZeta : it->second;
    }
};

} // namespace

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " input.list output_dir [train_files=65] [train_prescale=5] [max_files=0]"
                  << " [reference=before_tof|center] [global_zeta_Z2=1.9] [global_zeta_Z6=1.9] [global_zeta_Z8=1.9]\n";
        return 1;
    }
    const std::string listPath = argv[1];
    const fs::path outputDir = argv[2];
    int trainFiles = argc > 3 ? std::stoi(argv[3]) : 65;
    const int trainPrescale = argc > 4 ? std::stoi(argv[4]) : 5;
    const int maxFiles = argc > 5 ? std::stoi(argv[5]) : 0;
    const std::string reference = argc > 6 ? argv[6] : "before_tof";
    if (reference != "before_tof" && reference != "center")
        throw std::runtime_error("reference must be before_tof or center");
    gCenterReference = reference == "center";
    const std::map<int, double> mcGlobalZeta{{2, argc > 7 ? std::stod(argv[7]) : 1.9},
                                             {6, argc > 8 ? std::stod(argv[8]) : 1.9},
                                             {8, argc > 9 ? std::stod(argv[9]) : 1.9}};
    fs::create_directories(outputDir);

    auto files = readList(listPath, maxFiles);
    if (files.size() < 2) throw std::runtime_error("Need at least two files");
    trainFiles = std::clamp(trainFiles, 1, static_cast<int>(files.size()) - 1);

    std::cerr << "Computing independent high-beta layer-residual calibration...\n";
    const auto trainCalibration = computeHighBetaCalibration(files, 0, trainFiles);
    const auto validationCalibration = computeHighBetaCalibration(files, trainFiles, files.size());
    writeCalibrationCSV(outputDir / "high_beta_layer_calibration.csv", trainCalibration, validationCalibration);

    std::vector<double> zetas;
    for (int i = -10; i <= 120; ++i) zetas.push_back(i / 10.0);
    std::map<std::string, std::vector<ScanCell>> scan;
    for (const auto &g : {std::string("all"), std::string("Z2"), std::string("Z6"), std::string("Z8")})
        scan[g].resize(zetas.size());
    std::map<std::string, CurvatureStats> trainCurv;

    long long trainRead = 0, trainSelected = 0, trainScanned = 0;
    std::cerr << "Training files: [0," << trainFiles << "), validation files: [" << trainFiles << ',' << files.size() << ")\n";
    for (int fi = 0; fi < trainFiles; ++fi) {
        TFile f(files[fi].c_str(), "READ");
        auto *t = dynamic_cast<TTree *>(f.Get("amstreea"));
        if (!t) { std::cerr << "Skip missing tree: " << files[fi] << '\n'; continue; }
        BranchReader br(t);
        const Long64_t n = t->GetEntries();
        trainRead += n;
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            PreparedEvent e;
            if (!prepare(br.r, e)) continue;
            applyCalibration(e, trainCalibration);
            ++trainSelected;
            const auto lin = fitLinear(e);
            if (!lin.ok) continue;
            addCurvature(trainCurv, e, lin);
            if (!passPrescale(e, trainPrescale)) continue;
            ++trainScanned;
            const std::string zg = "Z" + std::to_string(e.Z);
            for (size_t iz = 0; iz < zetas.size(); ++iz) {
                const auto fit = fitAtZeta(e, zetas[iz], &lin);
                if (!fit.ok) continue;
                const double dchi2 = (fit.sse - lin.sse) / (kSigmaT * kSigmaT);
                scan["all"][iz].stats.add(dchi2, lin.sse, fit.sse);
                scan[zg][iz].stats.add(dchi2, lin.sse, fit.sse);
                const double linChi2 = lin.sse / (kSigmaT * kSigmaT);
                if (linChi2 < 10) {
                    const std::string zbeta = zg + "_" + betaBin(e.betaSeed) + "_core10";
                    for (const auto &key : {std::string("all_core10"), zg + "_core10", zbeta}) {
                        if (scan[key].empty()) scan[key].resize(zetas.size());
                        scan[key][iz].stats.add(dchi2, lin.sse, fit.sse);
                    }
                    if (e.betaSeed < 0.75) {
                        for (const auto &key : {std::string("all_lowbeta_core10"), zg + "_lowbeta_core10"}) {
                            if (scan[key].empty()) scan[key].resize(zetas.size());
                            scan[key][iz].stats.add(dchi2, lin.sse, fit.sse);
                        }
                    }
                }
            }
        }
        if ((fi + 1) % 10 == 0 || fi + 1 == trainFiles)
            std::cerr << "train file " << (fi + 1) << '/' << trainFiles << ", selected=" << trainSelected << ", scan=" << trainScanned << '\n';
    }

    writeScanCSV(outputDir / "train_zeta_scan.csv", zetas, scan);
    writeCurvatureCSV(outputDir / "train_curvature_projection.csv", trainCurv);
    auto best = chooseBest(zetas, scan);
    {
        std::ofstream out(outputDir / "chosen_zetas.txt");
        out << "train_files=" << trainFiles << "\ntrain_prescale=" << trainPrescale << "\n";
        out << "train_entries_read=" << trainRead << "\ntrain_selected=" << trainSelected << "\ntrain_scanned=" << trainScanned << "\n";
        for (const auto &[g, b] : best)
            out << g << " raw=" << b.raw << " clipped=" << b.clipped << " median=" << b.median << '\n';
    }

    std::vector<ModelSpec> models;
    models.push_back({"zeta_1.0", 1.0, false, {}});
    models.push_back({"zeta_1.9", 1.9, false, {}});
    models.push_back({"train_common_raw", best["all"].raw, false, {}});
    models.push_back({"train_common_clipped", best["all"].clipped, false, {}});
    ModelSpec perZ{"train_perZ_clipped", best["all"].clipped, true, {}};
    for (int z : kCharges) perZ.zetaByZ[z] = best["Z" + std::to_string(z)].clipped;
    models.push_back(perZ);
    models.push_back({"train_lowbeta_common_clipped", best["all_lowbeta_core10"].clipped, false, {}});
    ModelSpec perZLow{"train_lowbeta_perZ_clipped", best["all_lowbeta_core10"].clipped, true, {}};
    for (int z : kCharges) perZLow.zetaByZ[z] = best["Z" + std::to_string(z) + "_lowbeta_core10"].clipped;
    models.push_back(perZLow);
    ModelSpec mcGlobal{"mc_center_global", 1.9, true, mcGlobalZeta};
    models.push_back(mcGlobal);

    std::map<std::string, DeltaStats> validation;
    std::map<std::string, ShiftStats> validationShift;
    std::map<std::pair<int, std::string>, DeltaStats> validationByFile;
    std::map<std::string, CurvatureStats> validationCurv;
    long long valRead = 0, valSelected = 0;
    for (int fi = trainFiles; fi < static_cast<int>(files.size()); ++fi) {
        TFile f(files[fi].c_str(), "READ");
        auto *t = dynamic_cast<TTree *>(f.Get("amstreea"));
        if (!t) { std::cerr << "Skip missing tree: " << files[fi] << '\n'; continue; }
        BranchReader br(t);
        const Long64_t n = t->GetEntries();
        valRead += n;
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            PreparedEvent e;
            if (!prepare(br.r, e)) continue;
            applyCalibration(e, validationCalibration);
            ++valSelected;
            const auto lin = fitLinear(e);
            if (!lin.ok) continue;
            addCurvature(validationCurv, e, lin);
            const auto groups = groupKeysWithQuality(e, lin);
            for (const auto &model : models) {
                const double zeta = model.zeta(e.Z);
                const auto fit = fitAtZeta(e, zeta, &lin);
                if (!fit.ok) continue;
                const double dchi2 = (fit.sse - lin.sse) / (kSigmaT * kSigmaT);
                for (const auto &g : groups) {
                    const std::string key = model.name + "|" + g;
                    validation[key].add(dchi2, lin.sse, fit.sse);
                    if (model.name == "mc_center_global") {
                        validationShift[g].add(fit.invBeta - lin.invBeta, fit.beta - lin.beta);
                        if (lin.sse / (kSigmaT * kSigmaT) < 10 &&
                            (g == "all_core10" || g == "Z2_core10" || g == "Z6_core10" || g == "Z8_core10" ||
                             g.find("_core10") != std::string::npos)) {
                            validationByFile[{fi, g}].add(dchi2, lin.sse, fit.sse);
                        }
                    }
                }
            }
        }
        if ((fi - trainFiles + 1) % 10 == 0 || fi + 1 == static_cast<int>(files.size()))
            std::cerr << "validation file " << (fi - trainFiles + 1) << '/' << (files.size() - trainFiles)
                      << ", selected=" << valSelected << '\n';
    }

    {
        std::ofstream out(outputDir / "validation_delta_chi2.csv");
        out << "model,group,zeta_definition,n,mean_delta_chi2,rms_delta_chi2,clipped_mean_delta_chi2,median_delta_chi2,q05,q16,q84,q95,frac_better,mean_sse_linear_ns2,mean_sse_model_ns2\n";
        out << std::setprecision(10);
        for (auto &[key, stats] : validation) {
            const auto pos = key.find('|');
            const std::string modelName = key.substr(0, pos);
            const std::string group = key.substr(pos + 1);
            const auto &m = *std::find_if(models.begin(), models.end(), [&](const ModelSpec &x) { return x.name == modelName; });
            std::ostringstream zd;
            if (m.perZ) {
                zd << "perZ:";
                for (int z : kCharges) zd << "Z" << z << '=' << m.zeta(z) << ';';
            } else zd << m.commonZeta;
            auto s = summarize(stats);
            out << modelName << ',' << group << ',' << zd.str() << ',' << s.n << ',' << s.mean << ',' << s.rms << ',' << s.clippedMean << ','
                << s.median << ',' << s.q05 << ',' << s.q16 << ',' << s.q84 << ',' << s.q95 << ',' << s.fracBetter << ','
                << s.meanSSE0 << ',' << s.meanSSE1 << '\n';
        }
    }
    {
        std::ofstream out(outputDir / "validation_beta_shift_mc_center_global.csv");
        out << "group,n,mean_delta_invbeta,rms_delta_invbeta,median_delta_invbeta,q16_delta_invbeta,q84_delta_invbeta,mean_delta_beta,rms_delta_beta,median_delta_beta\n";
        out << std::setprecision(10);
        for (auto &[group, st] : validationShift) {
            if (!st.n) continue;
            std::sort(st.invValues.begin(), st.invValues.end());
            std::sort(st.betaValues.begin(), st.betaValues.end());
            const double mi = static_cast<double>(st.sumInv / st.n);
            const double mb = static_cast<double>(st.sumBeta / st.n);
            const double ri = std::sqrt(std::max(0.0, static_cast<double>(st.sumInv2 / st.n) - mi * mi));
            const double rb = std::sqrt(std::max(0.0, static_cast<double>(st.sumBeta2 / st.n) - mb * mb));
            out << group << ',' << st.n << ',' << mi << ',' << ri << ',' << quantileSorted(st.invValues, .5) << ','
                << quantileSorted(st.invValues, .16) << ',' << quantileSorted(st.invValues, .84) << ','
                << mb << ',' << rb << ',' << quantileSorted(st.betaValues, .5) << '\n';
        }
    }
    {
        std::ofstream out(outputDir / "validation_file_stability_mc_center_global.csv");
        out << "file_index,file,group,n,mean_delta_chi2,rms_delta_chi2,median_delta_chi2,frac_better,mean_sse_linear_ns2,mean_sse_model_ns2\n";
        out << std::setprecision(10);
        for (auto &[key, st] : validationByFile) {
            auto summary = summarize(st);
            const int fi = key.first;
            out << fi << ',' << files[fi] << ',' << key.second << ',' << summary.n << ',' << summary.mean << ',' << summary.rms << ','
                << summary.median << ',' << summary.fracBetter << ',' << summary.meanSSE0 << ',' << summary.meanSSE1 << '\n';
        }
    }
    writeCurvatureCSV(outputDir / "validation_curvature_projection.csv", validationCurv);

    {
        std::ofstream out(outputDir / "run_summary.txt");
        out << "files_total=" << files.size() << '\n';
        out << "train_files=" << trainFiles << '\n';
        out << "validation_files=" << files.size() - trainFiles << '\n';
        out << "train_entries_read=" << trainRead << '\n';
        out << "train_selected=" << trainSelected << '\n';
        out << "train_scanned=" << trainScanned << '\n';
        out << "validation_entries_read=" << valRead << '\n';
        out << "validation_selected=" << valSelected << '\n';
        out << "reference=" << reference << '\n';
        for (int z : kCharges) out << "mc_global_zeta_Z" << z << '=' << mcGlobalZeta.at(z) << '\n';
    }

    std::cerr << "Done. Results in " << outputDir << '\n';
    return 0;
}
