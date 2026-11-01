#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

namespace {
constexpr double kC = 29.9792458;
constexpr double kInf = 1e100;
constexpr int kBins = 14;
constexpr double kBetaMin = 0.25;
constexpr double kBetaMax = 0.95;
using Point = std::array<double, 3>;
struct MomentumPoint {
  double z{}, p{};
};
struct Event {
  double mcMass{}, fitMass{}, betaCenter{}, betaCheckpoint10{}, betaSeed{},
      centerTruthEnergy{};
  std::array<double, 4> timeMeasured{}, timeTruth{}, tofEdepGeV{},
      mtofEdepGeV{}, path{};
  std::array<Point, 4> tofPos{};
  std::array<double, 9> tkEdepGeV{};
  std::array<Point, 9> tkPos{};
  std::vector<MomentumPoint> momentumPoints;
};
struct FitResult {
  bool ok{ false };
  double invBeta{ NAN }, sse{ kInf };
};
struct Stat {
  std::vector<double> v;
  void add(double x) {
    if (std::isfinite(x))
      v.push_back(x);
  }
};
struct RowStats {
  long long n{};
  double mean{ NAN }, median{ NAN }, sigma68{ NAN }, rms{ NAN }, q05{ NAN },
      q95{ NAN }, coreMu{ NAN }, coreSigma{ NAN };
};

double distance(const Point &a, const Point &b) {
  double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}
double quantile(std::vector<double> v, double p) {
  if (v.empty())
    return NAN;
  std::sort(v.begin(), v.end());
  double x = p * (v.size() - 1);
  size_t i = (size_t)x;
  double f = x - i;
  return v[i] * (1 - f) + v[std::min(i + 1, v.size() - 1)] * f;
}
RowStats summarize(const Stat &s) {
  RowStats r;
  r.n = s.v.size();
  if (s.v.empty())
    return r;
  long double sum = 0, sum2 = 0;
  for (double x : s.v) {
    sum += x;
    sum2 += x * x;
  }
  r.mean = (double)(sum / s.v.size());
  r.rms =
      std::sqrt(std::max(0.0, (double)(sum2 / s.v.size()) - r.mean * r.mean));
  r.median = quantile(s.v, .5);
  r.sigma68 = .5 * (quantile(s.v, .84) - quantile(s.v, .16));
  r.q05 = quantile(s.v, .05);
  r.q95 = quantile(s.v, .95);
  if (s.v.size() >= 300) {
    static unsigned long counter = 0;
    TH1D h(("hCore" + std::to_string(counter)).c_str(), "", 800, -.5, 1.5);
    for (double x : s.v)
      h.Fill(x);
    double q[3]{}, p[3]{ .25, .5, .75 };
    h.GetQuantiles(3, q, p);
    double robustSigma = .741301109252801 * (q[2] - q[0]);
    if (robustSigma > 0 && std::isfinite(robustSigma)) {
      double lo = std::max(-.5, q[1] - 3 * robustSigma);
      double hi = std::min(1.5, q[1] + 3 * robustSigma);
      TF1 f(("fCore" + std::to_string(counter++)).c_str(), "gaus", lo, hi);
      f.SetParameters(h.GetMaximum(), q[1], robustSigma);
      TFitResultPtr fit = h.Fit(&f, "SQNR");
      if (fit.Get() && fit->Status() == 0 && std::isfinite(f.GetParameter(1)) &&
          std::isfinite(f.GetParameter(2)) && std::abs(f.GetParameter(2)) < .08 &&
          fit->ParError(1) < .02 &&
          fit->ParError(2) < .5 * std::abs(f.GetParameter(2))) {
        r.coreMu = f.GetParameter(1);
        r.coreSigma = std::abs(f.GetParameter(2));
      }
    }
  }
  return r;
}
int betaBin(double b) {
  if (!(b >= kBetaMin && b < kBetaMax))
    return -1;
  int i = (int)((b - kBetaMin) / (kBetaMax - kBetaMin) * kBins);
  return std::max(0, std::min(kBins - 1, i));
}
std::vector<std::string> readList(const std::string &p, int maxFiles) {
  std::ifstream in(p);
  std::vector<std::string> o;
  std::string l;
  while (std::getline(in, l)) {
    if (l.empty() || l[0] == '#')
      continue;
    o.push_back(l);
    if (maxFiles > 0 && (int)o.size() >= maxFiles)
      break;
  }
  return o;
}
std::string basenameOf(const std::string &p) {
  auto x = p.find_last_of('/');
  return x == std::string::npos ? p : p.substr(x + 1);
}
void normalize(std::vector<MomentumPoint> &p) {
  std::sort(p.begin(), p.end(), [](auto &a, auto &b) { return a.z > b.z; });
  std::vector<MomentumPoint> u;
  for (auto &x : p) {
    if (u.empty() || std::abs(x.z - u.back().z) > 1e-5)
      u.push_back(x);
    else {
      u.back().p = .5 * (u.back().p + x.p);
      u.back().z = .5 * (u.back().z + x.z);
    }
  }
  p.swap(u);
}
double momentumAt(const std::vector<MomentumPoint> &p, double z) {
  if (p.empty() || z > p.front().z || z < p.back().z)
    return NAN;
  for (size_t i = 0; i + 1 < p.size(); ++i) {
    if (p[i].z + 1e-9 < z || z < p[i + 1].z - 1e-9)
      continue;
    double dz = p[i].z - p[i + 1].z;
    if (std::abs(dz) < 1e-12)
      return .5 * (p[i].p + p[i + 1].p);
    double f = (p[i].z - z) / dz;
    return p[i].p * (1 - f) + p[i + 1].p * f;
  }
  return NAN;
}
double truthEnergyAt(const Event &e, double z) {
  double p = momentumAt(e.momentumPoints, z);
  return p > 0 ? std::sqrt(p * p + e.mcMass * e.mcMass) : NAN;
}
double candidateInvAtZ(const Event &e, double z, double ec) {
  double et = truthEnergyAt(e, z);
  if (!std::isfinite(et))
    return NAN;
  double en = ec + (et - e.centerTruthEnergy);
  if (!(en > e.mcMass))
    return NAN;
  double a = 1 - e.mcMass * e.mcMass / (en * en);
  return a > 0 && a < 1 ? 1 / std::sqrt(a) : NAN;
}
double averageCandidateInv(const Event &e, double z1, double z2, double ec) {
  if (z1 < z2)
    std::swap(z1, z2);
  if (!(z1 > z2))
    return NAN;
  std::vector<double> z{ z1, z2 };
  for (auto &p : e.momentumPoints)
    if (p.z < z1 && p.z > z2)
      z.push_back(p.z);
  std::sort(z.begin(), z.end(), std::greater<double>());
  z.erase(std::unique(z.begin(), z.end(), [](double a, double b) {
            return std::abs(a - b) < 1e-6;
          }),
          z.end());
  double in = 0;
  for (size_t i = 0; i + 1 < z.size(); ++i) {
    double a = z[i], b = z[i + 1], m = .5 * (a + b);
    double fa = candidateInvAtZ(e, a, ec), fm = candidateInvAtZ(e, m, ec),
           fb = candidateInvAtZ(e, b, ec);
    if (!std::isfinite(fa + fm + fb))
      return NAN;
    in += (a - b) * (fa + 4 * fm + fb) / 6;
  }
  return in / (z1 - z2);
}
bool oracleTimes(const Event &e, double inv, std::array<double, 4> &q) {
  if (!(inv > 1 && inv < 5))
    return false;
  double b = 1 / inv, ec = e.mcMass / std::sqrt(1 - b * b);
  q.fill(0);
  for (int i = 1; i < 4; ++i) {
    double a = averageCandidateInv(e, e.tofPos[i - 1][2], e.tofPos[i][2], ec),
           ds = e.path[i] - e.path[i - 1];
    if (!std::isfinite(a) || !(ds > 0 && ds < 300))
      return false;
    q[i] = q[i - 1] + ds * a / kC;
  }
  return true;
}
bool truthTimes(Event &e) {
  return oracleTimes(e, 1 / e.betaCenter, e.timeTruth);
}

bool currentTimes(const Event &e, double inv, double zeta, double trackerScale,
                  bool useMtof, std::array<double, 4> &q) {
  if (!(inv > 1 && inv < 5))
    return false;
  double b = 1 / inv, m = e.fitMass, m2 = m * m, en = m / std::sqrt(1 - b * b);
  if (!(en > m))
    return false;
  auto invAt = [&](double x)->double {
    if (!(x > m) || !std::isfinite(x))
      return NAN;
    double a = 1 - m2 / (x * x);
    return a > 0 && a < 1 ? 1 / std::sqrt(a) : NAN;
  };
  auto seg = [&](double l, double x) {
    double y = invAt(x);
    return std::isfinite(y) ? l * y / kC : NAN;
  };
  const auto &dep = useMtof ? e.mtofEdepGeV : e.tofEdepGeV;
  auto scalePath = [&](const Point &a, const Point &b, int l1, int l2,
                       double target) {
    double g = distance(a, e.tkPos[l1]);
    int st = l1 < l2 ? 1 : -1;
    for (int l = l1; l != l2; l += st)
      g += distance(e.tkPos[l], e.tkPos[l + st]);
    g += distance(e.tkPos[l2], b);
    return g > 0 ? std::abs(target) / g : 0;
  };
  double den = e.tkPos[4][2] - e.tkPos[5][2];
  if (!(den > 0 && e.tkPos[4][2] > 0 && e.tkPos[5][2] < 0))
    return false;
  double f = e.tkPos[4][2] / den;
  Point c{};
  for (int j = 0; j < 3; ++j)
    c[j] = e.tkPos[4][j] + f * (e.tkPos[5][j] - e.tkPos[4][j]);
  double us = scalePath(c, e.tofPos[1], 4, 1, e.path[1]),
         ls = scalePath(c, e.tofPos[2], 5, 7, e.path[2]);
  if (!(us > 0 && ls > 0))
    return false;
  q.fill(NAN);
  double ue = en, ut = 0;
  Point cur = c;
  for (int l = 4; l >= 1; --l) {
    double dt = seg(distance(cur, e.tkPos[l]) * us, ue);
    if (!std::isfinite(dt))
      return false;
    ut += dt;
    ue += trackerScale * e.tkEdepGeV[l];
    cur = e.tkPos[l];
  }
  double dt = seg(distance(cur, e.tofPos[1]) * us, ue);
  if (!std::isfinite(dt))
    return false;
  ut += dt;
  q[1] = -ut;
  double le = en, lt = 0;
  cur = c;
  for (int l = 5; l <= 7; ++l) {
    double d = seg(distance(cur, e.tkPos[l]) * ls, le);
    if (!std::isfinite(d))
      return false;
    lt += d;
    le -= trackerScale * e.tkEdepGeV[l];
    cur = e.tkPos[l];
  }
  dt = seg(distance(cur, e.tofPos[2]) * ls, le);
  if (!std::isfinite(dt))
    return false;
  lt += dt;
  q[2] = lt;
  q[0] = q[1] + seg(e.path[0] - e.path[1], ue + zeta * dep[1]);
  q[3] = q[2] + seg(e.path[3] - e.path[2], le - zeta * dep[2]);
  return std::all_of(q.begin(), q.end(),
                     [](double x) { return std::isfinite(x); });
}
bool legacyTimes(const Event &e, double inv, double zeta,
                 std::array<double, 4> &q) {
  if (!(inv > 1 && inv < 5))
    return false;
  double b = 1 / inv, m = e.fitMass, m2 = m * m, en = m / std::sqrt(1 - b * b);
  auto ia = [&](double x)->double {
    if (!(x > m))
      return NAN;
    double a = 1 - m2 / (x * x);
    return a > 0 ? 1 / std::sqrt(a) : NAN;
  };
  q[1] = e.path[1] * inv / kC;
  q[2] = e.path[2] * inv / kC;
  q[0] = q[1] + (e.path[0] - e.path[1]) * ia(en + zeta * e.tofEdepGeV[1]) / kC;
  q[3] = q[2] + (e.path[3] - e.path[2]) * ia(en - zeta * e.tofEdepGeV[2]) / kC;
  return std::all_of(q.begin(), q.end(),
                     [](double x) { return std::isfinite(x); });
}
double profileSSE(const std::array<double, 4> &t,
                  const std::array<double, 4> &q) {
  double a = 0;
  for (int i = 0; i < 4; ++i) {
    if (!std::isfinite(t[i] + q[i]))
      return kInf;
    a += t[i] - q[i];
  }
  a /= 4;
  double s = 0;
  for (int i = 0; i < 4; ++i) {
    double r = t[i] - (a + q[i]);
    s += r * r;
  }
  return s;
}
using Model = std::function<bool(double, std::array<double, 4> &)>;
FitResult fitModel(const std::array<double, 4> &t, const Model &m,
                   double seed) {
  constexpr double lo0 = 1.000001, hi0 = 4;
  constexpr int ng = 80;
  double bx = std::clamp(seed, lo0, hi0), bf = kInf;
  auto ev = [&](double x) {
    std::array<double, 4> q{};
    return m(x, q) ? profileSSE(t, q) : kInf;
  };
  bf = ev(bx);
  for (int j = 0; j <= ng; ++j) {
    double x = lo0 + (hi0 - lo0) * j / ng, f = ev(x);
    if (f < bf) {
      bf = f;
      bx = x;
    }
  }
  if (!(bf < kInf / 2))
    return {};
  double st = (hi0 - lo0) / ng, lo = std::max(lo0, bx - st),
         hi = std::min(hi0, bx + st), gr = .6180339887498948482,
         c = hi - gr * (hi - lo), d = lo + gr * (hi - lo), fc = ev(c),
         fd = ev(d);
  for (int it = 0; it < 40; ++it) {
    if (fc < fd) {
      hi = d;
      d = c;
      fd = fc;
      c = hi - gr * (hi - lo);
      fc = ev(c);
    } else {
      lo = c;
      c = d;
      fc = fd;
      d = lo + gr * (hi - lo);
      fd = ev(d);
    }
  }
  double x = .5 * (lo + hi), f = ev(x);
  FitResult r;
  if (f < bf) {
    r.invBeta = x;
    r.sse = f;
  } else {
    r.invBeta = bx;
    r.sse = bf;
  }
  r.ok = std::isfinite(r.invBeta) && r.sse < kInf / 2;
  return r;
}
FitResult fitLinear(const Event &e, const std::array<double, 4> &t) {
  std::array<double, 4> s{};
  for (int i = 1; i < 4; ++i)
    s[i] = s[i - 1] + e.path[i] - e.path[i - 1];
  double ms = 0, mt = 0;
  for (int i = 0; i < 4; ++i) {
    ms += s[i];
    mt += t[i];
  }
  ms /= 4;
  mt /= 4;
  double n = 0, d = 0;
  for (int i = 0; i < 4; ++i) {
    n += (s[i] - ms) * (t[i] - mt);
    d += (s[i] - ms) * (s[i] - ms);
  }
  FitResult r;
  if (!(d > 0))
    return r;
  r.invBeta = std::clamp(kC * n / d, 1.000001, 4.0);
  std::array<double, 4> q{};
  for (int i = 0; i < 4; ++i)
    q[i] = s[i] * r.invBeta / kC;
  r.sse = profileSSE(t, q);
  r.ok = std::isfinite(r.sse);
  return r;
}
double mcMassFor(int mpar, double mch) {
  return mpar == 69 ? 14.899169 : 2 * mch * .9314941;
}
struct OutputReader {
  TFile f;
  TTree *t{};
  Long64_t index{};
  double mcBeta{}, linear{}, legacy{}, nonlinear{};
  explicit OutputReader(const std::string &p) : f(p.c_str(), "READ") {
    if (!f.IsZombie())
      t = dynamic_cast<TTree *>(f.Get("betaTree"));
    if (t) {
      t->SetBranchAddress("mcBeta", &mcBeta);
      t->SetBranchAddress("linearBeta", &linear);
      t->SetBranchAddress("nonlinearBetaLegacy", &legacy);
      t->SetBranchAddress("nonlinearBeta", &nonlinear);
    }
  }
  bool next() {
    if (!t || index >= t->GetEntries())
      return false;
    t->GetEntry(index++);
    return true;
  }
};
}

int main(int argc, char **argv) {
  if (argc < 9) {
    std::cerr << "usage: bias_decomp input.list output_event_dir Z max_files "
                 "residual.csv timing.csv summary.txt physical.csv "
                 "[events.csv] [common_scale]\n";
    return 1;
  }
  std::string list = argv[1], outDir = argv[2], resPath = argv[5],
              timPath = argv[6], sumPath = argv[7], physPath = argv[8],
              eventPath = argc > 9 ? argv[9] : "";
  std::ofstream eventOut;
  if (!eventPath.empty()) {
    eventOut.open(eventPath);
    eventOut << "file_index,Z,beta_center,bar2,bar3,delta23_ns,raw_delta23_ns,"
                "correction23_ns,logA2,logA3,edep2_per_z2,edep3_per_z2,truth_"
                "dE_upper_MeV,truth_dE_lower_MeV,active_upper_MeV,active_lower_"
                "MeV\n";
  }
  int Z = std::stoi(argv[3]), maxFiles = std::stoi(argv[4]);
  double commonScale = argc > 10 ? std::stod(argv[10]) : 1.9;
  double fitMass = 2 * Z * .9314941;
  auto files = readList(list, maxFiles);
  std::map<std::string, std::array<Stat, kBins> > residual, timing, physical;
  long long raw = 0, selected = 0, valid = 0, paired = 0, mismatch = 0,
            cmFail = 0, ctFail = 0, omFail = 0, otFail = 0;
  for (size_t fi = 0; fi < files.size(); ++fi) {
    TFile f(files[fi].c_str(), "READ");
    auto *t = dynamic_cast<TTree *>(f.Get("amstreea"));
    if (!t)
      continue;
    int mpar = 0;
    float mch = 0, pm[21]{}, co[21][3]{}, tl[4]{}, rt[4][2]{}, aadc[4][2]{},
          ted[4]{}, mtd[4]{}, leng[4]{}, tpos[4][3]{}, ked[9]{}, kpos[9][3]{},
          tb = 0, tbmc = 0;
    int barid[4]{};
    t->SetBranchStatus("*", 0);
    for (const char *b :
         { "mpar",     "mch",       "mevmom1",   "mevcoo1",
           "tof_tl",   "tof_rtime", "tof_Aadc",  "tof_barid",
           "tof_edep", "mtofdep",   "tof_leng",  "tof_pos",
           "tk_edep",  "tk_pos",    "tof_betah", "tof_betahmc" })
      t->SetBranchStatus(b, 1);
    t->SetBranchAddress("mpar", &mpar);
    t->SetBranchAddress("mch", &mch);
    t->SetBranchAddress("mevmom1", pm);
    t->SetBranchAddress("mevcoo1", co);
    t->SetBranchAddress("tof_tl", tl);
    t->SetBranchAddress("tof_rtime", rt);
    t->SetBranchAddress("tof_Aadc", aadc);
    t->SetBranchAddress("tof_barid", barid);
    t->SetBranchAddress("tof_edep", ted);
    t->SetBranchAddress("mtofdep", mtd);
    t->SetBranchAddress("tof_leng", leng);
    t->SetBranchAddress("tof_pos", tpos);
    t->SetBranchAddress("tk_edep", ked);
    t->SetBranchAddress("tk_pos", kpos);
    t->SetBranchAddress("tof_betah", &tb);
    t->SetBranchAddress("tof_betahmc", &tbmc);
    OutputReader out(outDir + "/" + basenameOf(files[fi]));
    for (Long64_t ie = 0; ie < t->GetEntries(); ++ie) {
      ++raw;
      t->GetEntry(ie);
      if (pm[10] == -1000)
        continue;
      bool kc = true;
      for (int l = 1; l <= 7; ++l)
        kc = kc && ked[l] > 0 && std::isfinite(ked[l]) &&
             std::isfinite(kpos[l][2]);
      if (!kc || !(kpos[4][2] > 0) || !(kpos[5][2] < 0))
        continue;
      ++selected;
      bool has = out.next();
      if (has)
        ++paired;
      Event e;
      e.mcMass = mcMassFor(mpar, mch);
      e.fitMass = fitMass;
      e.betaSeed = tb;
      for (int i = 0; i < 21; ++i)
        if (pm[i] > 0 && std::isfinite(pm[i]) && std::isfinite(co[i][2]))
          e.momentumPoints.push_back({ co[i][2], pm[i] });
      normalize(e.momentumPoints);
      double pc = momentumAt(e.momentumPoints, 0);
      if (!(pc > 0))
        continue;
      e.centerTruthEnergy = std::sqrt(pc * pc + e.mcMass * e.mcMass);
      e.betaCenter = pc / e.centerTruthEnergy;
      e.betaCheckpoint10 =
          pm[10] / std::sqrt(pm[10] * pm[10] + e.mcMass * e.mcMass);
      if (has && std::abs(out.mcBeta - e.betaCheckpoint10) > 2e-5)
        ++mismatch;
      bool ok = true;
      for (int i = 0; i < 4; ++i) {
        e.timeMeasured[i] = tl[i];
        e.tofEdepGeV[i] = ted[i] * 1e-3;
        e.mtofEdepGeV[i] = mtd[i] * 1e-3;
        e.path[i] = leng[i];
        for (int c = 0; c < 3; ++c)
          e.tofPos[i][c] = tpos[i][c];
        ok = ok && std::isfinite(e.timeMeasured[i]) &&
             e.timeMeasured[i] != -1 && e.tofEdepGeV[i] > 0 &&
             std::isfinite(e.path[i]);
      }
      for (int l = 0; l < 9; ++l) {
        e.tkEdepGeV[l] = ked[l] * 1e-3;
        for (int c = 0; c < 3; ++c)
          e.tkPos[l][c] = kpos[l][c];
      }
      for (int i = 1; i < 4; ++i)
        ok = ok && (e.path[i] - e.path[i - 1] > 0) &&
             (e.path[i] - e.path[i - 1] < 300);
      if (!ok || !truthTimes(e))
        continue;
      int b = betaBin(e.betaCenter);
      if (b < 0)
        continue;
      ++valid;
      double it = 1 / e.betaCenter;
      residual["center_minus_checkpoint10_reference"][b]
          .add(it - 1 / e.betaCheckpoint10);
      auto add = [&](const std::string &n, const FitResult &r) {
        if (r.ok)
          residual[n][b].add(r.invBeta - it);
      };
      std::array<double, 4> rawTime{}, timeCorrection{};
      bool rawOk = true;
      for (int i = 0; i < 4; ++i) {
        rawTime[i] = -.5 * (rt[i][0] + rt[i][1]);
        timeCorrection[i] = e.timeMeasured[i] - rawTime[i];
        rawOk = rawOk && std::isfinite(rawTime[i]);
        physical["time_correction_S" + std::to_string(i + 1) + "_ns"][b]
            .add(timeCorrection[i]);
        physical["log10_geomAadc_S" + std::to_string(i + 1)][b]
            .add(.5 * (std::log10(std::max(1.0, (double)aadc[i][0])) +
                       std::log10(std::max(1.0, (double)aadc[i][1]))));
      }
      auto lm = fitLinear(e, e.timeMeasured);
      add("linear_measured_robust", lm);
      if (rawOk) {
        add("linear_raw_rtime", fitLinear(e, rawTime));
        Model rawOra = [&](double x, std::array<double, 4> &q) {
          return oracleTimes(e, x, q);
        };
        add("oracle_raw_rtime",
            fitModel(rawTime, rawOra,
                     1 / std::clamp((double)tb, .26, .999999)));
      }
      if (tbmc > 0 && std::isfinite(tbmc))
        residual["tof_betahmc"][b].add(1.0 / tbmc - it);
      if (tb > 0 && std::isfinite(tb))
        residual["tof_betah_branch"][b].add(1.0 / tb - it);
      add("linear_truth", fitLinear(e, e.timeTruth));
      Model cur = [&](double x, std::array<double, 4> &q) {
        return currentTimes(e, x, 1.9, 1, false, q);
      },
          common = [&](double x, std::array<double, 4> &q) {
            return currentTimes(e, x, commonScale, commonScale, false, q);
          },
          z43 = [&](double x, std::array<double, 4> &q) {
            return currentTimes(e, x, 4.3, 1, false, q);
          },
          mt = [&](double x, std::array<double, 4> &q) {
            return currentTimes(e, x, 1.9, 1, true, q);
          },
          ora = [&](double x,
                    std::array<double, 4> &q) { return oracleTimes(e, x, q); },
          leg = [&](double x, std::array<double, 4> &q) {
        return legacyTimes(e, x, 1.9, q);
      };
      double seed = 1 / std::clamp((double)tb, .26, .999999);
      auto cm = fitModel(e.timeMeasured, cur, seed);
      if (!cm.ok)
        ++cmFail;
      add("current_measured_robust", cm);
      auto commonMeasured = fitModel(e.timeMeasured, common, seed);
      add("common_scale_measured", commonMeasured);
      if (commonMeasured.ok)
        residual["common_scale_measured_checkpoint10"][b]
            .add(commonMeasured.invBeta - 1 / e.betaCheckpoint10);
      auto ct = fitModel(e.timeTruth, cur, it);
      if (!ct.ok)
        ++ctFail;
      add("current_truth", ct);
      add("common_scale_truth", fitModel(e.timeTruth, common, it));
      add("current_zeta4p3_truth", fitModel(e.timeTruth, z43, it));
      for (double ts : { 0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0 }) {
        Model sm = [&](double x, std::array<double, 4> &q) {
          return currentTimes(e, x, 1.9, ts, false, q);
        };
        add("current_truth_trackerScale_" + std::to_string((int)ts),
            fitModel(e.timeTruth, sm, it));
      }
      auto om = fitModel(e.timeMeasured, ora, seed);
      if (!om.ok)
        ++omFail;
      add("oracle_measured", om);
      auto ot = fitModel(e.timeTruth, ora, it);
      if (!ot.ok)
        ++otFail;
      add("oracle_truth", ot);
      add("current_zeta4p3_measured", fitModel(e.timeMeasured, z43, seed));
      add("current_mtof_measured", fitModel(e.timeMeasured, mt, seed));
      add("legacy_measured_robust", fitModel(e.timeMeasured, leg, seed));
      if (has) {
        if (out.nonlinear > 0)
          residual["production_nonlinear"][b].add(1 / out.nonlinear - it);
        if (out.nonlinear > 0 && out.mcBeta > 0)
          residual["production_nonlinear_checkpoint10"][b]
              .add(1 / out.nonlinear - 1 / out.mcBeta);
        if (out.legacy > 0)
          residual["production_legacy"][b].add(1 / out.legacy - it);
        if (out.linear > 0)
          residual["production_linear"][b].add(1 / out.linear - it);
        if (cm.ok && out.nonlinear > 0)
          residual["robust_minus_production_nonlinear"][b]
              .add(cm.invBeta - 1 / out.nonlinear);
      }
      double eu = truthEnergyAt(e, e.tofPos[1][2]),
             ed = truthEnergyAt(e, e.tofPos[2][2]), au = 0, ad = 0;
      for (int l = 1; l <= 4; ++l)
        au += e.tkEdepGeV[l];
      for (int l = 5; l <= 7; ++l)
        ad += e.tkEdepGeV[l];
      if (std::isfinite(eu + ed) && eu > e.centerTruthEnergy &&
          e.centerTruthEnergy > ed && au > 0 && ad > 0) {
        physical["truth_dE_upper_MeV"][b]
            .add((eu - e.centerTruthEnergy) * 1000);
        physical["truth_dE_lower_MeV"][b]
            .add((e.centerTruthEnergy - ed) * 1000);
        physical["tk_active_upper_MeV"][b].add(au * 1000);
        physical["tk_active_lower_MeV"][b].add(ad * 1000);
        physical["truth_over_active_upper"][b]
            .add((eu - e.centerTruthEnergy) / au);
        physical["truth_over_active_lower"][b]
            .add((e.centerTruthEnergy - ed) / ad);
        if (eventOut) {
          double logA2 = .5 * (std::log10(std::max(1.0, (double)aadc[1][0])) +
                               std::log10(std::max(1.0, (double)aadc[1][1]))),
                 logA3 = .5 * (std::log10(std::max(1.0, (double)aadc[2][0])) +
                               std::log10(std::max(1.0, (double)aadc[2][1])));
          eventOut << fi << ',' << Z << ',' << std::setprecision(10)
                   << e.betaCenter << ',' << barid[1] << ',' << barid[2] << ','
                   << ((e.timeMeasured[2] - e.timeMeasured[1]) -
                       (e.timeTruth[2] - e.timeTruth[1])) << ','
                   << ((rawTime[2] - rawTime[1]) -
                       (e.timeTruth[2] - e.timeTruth[1])) << ','
                   << (timeCorrection[2] - timeCorrection[1]) << ',' << logA2
                   << ',' << logA3 << ',' << (ted[1] / (Z * Z)) << ','
                   << (ted[2] / (Z * Z)) << ','
                   << ((eu - e.centerTruthEnergy) * 1000) << ','
                   << ((e.centerTruthEnergy - ed) * 1000) << ',' << (au * 1000)
                   << ',' << (ad * 1000) << '\n';
        }
      }
      std::array<double, 4> qc{}, qo{};
      currentTimes(e, it, 1.9, 1, false, qc);
      std::array<double, 4> qcommon{};
      currentTimes(e, it, commonScale, commonScale, false, qcommon);
      oracleTimes(e, it, qo);
      const char *iv[3] = { "S1S2", "S2S3", "S3S4" };
      for (int i = 0; i < 3; ++i) {
        timing[std::string("measured_minus_truth_") + iv[i]][b]
            .add((e.timeMeasured[i + 1] - e.timeMeasured[i]) -
                 (e.timeTruth[i + 1] - e.timeTruth[i]));
        timing[std::string("raw_rtime_minus_truth_") + iv[i]][b]
            .add((rawTime[i + 1] - rawTime[i]) -
                 (e.timeTruth[i + 1] - e.timeTruth[i]));
        timing[std::string("calibration_correction_") + iv[i]][b]
            .add(timeCorrection[i + 1] - timeCorrection[i]);
        timing[std::string("current_minus_truth_") + iv[i]][b]
            .add((qc[i + 1] - qc[i]) - (e.timeTruth[i + 1] - e.timeTruth[i]));
        timing[std::string("common_scale_minus_truth_") + iv[i]][b]
            .add((qcommon[i + 1] - qcommon[i]) -
                 (e.timeTruth[i + 1] - e.timeTruth[i]));
        timing[std::string("oracle_minus_truth_") + iv[i]][b]
            .add((qo[i + 1] - qo[i]) - (e.timeTruth[i + 1] - e.timeTruth[i]));
      }
    }
    std::cerr << "file " << fi + 1 << "/" << files.size() << " valid=" << valid
              << "\n";
  }
  {
    std::ofstream o(resPath);
    o << "method,beta_center,n,mean,median,sigma68,rms,q05,q95,core_mu,core_sigma\n"
      << std::setprecision(10);
    for (auto &kv : residual)
      for (int b = 0; b < kBins; ++b) {
        auto s = summarize(kv.second[b]);
        if (!s.n)
          continue;
        double bc = kBetaMin + (b + .5) * (kBetaMax - kBetaMin) / kBins;
        o << kv.first << ',' << bc << ',' << s.n << ',' << s.mean << ','
          << s.median << ',' << s.sigma68 << ',' << s.rms << ',' << s.q05 << ','
          << s.q95 << ',' << s.coreMu << ',' << s.coreSigma << '\n';
      }
  }
  {
    std::ofstream o(timPath);
    o << "source_interval,beta_center,n,mean_ns,median_ns,sigma68_ns,rms_ns,"
         "q05_ns,q95_ns\n" << std::setprecision(10);
    for (auto &kv : timing)
      for (int b = 0; b < kBins; ++b) {
        auto s = summarize(kv.second[b]);
        if (!s.n)
          continue;
        double bc = kBetaMin + (b + .5) * (kBetaMax - kBetaMin) / kBins;
        o << kv.first << ',' << bc << ',' << s.n << ',' << s.mean << ','
          << s.median << ',' << s.sigma68 << ',' << s.rms << ',' << s.q05 << ','
          << s.q95 << '\n';
      }
  }
  {
    std::ofstream o(physPath);
    o << "quantity,beta_center,n,mean,median,sigma68,rms,q05,q95\n"
      << std::setprecision(10);
    for (auto &kv : physical)
      for (int b = 0; b < kBins; ++b) {
        auto s = summarize(kv.second[b]);
        if (!s.n)
          continue;
        double bc = kBetaMin + (b + .5) * (kBetaMax - kBetaMin) / kBins;
        o << kv.first << ',' << bc << ',' << s.n << ',' << s.mean << ','
          << s.median << ',' << s.sigma68 << ',' << s.rms << ',' << s.q05 << ','
          << s.q95 << '\n';
      }
  }
  {
    std::ofstream o(sumPath);
    o << "files=" << files.size() << "\nraw=" << raw
      << "\nselected=" << selected << "\nvalid=" << valid
      << "\npaired=" << paired << "\npair_mismatch=" << mismatch
      << "\ncurrent_measured_fail=" << cmFail
      << "\ncurrent_truth_fail=" << ctFail
      << "\noracle_measured_fail=" << omFail << "\noracle_truth_fail=" << otFail
      << "\n";
  }
  std::cerr << "done valid=" << valid << "\n";
}
