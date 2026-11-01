class KernelMaterialModel {
public:
  void add(const Event &event) {
    const int bin = learnBin(event.betaCenter);
    if (bin < 0)
      return;
    const double centerInvBeta = 1.0 / event.betaCenter;
    for (int interval = 0; interval < kIntervals; ++interval) {
      const double length = event.path[interval + 1] - event.path[interval];
      if (!(length > 0 && length < 300))
        continue;
      const double averageInvBeta =
          (event.timeTruth[interval + 1] - event.timeTruth[interval]) * kC /
          length;
      const double correction = averageInvBeta - centerInvBeta;
      if (std::isfinite(correction)) {
        samples_[interval][bin].add(correction);
        betaSamples_[interval][bin].add(event.betaCenter);
      }
    }
  }

  void finalize() {
    for (int interval = 0; interval < kIntervals; ++interval) {
      for (int bin = 0; bin < kLearnBins; ++bin) {
        counts_[interval][bin] = samples_[interval][bin].v.size();
        knotBeta_[interval][bin] = learnCenter(bin);
        if (counts_[interval][bin] >= 20) {
          correction_[interval][bin] = median(samples_[interval][bin].v);
          knotBeta_[interval][bin] = median(betaSamples_[interval][bin].v);
        }
      }
      fillMissing(correction_[interval]);
    }
  }

  bool times(const Event &event, double invBeta,
             std::array<double, 4> &times) const {
    if (!(invBeta > 1.0 && invBeta < 5.0))
      return false;
    const double beta = 1.0 / invBeta;
    times.fill(0.0);
    for (int interval = 0; interval < kIntervals; ++interval) {
      const double length = event.path[interval + 1] - event.path[interval];
      const double averageInvBeta =
          invBeta + lookup(correction_[interval], knotBeta_[interval], beta);
      if (!(length > 0 && length < 300) || !(averageInvBeta > 0) ||
          !std::isfinite(averageInvBeta))
        return false;
      times[interval + 1] = times[interval] + length * averageInvBeta / kC;
    }
    return true;
  }

  void write(const std::string &path) const {
    std::ofstream output(path);
    output << "interval,beta_knot,n,delta_average_invbeta\n";
    output << std::setprecision(10);
    for (int interval = 0; interval < kIntervals; ++interval)
      for (int bin = 0; bin < kLearnBins; ++bin)
        output << interval << ',' << knotBeta_[interval][bin] << ','
               << counts_[interval][bin] << ',' << correction_[interval][bin]
               << '\n';
  }

private:
  std::array<std::array<Stat, kLearnBins>, kIntervals> samples_{};
  std::array<std::array<Stat, kLearnBins>, kIntervals> betaSamples_{};
  std::array<std::array<long long, kLearnBins>, kIntervals> counts_{};
  std::array<std::array<double, kLearnBins>, kIntervals> correction_{};
  std::array<std::array<double, kLearnBins>, kIntervals> knotBeta_{};

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

  static double lookup(const std::array<double, kLearnBins> &table,
                       const std::array<double, kLearnBins> &knots,
                       double beta) {
    if (beta <= knots.front())
      return table.front();
    if (beta >= knots.back())
      return table.back();
    int upper = 1;
    while (upper < kLearnBins && knots[upper] < beta)
      ++upper;
    upper = std::min(upper, kLearnBins - 1);
    const int lower = upper - 1;
    const double width = knots[upper] - knots[lower];
    if (!(width > 0))
      return table[lower];
    const double fraction = std::clamp((beta - knots[lower]) / width, 0.0, 1.0);
    return table[lower] * (1.0 - fraction) + table[upper] * fraction;
  }
};
