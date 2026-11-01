class ConstrainedMaterialModel {
public:
  void add(const Event &event) {
    const int bin = learnBin(event.betaCenter);
    if (bin < 0)
      return;
    SegmentGeometry g;
    if (!segmentGeometry(event, g))
      return;
    std::array<double, 7> upper{};
    std::array<double, 6> lower{};
    for (int index = 0; index < 7; ++index)
      upper[index] = truthEnergyAt(event, g.upperZ[index]);
    for (int index = 0; index < 6; ++index)
      lower[index] = truthEnergyAt(event, g.lowerZ[index]);
    if (!std::all_of(upper.begin(), upper.end(),
                     [](double value) { return std::isfinite(value); }) ||
        !std::all_of(lower.begin(), lower.end(),
                     [](double value) { return std::isfinite(value); }))
      return;

    const double upperLength =
        std::accumulate(g.length.begin(), g.length.begin() + 5, 0.0);
    const double lowerLength =
        std::accumulate(g.length.begin() + 6, g.length.begin() + 10, 0.0);
    const double topLoss = upper[6] - upper[5];
    const double upperLoss = upper[5] - upper[0];
    const double lowerLoss = lower[0] - lower[4];
    const double bottomLoss = lower[4] - lower[5];
    const std::array<double, 4> rate{ topLoss / g.length[5],
                                      upperLoss / upperLength,
                                      lowerLoss / lowerLength,
                                      bottomLoss / g.length[10] };
    for (int region = 0; region < 4; ++region)
      if (std::isfinite(rate[region]))
        rateSamples_[region][bin].add(rate[region]);

    if (upperLoss > 0) {
      for (int node = 1; node <= 4; ++node) {
        const double fraction = (upper[node] - upper[0]) / upperLoss;
        if (std::isfinite(fraction))
          upperFractionSamples_[node - 1][bin].add(fraction);
      }
    }
    if (lowerLoss > 0) {
      for (int node = 1; node <= 3; ++node) {
        const double fraction = (lower[0] - lower[node]) / lowerLoss;
        if (std::isfinite(fraction))
          lowerFractionSamples_[node - 1][bin].add(fraction);
      }
    }
  }

  void finalize() {
    for (int region = 0; region < 4; ++region)
      buildTable(rateSamples_[region], rate_[region], rateCounts_[region],
                 true);
    for (int node = 0; node < 4; ++node)
      buildTable(upperFractionSamples_[node], upperFraction_[node],
                 upperFractionCounts_[node], false);
    for (int node = 0; node < 3; ++node)
      buildTable(lowerFractionSamples_[node], lowerFraction_[node],
                 lowerFractionCounts_[node], false);

    for (int bin = 0; bin < kLearnBins; ++bin) {
      double previous = 0;
      for (int node = 0; node < 4; ++node) {
        upperFraction_[node][bin] =
            std::clamp(std::max(previous, upperFraction_[node][bin]), 0.0, 1.0);
        previous = upperFraction_[node][bin];
      }
      previous = 0;
      for (int node = 0; node < 3; ++node) {
        lowerFraction_[node][bin] =
            std::clamp(std::max(previous, lowerFraction_[node][bin]), 0.0, 1.0);
        previous = lowerFraction_[node][bin];
      }
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

    const double upperLength =
        std::accumulate(g.length.begin(), g.length.begin() + 5, 0.0);
    const double lowerLength =
        std::accumulate(g.length.begin() + 6, g.length.begin() + 10, 0.0);
    const double topLoss = lookup(rate_[0], beta) * g.length[5];
    const double upperLoss = lookup(rate_[1], beta) * upperLength;
    const double lowerLoss = lookup(rate_[2], beta) * lowerLength;
    const double bottomLoss = lookup(rate_[3], beta) * g.length[10];
    if (!(topLoss >= 0 && upperLoss >= 0 && lowerLoss >= 0 && bottomLoss >= 0))
      return false;

    std::array<double, 7> upperEnergy{};
    upperEnergy[0] = center;
    for (int node = 1; node <= 4; ++node)
      upperEnergy[node] =
          center + lookup(upperFraction_[node - 1], beta) * upperLoss;
    upperEnergy[5] = center + upperLoss;
    upperEnergy[6] = upperEnergy[5] + topLoss;

    std::array<double, 6> lowerEnergy{};
    lowerEnergy[0] = center;
    for (int node = 1; node <= 3; ++node)
      lowerEnergy[node] =
          center - lookup(lowerFraction_[node - 1], beta) * lowerLoss;
    lowerEnergy[4] = center - lowerLoss;
    lowerEnergy[5] = lowerEnergy[4] - bottomLoss;
    if (!(lowerEnergy[5] > mass))
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

    double upperTime = 0;
    for (int segment = 0; segment < 5; ++segment) {
      const double dt = integrate(g.length[segment], upperEnergy[segment],
                                  upperEnergy[segment + 1]);
      if (!std::isfinite(dt))
        return false;
      upperTime += dt;
    }
    const double topTime =
        integrate(g.length[5], upperEnergy[5], upperEnergy[6]);
    double lowerTime = 0;
    for (int segment = 0; segment < 4; ++segment) {
      const double dt = integrate(g.length[6 + segment], lowerEnergy[segment],
                                  lowerEnergy[segment + 1]);
      if (!std::isfinite(dt))
        return false;
      lowerTime += dt;
    }
    const double bottomTime =
        integrate(g.length[10], lowerEnergy[4], lowerEnergy[5]);
    if (!std::isfinite(topTime + bottomTime))
      return false;
    times[1] = -upperTime;
    times[0] = times[1] - topTime;
    times[2] = lowerTime;
    times[3] = times[2] + bottomTime;
    return true;
  }

  void write(const std::string &path) const {
    std::ofstream output(path);
    output << "kind,index,beta_center,n,value\n";
    output << std::setprecision(10);
    for (int region = 0; region < 4; ++region)
      for (int bin = 0; bin < kLearnBins; ++bin)
        output << "rate," << region << ',' << learnCenter(bin) << ','
               << rateCounts_[region][bin] << ',' << rate_[region][bin] << '\n';
    for (int node = 0; node < 4; ++node)
      for (int bin = 0; bin < kLearnBins; ++bin)
        output << "upper_fraction," << node << ',' << learnCenter(bin) << ','
               << upperFractionCounts_[node][bin] << ','
               << upperFraction_[node][bin] << '\n';
    for (int node = 0; node < 3; ++node)
      for (int bin = 0; bin < kLearnBins; ++bin)
        output << "lower_fraction," << node << ',' << learnCenter(bin) << ','
               << lowerFractionCounts_[node][bin] << ','
               << lowerFraction_[node][bin] << '\n';
  }

private:
  using SampleTable = std::array<Stat, kLearnBins>;
  using ValueTable = std::array<double, kLearnBins>;
  using CountTable = std::array<long long, kLearnBins>;
  std::array<SampleTable, 4> rateSamples_{};
  std::array<ValueTable, 4> rate_{};
  std::array<CountTable, 4> rateCounts_{};
  std::array<SampleTable, 4> upperFractionSamples_{};
  std::array<ValueTable, 4> upperFraction_{};
  std::array<CountTable, 4> upperFractionCounts_{};
  std::array<SampleTable, 3> lowerFractionSamples_{};
  std::array<ValueTable, 3> lowerFraction_{};
  std::array<CountTable, 3> lowerFractionCounts_{};

  static void buildTable(const SampleTable &samples, ValueTable &values,
                         CountTable &counts, bool nonnegative) {
    for (int bin = 0; bin < kLearnBins; ++bin) {
      counts[bin] = samples[bin].v.size();
      if (counts[bin] >= 20) {
        values[bin] = median(samples[bin].v);
        if (nonnegative)
          values[bin] = std::max(0.0, values[bin]);
      }
    }
    fillMissing(values);
    smooth(values);
  }

  static void fillMissing(ValueTable &table) {
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

  static void smooth(ValueTable &table) {
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

  static double lookup(const ValueTable &table, double beta) {
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
