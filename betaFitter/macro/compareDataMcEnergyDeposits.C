#include <TChain.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct Values {
    std::vector<double> beta;
    std::vector<double> tof[4];
    std::vector<double> tracker[7];
    std::vector<double> tofSum;
    std::vector<double> trackerSum;
};

double quantile(std::vector<double> values, double probability)
{
    if (values.empty())
        return 0;
    std::sort(values.begin(), values.end());
    const double position = probability * (values.size() - 1);
    const size_t low = static_cast<size_t>(position);
    const size_t high = std::min(low + 1, values.size() - 1);
    const double fraction = position - low;
    return values[low] * (1.0 - fraction) + values[high] * fraction;
}

void writeMetric(std::ofstream &output, const char *sample, int charge,
                 const char *metric, const std::vector<double> &values)
{
    output << sample << '\t' << charge << '\t' << metric << '\t'
           << values.size() << '\t' << quantile(values, 0.16) << '\t'
           << quantile(values, 0.50) << '\t' << quantile(values, 0.84)
           << '\n';
}

Values readSample(const char *pattern, int requestedCharge, bool data,
                  double betaLow, double betaHigh)
{
    TChain chain("amstreea");
    chain.Add(pattern);
    float tofBeta = 0;
    float tofEdep[4]{};
    float trackerEdep[9]{};
    float trackerQ[2]{};
    float trackerInnerQ[2][3]{};
    chain.SetBranchAddress("tof_betah", &tofBeta);
    chain.SetBranchAddress("tof_edep", tofEdep);
    chain.SetBranchAddress("tk_edep", trackerEdep);
    chain.SetBranchAddress("tk_q", trackerQ);
    chain.SetBranchAddress("tk_qin", trackerInnerQ);

    Values result;
    const double chargeSquared = requestedCharge * requestedCharge;
    // Deterministic prescale is sufficient for response-scale comparisons.
    constexpr Long64_t mcStride = 5;
    for (Long64_t entry = 0; entry < chain.GetEntries(); entry += mcStride) {
        chain.GetEntry(entry);
        const int charge = static_cast<int>(
            (trackerInnerQ[0][2] < 2.5 ? trackerQ[1]
                                      : trackerInnerQ[0][2]) +
            0.5);
        if (charge != requestedCharge || !std::isfinite(tofBeta) ||
            tofBeta <= betaLow || tofBeta >= betaHigh)
            continue;

        bool valid = true;
        double tofSum = 0;
        double trackerSum = 0;
        for (int station = 0; station < 4; ++station) {
            valid = valid && std::isfinite(tofEdep[station]) &&
                    tofEdep[station] > 0;
            tofSum += tofEdep[station];
        }
        for (int layer = 1; layer <= 7; ++layer) {
            valid = valid && std::isfinite(trackerEdep[layer]) &&
                    trackerEdep[layer] > 0;
            trackerSum += trackerEdep[layer];
        }
        if (!valid)
            continue;

        result.beta.push_back(tofBeta);
        result.tofSum.push_back(tofSum / chargeSquared);
        result.trackerSum.push_back(trackerSum / chargeSquared);
        for (int station = 0; station < 4; ++station)
            result.tof[station].push_back(tofEdep[station] / chargeSquared);
        for (int layer = 1; layer <= 7; ++layer)
            result.tracker[layer - 1].push_back(
                trackerEdep[layer] / chargeSquared);
    }
    std::cout << (data ? "data" : "mc") << " Z=" << requestedCharge
              << " input=" << chain.GetEntries()
              << " selected=" << result.beta.size() << std::endl;
    return result;
}

std::array<Values, 3> readDataSamples(const char *pattern)
{
    TChain chain("amstreea");
    chain.Add(pattern);
    float tofBeta = 0;
    float tofEdep[4]{};
    float trackerEdep[9]{};
    float trackerQ[2]{};
    float trackerInnerQ[2][3]{};
    chain.SetBranchAddress("tof_betah", &tofBeta);
    chain.SetBranchAddress("tof_edep", tofEdep);
    chain.SetBranchAddress("tk_edep", trackerEdep);
    chain.SetBranchAddress("tk_q", trackerQ);
    chain.SetBranchAddress("tk_qin", trackerInnerQ);

    std::array<Values, 3> results;
    constexpr Long64_t dataStride = 1;
    for (Long64_t entry = 0; entry < chain.GetEntries(); entry += dataStride) {
        chain.GetEntry(entry);
        const int charge = static_cast<int>(
            (trackerInnerQ[0][2] < 2.5 ? trackerQ[1] : trackerInnerQ[0][2]) +
            0.5);
        const int index = charge == 2 ? 0 : (charge == 6 ? 1 : (charge == 8 ? 2 : -1));
        if (index < 0 || !std::isfinite(tofBeta))
            continue;
        const double betaLow = charge == 2 ? 0.35 : (charge == 6 ? 0.40 : 0.45);
        const double betaHigh = charge == 2 ? 0.50 : 0.60;
        if (tofBeta <= betaLow || tofBeta >= betaHigh)
            continue;

        bool valid = true;
        double tofSum = 0;
        double trackerSum = 0;
        for (int station = 0; station < 4; ++station) {
            valid = valid && std::isfinite(tofEdep[station]) &&
                    tofEdep[station] > 0;
            tofSum += tofEdep[station];
        }
        for (int layer = 1; layer <= 7; ++layer) {
            valid = valid && std::isfinite(trackerEdep[layer]) &&
                    trackerEdep[layer] > 0;
            trackerSum += trackerEdep[layer];
        }
        if (!valid)
            continue;

        Values &result = results[index];
        const double chargeSquared = charge * charge;
        result.beta.push_back(tofBeta);
        result.tofSum.push_back(tofSum / chargeSquared);
        result.trackerSum.push_back(trackerSum / chargeSquared);
        for (int station = 0; station < 4; ++station)
            result.tof[station].push_back(tofEdep[station] / chargeSquared);
        for (int layer = 1; layer <= 7; ++layer)
            result.tracker[layer - 1].push_back(
                trackerEdep[layer] / chargeSquared);
    }
    std::cout << "data input=" << chain.GetEntries()
              << " selected Z2/Z6/Z8=" << results[0].beta.size() << "/"
              << results[1].beta.size() << "/" << results[2].beta.size()
              << std::endl;
    return results;
}

void writeSample(std::ofstream &output, const char *sample, int charge,
                 const Values &values)
{
    writeMetric(output, sample, charge, "tof_beta", values.beta);
    writeMetric(output, sample, charge, "tof_edep_sum_over_Z2",
                values.tofSum);
    writeMetric(output, sample, charge, "tracker_edep_sum_over_Z2",
                values.trackerSum);
    for (int station = 0; station < 4; ++station)
        writeMetric(output, sample, charge,
                    Form("tof_edep_S%d_over_Z2", station + 1),
                    values.tof[station]);
    for (int layer = 0; layer < 7; ++layer)
        writeMetric(output, sample, charge,
                    Form("tracker_edep_L%d_over_Z2", layer + 2),
                    values.tracker[layer]);
}

} // namespace

void compareDataMcEnergyDeposits(
    const char *dataPattern =
        "/lustre02/user/yachen/public/"
        "amsd69n_B1130_Data_TofCalibN2/*.root",
    const char *mcRoot =
        "/lustre02/user/yachen/public/hxwu/"
        "betaFitter_rdst_amsd69n_B1308_20260731",
    const char *outputName =
        "macro/results/data_amsd69n_center_common_20260805_z268_lowR/"
        "data_mc_edep_summary.tsv")
{
    std::ofstream output(outputName);
    output << std::setprecision(10);
    output << "sample\tZ\tmetric\tcount\tq16\tmedian\tq84\n";
    const std::array<Values, 3> dataSamples = readDataSamples(dataPattern);
    int sampleIndex = 0;
    for (int charge : {2, 6, 8}) {
        const double betaLow = charge == 2 ? 0.35 : (charge == 6 ? 0.40 : 0.45);
        const double betaHigh = charge == 2 ? 0.50 : 0.60;
        const char *directory =
            charge == 2 ? "He4" : (charge == 6 ? "C12" : "O16");
        const std::string mcPattern =
            std::string(mcRoot) + "/" + directory + "/*.root";
        const Values mc =
            readSample(mcPattern.c_str(), charge, false, betaLow, betaHigh);
        writeSample(output, "data_period0", charge,
                    dataSamples[sampleIndex]);
        writeSample(output, "mc_prescale5", charge, mc);
        ++sampleIndex;
    }
    std::cout << "Wrote " << outputName << std::endl;
}
