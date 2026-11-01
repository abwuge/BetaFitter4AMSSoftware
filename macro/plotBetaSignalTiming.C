#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TChain.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{
constexpr double kC = 29.9792458; // cm/ns
constexpr double kIqrToSigma = 0.741301109252801;
constexpr double kBetaMin = 0.30;
constexpr double kBetaMax = 1.00;
constexpr int kNBetaBins = 28;
constexpr int kMinEntries = 100;
constexpr double kStoredTofTimeErrorNs = 0.154481;

struct BinData
{
    std::vector<double> rawInvBetaResidual;
    std::vector<double> betaHInvBetaResidual;
    std::vector<double> timingEquivalent;
    std::vector<double> meanEdep;
    std::vector<double> meanAdc;
};

struct ChargeData
{
    int charge = 0;
    Long64_t selected = 0;
    std::array<BinData, kNBetaBins> bins;
};

struct Series
{
    std::vector<double> x;
    std::vector<double> y;
};

double quantileSorted(const std::vector<double> &values, double probability)
{
    if (values.empty())
        return std::numeric_limits<double>::quiet_NaN();
    const double position = probability * (values.size() - 1);
    const size_t lower = static_cast<size_t>(std::floor(position));
    const size_t upper = static_cast<size_t>(std::ceil(position));
    if (lower == upper)
        return values[lower];
    const double fraction = position - lower;
    return values[lower] * (1.0 - fraction) + values[upper] * fraction;
}

double median(std::vector<double> values)
{
    if (values.empty())
        return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    return quantileSorted(values, 0.5);
}

double robustSigma(std::vector<double> values)
{
    values.erase(std::remove_if(values.begin(), values.end(),
                                [](double value) { return !std::isfinite(value); }),
                 values.end());
    if (values.size() < 8)
        return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    return kIqrToSigma *
           (quantileSorted(values, 0.75) - quantileSorted(values, 0.25));
}

std::vector<std::string> validationFiles(const char *sourceFile,
                                         int charge,
                                         int maxFiles)
{
    std::ifstream input(sourceFile);
    std::vector<std::string> files;
    if (!input)
    {
        std::cerr << "Cannot open source list: " << sourceFile << std::endl;
        return files;
    }

    std::string line;
    while (std::getline(input, line))
    {
        const std::string::size_type first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#')
            continue;
        std::istringstream fields(line);
        std::string role;
        int z = 0;
        std::string path;
        fields >> role >> z >> path;
        if (!fields || role != "validation" || z != charge)
            continue;
        files.push_back(path);
        if (maxFiles > 0 && static_cast<int>(files.size()) >= maxFiles)
            break;
    }
    return files;
}

ChargeData analyzeCharge(const char *sourceFile, int charge, int maxFiles)
{
    ChargeData result;
    result.charge = charge;

    const auto files = validationFiles(sourceFile, charge, maxFiles);
    if (files.empty())
        return result;

    TChain tree("amstreea");
    for (const auto &file : files)
        tree.Add(file.c_str());

    float mch = 0;
    float mmass = 0;
    float mevmom1[21] = {};
    float tofBetah = 0;
    float tofTl[4] = {};
    float tofLeng[4] = {};
    float tofEdep[4] = {};
    float tofAadc[4][2] = {};
    float tkQ[2] = {};
    float tkQin[2][3] = {};
    float tkPos[9][3] = {};
    float tkEdep[9] = {};

    const bool hasAdc = tree.GetBranch("tof_Aadc") != nullptr;
    tree.SetBranchStatus("*", 0);
    const char *branches[] = {"mch", "mmass", "mevmom1", "tof_betah",
                              "tof_tl", "tof_leng", "tof_edep", "tk_q",
                              "tk_qin", "tk_pos", "tk_edep"};
    for (const char *branch : branches)
        tree.SetBranchStatus(branch, 1);
    if (hasAdc)
        tree.SetBranchStatus("tof_Aadc", 1);

    tree.SetBranchAddress("mch", &mch);
    tree.SetBranchAddress("mmass", &mmass);
    tree.SetBranchAddress("mevmom1", mevmom1);
    tree.SetBranchAddress("tof_betah", &tofBetah);
    tree.SetBranchAddress("tof_tl", tofTl);
    tree.SetBranchAddress("tof_leng", tofLeng);
    tree.SetBranchAddress("tof_edep", tofEdep);
    if (hasAdc)
        tree.SetBranchAddress("tof_Aadc", tofAadc);
    tree.SetBranchAddress("tk_q", tkQ);
    tree.SetBranchAddress("tk_qin", tkQin);
    tree.SetBranchAddress("tk_pos", tkPos);
    tree.SetBranchAddress("tk_edep", tkEdep);

    const double binWidth = (kBetaMax - kBetaMin) / kNBetaBins;
    const Long64_t entries = tree.GetEntries();
    std::cout << "Z=" << charge << ": " << files.size()
              << " validation files, " << entries << " raw entries" << std::endl;

    for (Long64_t entry = 0; entry < entries; ++entry)
    {
        tree.GetEntry(entry);

        if (!(mevmom1[10] > 0) || !std::isfinite(mevmom1[10]) ||
            !(mmass > 0) || !std::isfinite(mmass))
            continue;

        const int mcCharge = static_cast<int>(std::floor(mch + 0.5));
        const double measuredChargeValue = tkQin[0][2] < 2.5 ? tkQ[1] : tkQin[0][2];
        const int measuredCharge = static_cast<int>(std::floor(measuredChargeValue + 0.5));
        if (mcCharge != charge || measuredCharge != charge)
            continue;

        bool trackerComplete = true;
        for (int layer = 1; layer <= 7; ++layer)
        {
            trackerComplete = trackerComplete && tkEdep[layer] > 0 &&
                              std::isfinite(tkEdep[layer]);
            for (int coordinate = 0; coordinate < 3; ++coordinate)
                trackerComplete = trackerComplete &&
                                  std::isfinite(tkPos[layer][coordinate]);
        }
        if (!trackerComplete || !(tkPos[4][2] > 0) || !(tkPos[5][2] < 0))
            continue;

        const double momentum = mevmom1[10];
        const double beta = momentum / std::sqrt(momentum * momentum + mmass * mmass);
        if (!(beta >= kBetaMin && beta < kBetaMax) || !std::isfinite(beta))
            continue;

        bool fourHits = true;
        for (int layer = 0; layer < 4; ++layer)
            fourHits = fourHits && tofTl[layer] != -1 && std::isfinite(tofTl[layer]) &&
                       std::isfinite(tofLeng[layer]) && tofEdep[layer] > 0 &&
                       std::isfinite(tofEdep[layer]);
        if (!fourHits)
            continue;

        double meanLength = 0;
        double meanTime = 0;
        for (int layer = 0; layer < 4; ++layer)
        {
            meanLength += tofLeng[layer];
            meanTime += tofTl[layer];
        }
        meanLength /= 4.0;
        meanTime /= 4.0;

        double sxx = 0;
        double sxy = 0;
        for (int layer = 0; layer < 4; ++layer)
        {
            const double dx = tofLeng[layer] - meanLength;
            const double dy = tofTl[layer] - meanTime;
            sxx += dx * dx;
            sxy += dx * dy;
        }
        if (!(sxx > 0))
            continue;

        const double rawInvBeta = kC * sxy / sxx;
        const double trueInvBeta = 1.0 / beta;
        const double residual = rawInvBeta - trueInvBeta;
        const double timingEquivalent = residual * std::sqrt(sxx) / kC;

        double edep = 0;
        for (int layer = 0; layer < 4; ++layer)
            edep += tofEdep[layer];
        edep /= 4.0;

        double adc = std::numeric_limits<double>::quiet_NaN();
        if (hasAdc)
        {
            double sum = 0;
            int count = 0;
            for (int layer = 0; layer < 4; ++layer)
                for (int side = 0; side < 2; ++side)
                    if (tofAadc[layer][side] > 0 && std::isfinite(tofAadc[layer][side]))
                    {
                        sum += tofAadc[layer][side];
                        ++count;
                    }
            if (count == 8)
                adc = sum / count;
        }

        int bin = static_cast<int>((beta - kBetaMin) / binWidth);
        bin = std::max(0, std::min(kNBetaBins - 1, bin));
        BinData &data = result.bins[bin];
        data.rawInvBetaResidual.push_back(residual);
        data.timingEquivalent.push_back(timingEquivalent);
        data.meanEdep.push_back(edep);
        if (std::isfinite(adc))
            data.meanAdc.push_back(adc);
        if (tofBetah > 0 && std::isfinite(tofBetah))
            data.betaHInvBetaResidual.push_back(1.0 / tofBetah - trueInvBeta);
        ++result.selected;
    }

    std::cout << "Z=" << charge << ": selected four-hit N=" << result.selected << std::endl;
    return result;
}

Series sigmaSeries(const ChargeData &data, bool betaH)
{
    Series result;
    const double width = (kBetaMax - kBetaMin) / kNBetaBins;
    for (int bin = 0; bin < kNBetaBins; ++bin)
    {
        const auto &values = betaH ? data.bins[bin].betaHInvBetaResidual
                                   : data.bins[bin].rawInvBetaResidual;
        if (static_cast<int>(values.size()) < kMinEntries)
            continue;
        const double sigma = robustSigma(values);
        if (!std::isfinite(sigma))
            continue;
        result.x.push_back(kBetaMin + (bin + 0.5) * width);
        result.y.push_back(sigma);
    }
    return result;
}

Series timingSeries(const ChargeData &data)
{
    Series result;
    const double width = (kBetaMax - kBetaMin) / kNBetaBins;
    for (int bin = 0; bin < kNBetaBins; ++bin)
    {
        if (static_cast<int>(data.bins[bin].timingEquivalent.size()) < kMinEntries)
            continue;
        const double sigma = robustSigma(data.bins[bin].timingEquivalent);
        if (!std::isfinite(sigma))
            continue;
        result.x.push_back(kBetaMin + (bin + 0.5) * width);
        result.y.push_back(sigma);
    }
    return result;
}

Series medianSeries(const ChargeData &data, bool adc)
{
    Series result;
    const double width = (kBetaMax - kBetaMin) / kNBetaBins;
    for (int bin = 0; bin < kNBetaBins; ++bin)
    {
        const auto &values = adc ? data.bins[bin].meanAdc : data.bins[bin].meanEdep;
        if (static_cast<int>(values.size()) < kMinEntries)
            continue;
        const double value = median(values);
        if (!std::isfinite(value))
            continue;
        result.x.push_back(kBetaMin + (bin + 0.5) * width);
        result.y.push_back(value);
    }
    return result;
}

Series timingVsSignal(const ChargeData &data, bool adc)
{
    Series result;
    for (int bin = 0; bin < kNBetaBins; ++bin)
    {
        const auto &signal = adc ? data.bins[bin].meanAdc : data.bins[bin].meanEdep;
        if (static_cast<int>(data.bins[bin].timingEquivalent.size()) < kMinEntries ||
            static_cast<int>(signal.size()) < kMinEntries)
            continue;
        const double x = median(signal);
        const double y = robustSigma(data.bins[bin].timingEquivalent);
        if (!std::isfinite(x) || !std::isfinite(y))
            continue;
        result.x.push_back(x);
        result.y.push_back(y);
    }
    return result;
}

TGraph *makeGraph(const Series &series, int color, int marker, const char *name)
{
    auto *graph = new TGraph(series.x.size(), series.x.data(), series.y.data());
    graph->SetName(name);
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(1.05);
    graph->SetLineWidth(3);
    return graph;
}

void configurePad()
{
    gPad->SetLeftMargin(0.145);
    gPad->SetRightMargin(0.045);
    gPad->SetBottomMargin(0.14);
    gPad->SetTopMargin(0.08);
    gPad->SetGridx();
    gPad->SetGridy();
}

void styleAxis(TH1 *frame, const char *xTitle, const char *yTitle)
{
    frame->GetXaxis()->SetTitle(xTitle);
    frame->GetYaxis()->SetTitle(yTitle);
    frame->GetXaxis()->SetTitleSize(0.050);
    frame->GetYaxis()->SetTitleSize(0.050);
    frame->GetXaxis()->SetLabelSize(0.041);
    frame->GetYaxis()->SetLabelSize(0.041);
    frame->GetYaxis()->SetTitleOffset(1.35);
}

std::pair<double, double> rangeFor(const std::vector<Series> &series,
                                   double floorValue = 0.0)
{
    double minimum = std::numeric_limits<double>::max();
    double maximum = -std::numeric_limits<double>::max();
    for (const auto &item : series)
        for (double value : item.y)
        {
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
        }
    if (!std::isfinite(minimum) || !std::isfinite(maximum))
        return std::make_pair(0.0, 1.0);
    const double span = std::max(1e-6, maximum - minimum);
    return std::make_pair(std::min(floorValue, minimum - 0.10 * span),
                          maximum + 0.18 * span);
}
} // namespace

void plotBetaSignalTiming(
    const char *sourceFile = "config/hybrid_material_fine.sources",
    const char *outputDir = "macro/results/beta_signal_timing",
    int maxFilesPerCharge = -1)
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);
    gSystem->mkdir(outputDir, true);

    const std::array<int, 4> charges = {1, 2, 6, 8};
    const std::array<int, 4> colors = {
        TColor::GetColor("#A63D40"),
        TColor::GetColor("#D2691E"),
        TColor::GetColor("#C18F00"),
        TColor::GetColor("#7F5539")};
    const std::array<int, 4> markers = {20, 21, 22, 33};

    std::array<ChargeData, 4> data;
    for (size_t i = 0; i < charges.size(); ++i)
        data[i] = analyzeCharge(sourceFile, charges[i], maxFilesPerCharge);

    // Numerical summary first, so the data remain useful even if drawing fails.
    std::ofstream summary(Form("%s/beta_signal_timing_summary.tsv", outputDir));
    summary << "Z\tbeta_center\tN\tsigma_invbeta_raw\tsigma_invbeta_betaH"
               "\tsigma_t_eff_ns\tmedian_edep_MeV_per_layer\tmedian_anode_ADC_per_side\n";
    const double width = (kBetaMax - kBetaMin) / kNBetaBins;
    for (const auto &chargeData : data)
        for (int bin = 0; bin < kNBetaBins; ++bin)
        {
            const auto &item = chargeData.bins[bin];
            if (static_cast<int>(item.rawInvBetaResidual.size()) < kMinEntries)
                continue;
            summary << chargeData.charge << '\t'
                    << kBetaMin + (bin + 0.5) * width << '\t'
                    << item.rawInvBetaResidual.size() << '\t'
                    << robustSigma(item.rawInvBetaResidual) << '\t'
                    << robustSigma(item.betaHInvBetaResidual) << '\t'
                    << robustSigma(item.timingEquivalent) << '\t'
                    << median(item.meanEdep) << '\t'
                    << median(item.meanAdc) << '\n';
        }
    summary.close();

    std::array<Series, 4> rawSigma, betaHSigma, timing, edep, adc, timingEdep, timingAdc;
    for (size_t i = 0; i < charges.size(); ++i)
    {
        rawSigma[i] = sigmaSeries(data[i], false);
        betaHSigma[i] = sigmaSeries(data[i], true);
        timing[i] = timingSeries(data[i]);
        edep[i] = medianSeries(data[i], false);
        adc[i] = medianSeries(data[i], true);
        timingEdep[i] = timingVsSignal(data[i], false);
        timingAdc[i] = timingVsSignal(data[i], true);
    }

    TString report = TString(outputDir) + "/beta_signal_timing_report.pdf";

    // Page 1: confirm that the beta dependence is already present in raw hit times.
    TCanvas residualCanvas("residualCanvas", "", 1500, 1050);
    residualCanvas.Divide(2, 2, 0.003, 0.003);
    for (size_t i = 0; i < charges.size(); ++i)
    {
        residualCanvas.cd(i + 1);
        configurePad();
        const std::vector<Series> local = {rawSigma[i], betaHSigma[i]};
        const auto yrange = rangeFor(local, 0.0);
        auto *frame = gPad->DrawFrame(kBetaMin, 0.0, kBetaMax, yrange.second);
        styleAxis(frame, "#beta_{MC} at center",
                  "core #sigma(1/#beta - 1/#beta_{MC})");
        frame->SetTitle(Form("Z = %d, four-hit N = %lld", charges[i], data[i].selected));

        auto *raw = makeGraph(rawSigma[i], colors[i], 20,
                              Form("raw_sigma_Z%d", charges[i]));
        auto *betah = makeGraph(betaHSigma[i], TColor::GetColor("#333333"), 24,
                                Form("betah_sigma_Z%d", charges[i]));
        raw->Draw("LP SAME");
        betah->SetLineStyle(2);
        betah->Draw("LP SAME");
        if (i == 0)
        {
            auto *legend = new TLegend(0.52, 0.69, 0.92, 0.86);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->SetTextSize(0.037);
            legend->AddEntry(raw, "4-hit raw-time refit", "lp");
            legend->AddEntry(betah, "stored BetaH", "lp");
            legend->Draw();
        }
    }
    residualCanvas.SaveAs(Form("%s/sigma_invbeta_vs_beta.pdf", outputDir));
    residualCanvas.SaveAs(Form("%s/sigma_invbeta_vs_beta.png", outputDir));

    // Page 2: inferred timing resolution versus beta, with the stored ETime reference.
    TCanvas timingCanvas("timingCanvas", "", 1200, 900);
    configurePad();
    std::vector<Series> timingVector(timing.begin(), timing.end());
    const auto timingRange = rangeFor(timingVector, 0.0);
    auto *timingFrame = gPad->DrawFrame(kBetaMin, 0.0, kBetaMax,
                                        std::max(timingRange.second, 0.18));
    styleAxis(timingFrame, "#beta_{MC} at center", "effective #sigma_{t} [ns]");
    timingFrame->SetTitle("Effective TOF timing width from raw four-hit slope");
    auto *timingLegend = new TLegend(0.66, 0.59, 0.91, 0.86);
    timingLegend->SetBorderSize(0);
    timingLegend->SetFillStyle(0);
    timingLegend->SetTextSize(0.038);
    for (size_t i = 0; i < charges.size(); ++i)
    {
        auto *graph = makeGraph(timing[i], colors[i], markers[i],
                                Form("timing_Z%d", charges[i]));
        graph->Draw("LP SAME");
        timingLegend->AddEntry(graph, Form("Z = %d", charges[i]), "lp");
    }
    TLine storedError(kBetaMin, kStoredTofTimeErrorNs, kBetaMax, kStoredTofTimeErrorNs);
    storedError.SetLineColor(TColor::GetColor("#555555"));
    storedError.SetLineStyle(7);
    storedError.SetLineWidth(2);
    storedError.DrawClone("SAME");
    timingLegend->AddEntry(&storedError, "stored ETime = 154.481 ps", "l");
    timingLegend->Draw();
    timingCanvas.SaveAs(Form("%s/sigma_t_eff_vs_beta.pdf", outputDir));
    timingCanvas.SaveAs(Form("%s/sigma_t_eff_vs_beta.png", outputDir));

    // Page 3: signal scale versus beta.
    TCanvas signalCanvas("signalCanvas", "", 1500, 650);
    signalCanvas.Divide(2, 1, 0.003, 0.003);
    signalCanvas.cd(1);
    configurePad();
    std::vector<Series> edepVector(edep.begin(), edep.end());
    const auto edepRange = rangeFor(edepVector, 0.0);
    auto *edepFrame = gPad->DrawFrame(kBetaMin, 0.0, kBetaMax, edepRange.second);
    styleAxis(edepFrame, "#beta_{MC} at center", "median TOF E_{dep} / layer [MeV]");
    edepFrame->SetTitle("Reconstructed TOF energy deposit");
    auto *signalLegend = new TLegend(0.69, 0.63, 0.91, 0.87);
    signalLegend->SetBorderSize(0);
    signalLegend->SetFillStyle(0);
    signalLegend->SetTextSize(0.040);
    for (size_t i = 0; i < charges.size(); ++i)
    {
        auto *graph = makeGraph(edep[i], colors[i], markers[i],
                                Form("edep_Z%d", charges[i]));
        graph->Draw("LP SAME");
        signalLegend->AddEntry(graph, Form("Z = %d", charges[i]), "lp");
    }
    signalLegend->Draw();

    signalCanvas.cd(2);
    configurePad();
    std::vector<Series> adcVector(adc.begin(), adc.end());
    const auto adcRange = rangeFor(adcVector, 0.0);
    auto *adcFrame = gPad->DrawFrame(kBetaMin, 0.0, kBetaMax, adcRange.second);
    styleAxis(adcFrame, "#beta_{MC} at center", "median anode ADC / side");
    adcFrame->SetTitle("Observed TOF anode signal");
    for (size_t i = 0; i < charges.size(); ++i)
        makeGraph(adc[i], colors[i], markers[i], Form("adc_Z%d", charges[i]))->Draw("LP SAME");
    signalLegend->DrawClone();
    signalCanvas.SaveAs(Form("%s/signal_vs_beta.pdf", outputDir));
    signalCanvas.SaveAs(Form("%s/signal_vs_beta.png", outputDir));

    // Page 4: timing resolution against signal itself.
    TCanvas relationCanvas("relationCanvas", "", 1500, 650);
    relationCanvas.Divide(2, 1, 0.003, 0.003);
    relationCanvas.cd(1);
    configurePad();
    std::vector<Series> timingEdepVector(timingEdep.begin(), timingEdep.end());
    double edepXmax = 0;
    for (const auto &item : timingEdep)
        for (double value : item.x)
            edepXmax = std::max(edepXmax, value);
    auto *relationEdepFrame = gPad->DrawFrame(0.0, 0.0, edepXmax * 1.08, 0.19);
    styleAxis(relationEdepFrame, "median TOF E_{dep} / layer [MeV]", "effective #sigma_{t} [ns]");
    relationEdepFrame->SetTitle("Timing width vs reconstructed TOF deposit");
    auto *relationLegend = new TLegend(0.68, 0.63, 0.91, 0.87);
    relationLegend->SetBorderSize(0);
    relationLegend->SetFillStyle(0);
    relationLegend->SetTextSize(0.040);
    for (size_t i = 0; i < charges.size(); ++i)
    {
        auto *graph = makeGraph(timingEdep[i], colors[i], markers[i],
                                Form("timing_edep_Z%d", charges[i]));
        graph->Draw("LP SAME");
        relationLegend->AddEntry(graph, Form("Z = %d", charges[i]), "lp");
    }
    relationLegend->Draw();

    relationCanvas.cd(2);
    configurePad();
    double adcXmax = 0;
    for (const auto &item : timingAdc)
        for (double value : item.x)
            adcXmax = std::max(adcXmax, value);
    auto *relationAdcFrame = gPad->DrawFrame(0.0, 0.0, adcXmax * 1.08, 0.19);
    styleAxis(relationAdcFrame, "median anode ADC / side", "effective #sigma_{t} [ns]");
    relationAdcFrame->SetTitle("Timing width vs observed anode signal");
    for (size_t i = 0; i < charges.size(); ++i)
        makeGraph(timingAdc[i], colors[i], markers[i],
                  Form("timing_adc_Z%d", charges[i]))->Draw("LP SAME");
    relationLegend->DrawClone();
    relationCanvas.SaveAs(Form("%s/sigma_t_eff_vs_signal.pdf", outputDir));
    relationCanvas.SaveAs(Form("%s/sigma_t_eff_vs_signal.png", outputDir));

    // Assemble a clean multi-page report only after the individual files are closed.
    residualCanvas.Print(report + "[");
    residualCanvas.Print(report);
    timingCanvas.Print(report);
    signalCanvas.Print(report);
    relationCanvas.Print(report);
    relationCanvas.Print(report + "]");

    TFile rootOutput(Form("%s/beta_signal_timing.root", outputDir), "RECREATE");
    for (size_t i = 0; i < charges.size(); ++i)
    {
        makeGraph(rawSigma[i], colors[i], markers[i], Form("sigma_invbeta_raw_Z%d", charges[i]))->Write();
        makeGraph(betaHSigma[i], colors[i], markers[i], Form("sigma_invbeta_betaH_Z%d", charges[i]))->Write();
        makeGraph(timing[i], colors[i], markers[i], Form("sigma_t_eff_Z%d", charges[i]))->Write();
        makeGraph(edep[i], colors[i], markers[i], Form("median_edep_Z%d", charges[i]))->Write();
        makeGraph(adc[i], colors[i], markers[i], Form("median_adc_Z%d", charges[i]))->Write();
        makeGraph(timingEdep[i], colors[i], markers[i], Form("sigma_t_vs_edep_Z%d", charges[i]))->Write();
        makeGraph(timingAdc[i], colors[i], markers[i], Form("sigma_t_vs_adc_Z%d", charges[i]))->Write();
    }
    rootOutput.Close();

    std::cout << "Saved report and numerical outputs under " << outputDir << std::endl;
}
