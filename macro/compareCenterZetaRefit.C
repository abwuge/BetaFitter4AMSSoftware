#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace
{
constexpr int kNMethods = 3;
constexpr int kNBetaBins = 10;
constexpr double kResidualMin = -0.5;
constexpr double kResidualMax = 1.5;
constexpr double kIqrToSigma = 0.741301109252801;

const char *kMethodNames[kNMethods] = {"linear", "nonlinear #zeta=1.9", "nonlinear refit #zeta"};

double quantile(std::vector<double> values, double probability)
{
    if (values.empty())
        return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    const double position = probability * (values.size() - 1);
    const size_t lower = static_cast<size_t>(position);
    const size_t upper = std::min(lower + 1, values.size() - 1);
    const double fraction = position - lower;
    return values[lower] * (1 - fraction) + values[upper] * fraction;
}

struct Metrics
{
    bool fitValid = false;
    double mean = 0;
    double meanError = 0;
    double sigma = 0;
    double sigmaError = 0;
    double median = 0;
    double width68 = 0;
    double width95 = 0;
};

Metrics summarize(TH1D *hist, const std::vector<double> &values)
{
    Metrics output;
    if (!hist || values.size() < 40)
        return output;

    const double q16 = quantile(values, 0.16);
    const double q25 = quantile(values, 0.25);
    const double q50 = quantile(values, 0.50);
    const double q75 = quantile(values, 0.75);
    const double q84 = quantile(values, 0.84);
    const double q025 = quantile(values, 0.025);
    const double q975 = quantile(values, 0.975);
    output.median = q50;
    output.width68 = (q84 - q16) / 2;
    output.width95 = (q975 - q025) / 2;

    const double robustSigma = kIqrToSigma * (q75 - q25);
    if (!(robustSigma > 0) || !std::isfinite(robustSigma))
        return output;
    const double fitMin = std::max(kResidualMin, q50 - 3 * robustSigma);
    const double fitMax = std::min(kResidualMax, q50 + 3 * robustSigma);
    TF1 gaussian(Form("gaussian_%s", hist->GetName()), "gaus", fitMin, fitMax);
    gaussian.SetParameters(hist->GetMaximum(), q50, robustSigma);
    const TFitResultPtr fit = hist->Fit(&gaussian, "SQNR");
    if (!fit.Get() || fit->Status() != 0)
        return output;
    output.mean = fit->Parameter(1);
    output.meanError = fit->ParError(1);
    output.sigma = std::abs(fit->Parameter(2));
    output.sigmaError = fit->ParError(2);
    output.fitValid = std::isfinite(output.mean) && std::isfinite(output.sigma) &&
                      output.sigma < 0.08 && output.meanError < 0.02 &&
                      output.sigmaError < 0.5 * output.sigma;
    return output;
}

struct Series
{
    std::vector<double> beta;
    std::vector<double> mean;
    std::vector<double> meanError;
    std::vector<double> sigma;
    std::vector<double> sigmaError;
};

struct ChargeResult
{
    int charge = 0;
    double betaLow = 0;
    double betaHigh = 0;
    double fittedZeta = 0;
    Long64_t selected = 0;
    Long64_t common = 0;
    std::array<Long64_t, kNMethods> failures = {};
    std::array<TH1D *, kNMethods> residualHistograms = {};
    std::array<std::vector<double>, kNMethods> residuals;
    std::array<Metrics, kNMethods> metrics;
    std::array<Series, kNMethods> binned;
};

ChargeResult analyzeCharge(const char *baselinePath, const char *refitPath,
                           int charge, double betaLow, double betaHigh,
                           double fittedZeta)
{
    ChargeResult output;
    output.charge = charge;
    output.betaLow = betaLow;
    output.betaHigh = betaHigh;
    output.fittedZeta = fittedZeta;

    TFile baselineFile(baselinePath, "READ");
    TFile refitFile(refitPath, "READ");
    TTree *baseline = static_cast<TTree *>(baselineFile.Get("betaTree"));
    TTree *refit = static_cast<TTree *>(refitFile.Get("betaTree"));
    if (!baseline || !refit || baseline->GetEntries() != refit->GetEntries())
    {
        std::cerr << "Error: missing or unaligned betaTree for Z=" << charge << std::endl;
        return output;
    }

    double mcBaseline = 0, linearBaseline = 0, nonlinearBaseline = 0;
    double mcRefit = 0, linearRefit = 0, nonlinearRefit = 0;
    baseline->SetBranchAddress("mcBeta", &mcBaseline);
    baseline->SetBranchAddress("linearBeta", &linearBaseline);
    baseline->SetBranchAddress("nonlinearBeta", &nonlinearBaseline);
    refit->SetBranchAddress("mcBeta", &mcRefit);
    refit->SetBranchAddress("linearBeta", &linearRefit);
    refit->SetBranchAddress("nonlinearBeta", &nonlinearRefit);

    std::array<std::array<TH1D *, kNBetaBins>, kNMethods> binHistograms;
    std::array<std::array<std::vector<double>, kNBetaBins>, kNMethods> binValues;
    for (int method = 0; method < kNMethods; ++method)
    {
        output.residualHistograms[method] = new TH1D(
            Form("hResidual_Z%d_M%d", charge, method), "", 1000,
            kResidualMin, kResidualMax);
        output.residualHistograms[method]->SetDirectory(nullptr);
        for (int bin = 0; bin < kNBetaBins; ++bin)
        {
            binHistograms[method][bin] = new TH1D(
                Form("hResidual_Z%d_M%d_B%d", charge, method, bin), "", 600,
                kResidualMin, kResidualMax);
            binHistograms[method][bin]->SetDirectory(nullptr);
        }
    }

    const double plotBetaHigh = 0.90;
    const double betaWidth = (plotBetaHigh - betaLow) / kNBetaBins;
    for (Long64_t entry = 0; entry < baseline->GetEntries(); ++entry)
    {
        baseline->GetEntry(entry);
        refit->GetEntry(entry);
        if (std::abs(mcBaseline - mcRefit) > 1e-12 ||
            std::abs(linearBaseline - linearRefit) > 1e-12)
        {
            std::cerr << "Error: entry mismatch for Z=" << charge
                      << " at " << entry << std::endl;
            return output;
        }
        if (!std::isfinite(mcBaseline) || mcBaseline < betaLow || mcBaseline >= plotBetaHigh)
            continue;

        const double reconstructed[kNMethods] = {linearBaseline, nonlinearBaseline, nonlinearRefit};
        bool valid[kNMethods] = {};
        for (int method = 0; method < kNMethods; ++method)
            valid[method] = reconstructed[method] > 0 && std::isfinite(reconstructed[method]);
        const bool inSummaryRange = mcBaseline < betaHigh;
        if (inSummaryRange)
        {
            ++output.selected;
            for (int method = 0; method < kNMethods; ++method)
                if (!valid[method])
                    ++output.failures[method];
        }
        if (!(valid[0] && valid[1] && valid[2]))
            continue;
        if (inSummaryRange)
            ++output.common;

        int bin = std::min(kNBetaBins - 1,
                           static_cast<int>((mcBaseline - betaLow) / betaWidth));
        for (int method = 0; method < kNMethods; ++method)
        {
            const double residual = 1.0 / reconstructed[method] - 1.0 / mcBaseline;
            if (!std::isfinite(residual))
                continue;
            binHistograms[method][bin]->Fill(residual);
            binValues[method][bin].push_back(residual);
            if (inSummaryRange)
            {
                output.residualHistograms[method]->Fill(residual);
                output.residuals[method].push_back(residual);
            }
        }
    }

    for (int method = 0; method < kNMethods; ++method)
    {
        output.metrics[method] = summarize(output.residualHistograms[method],
                                           output.residuals[method]);
        for (int bin = 0; bin < kNBetaBins; ++bin)
        {
            const Metrics metrics = summarize(binHistograms[method][bin],
                                              binValues[method][bin]);
            if (!metrics.fitValid)
                continue;
            Series &series = output.binned[method];
            series.beta.push_back(betaLow + (bin + 0.5) * betaWidth);
            series.mean.push_back(metrics.mean);
            series.meanError.push_back(metrics.meanError);
            series.sigma.push_back(metrics.sigma);
            series.sigmaError.push_back(metrics.sigmaError);
        }
    }
    return output;
}

TGraphErrors *makeGraph(const Series &series, bool sigma, int color, int marker)
{
    const std::vector<double> &values = sigma ? series.sigma : series.mean;
    const std::vector<double> &errors = sigma ? series.sigmaError : series.meanError;
    TGraphErrors *graph = new TGraphErrors(series.beta.size(), series.beta.data(),
                                           values.data(), nullptr, errors.data());
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(0.8);
    graph->SetLineWidth(2);
    return graph;
}
}

void compareCenterZetaRefit(
    const char *z2Baseline = "results/validation_center_zeta/beta_Z2_zeta1p9.root",
    const char *z2Refit = "results/validation_center_zeta/beta_Z2_zetafit.root",
    const char *z6Baseline = "results/validation_center_zeta/beta_Z6_zeta1p9.root",
    const char *z6Refit = "results/validation_center_zeta/beta_Z6_zetafit.root",
    const char *z8Baseline = "results/validation_center_zeta/beta_Z8_zeta1p9.root",
    const char *z8Refit = "results/validation_center_zeta/beta_Z8_zetafit.root",
    const char *outputDirectory = "results/validation_center_zeta")
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);
    gSystem->mkdir(outputDirectory, true);

    const std::array<const char *, 3> baselinePaths = {z2Baseline, z6Baseline, z8Baseline};
    const std::array<const char *, 3> refitPaths = {z2Refit, z6Refit, z8Refit};
    const std::array<int, 3> charges = {2, 6, 8};
    const std::array<double, 3> betaLow = {0.35, 0.40, 0.45};
    const std::array<double, 3> betaHigh = {0.50, 0.60, 0.60};
    const std::array<double, 3> fittedZeta = {3.265, 2.685, 2.275};
    std::array<ChargeResult, 3> results;
    for (int i = 0; i < 3; ++i)
        results[i] = analyzeCharge(baselinePaths[i], refitPaths[i], charges[i],
                                   betaLow[i], betaHigh[i], fittedZeta[i]);

    std::ofstream table((std::string(outputDirectory) + "/comparison_summary.tsv").c_str());
    table << "Z\tbeta_low\tbeta_high\tfitted_zeta\tmethod\tselected\tcommon\tfailures\tfailure_rate\tcore_mu\tcore_mu_error\tcore_sigma\tcore_sigma_error\tmedian\twidth68\twidth95\n";
    table << std::setprecision(8);
    for (const ChargeResult &result : results)
        for (int method = 0; method < kNMethods; ++method)
        {
            const Metrics &metrics = result.metrics[method];
            table << result.charge << '\t' << result.betaLow << '\t' << result.betaHigh << '\t'
                  << result.fittedZeta << '\t' << kMethodNames[method] << '\t'
                  << result.selected << '\t' << result.common << '\t' << result.failures[method] << '\t'
                  << (result.selected ? 1.0 * result.failures[method] / result.selected : 0) << '\t'
                  << metrics.mean << '\t' << metrics.meanError << '\t'
                  << metrics.sigma << '\t' << metrics.sigmaError << '\t'
                  << metrics.median << '\t' << metrics.width68 << '\t' << metrics.width95 << '\n';
        }

    const int colors[kNMethods] = {
        TColor::GetColor("#45364B"), TColor::GetColor("#5F6077"),
        TColor::GetColor("#788AA3")};
    const int markers[kNMethods] = {20, 21, 22};
    TCanvas canvas("centerZetaComparison", "", 1800, 1350);
    canvas.Divide(3, 3, 0.002, 0.002);
    for (int chargeIndex = 0; chargeIndex < 3; ++chargeIndex)
    {
        const ChargeResult &result = results[chargeIndex];
        canvas.cd(chargeIndex + 1);
        gPad->SetLeftMargin(0.17);
        gPad->SetRightMargin(0.035);
        gPad->SetTopMargin(0.12);
        gPad->SetBottomMargin(0.16);
        double maximum = 0;
        for (int method = 0; method < kNMethods; ++method)
        {
            TH1D *hist = result.residualHistograms[method];
            if (hist->Integral() > 0)
                hist->Scale(1.0 / hist->Integral(), "width");
            hist->SetLineColor(colors[method]);
            hist->SetLineWidth(3);
            maximum = std::max(maximum, hist->GetMaximum());
        }
        TH1F *frame = gPad->DrawFrame(-0.12, 0, 0.18, maximum * 1.25);
        frame->SetTitle(Form("Z = %d, %.2f #leq #beta_{MC} < %.2f",
                             result.charge, result.betaLow, result.betaHigh));
        frame->GetXaxis()->SetTitle("1/#beta_{rec} - 1/#beta_{MC}");
        frame->GetYaxis()->SetTitle("Density");
        frame->GetXaxis()->SetTitleSize(0.042);
        frame->GetYaxis()->SetTitleSize(0.042);
        frame->GetXaxis()->SetLabelSize(0.034);
        frame->GetYaxis()->SetLabelSize(0.034);
        frame->GetYaxis()->SetTitleOffset(1.55);
        for (int method = 0; method < kNMethods; ++method)
            result.residualHistograms[method]->Draw("HIST SAME");
        TLatex chargeLabel;
        chargeLabel.SetNDC(true);
        chargeLabel.SetTextSize(0.032);
        chargeLabel.DrawLatex(0.19, 0.93,
                              Form("Z = %d, %.2f #leq #beta_{MC} < %.2f",
                                   result.charge, result.betaLow, result.betaHigh));
        TLegend *legend = new TLegend(0.36, 0.57, 0.95, 0.88);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.027);
        legend->SetHeader(Form("common N = %lld", result.common), "C");
        for (int method = 0; method < kNMethods; ++method)
        {
            const Metrics &m = result.metrics[method];
            legend->AddEntry(result.residualHistograms[method],
                             Form("%s: #mu=%.4f, #sigma=%.4f",
                                  method == 0 ? "linear" : (method == 1 ? "nonlinear (1.9)" : "nonlinear (refit)"),
                                  m.mean, m.sigma), "l");
        }
        legend->Draw();

        for (int row = 1; row < 3; ++row)
        {
            canvas.cd(row * 3 + chargeIndex + 1);
            gPad->SetLeftMargin(0.17);
            gPad->SetRightMargin(0.035);
            gPad->SetTopMargin(0.06);
            gPad->SetBottomMargin(row == 2 ? 0.17 : 0.13);
            gPad->SetGridx();
            gPad->SetGridy();
            const bool sigma = row == 2;
            double minimum = sigma ? 0 : std::numeric_limits<double>::max();
            double maximumValue = sigma ? 0 : -std::numeric_limits<double>::max();
            for (int method = 0; method < kNMethods; ++method)
            {
                const Series &series = result.binned[method];
                const std::vector<double> &values = sigma ? series.sigma : series.mean;
                const std::vector<double> &errors = sigma ? series.sigmaError : series.meanError;
                for (size_t i = 0; i < values.size(); ++i)
                {
                    minimum = std::min(minimum, values[i] - errors[i]);
                    maximumValue = std::max(maximumValue, values[i] + errors[i]);
                }
            }
            if (!sigma)
            {
                minimum = std::min(minimum, 0.0);
                maximumValue = std::max(maximumValue, 0.0);
                const double padding = std::max(0.005, 0.15 * (maximumValue - minimum));
                minimum -= padding;
                maximumValue += padding;
            }
            else
                maximumValue = std::max(0.02, 1.18 * maximumValue);
            TH1F *axis = gPad->DrawFrame(result.betaLow, minimum, 0.90, maximumValue);
            axis->GetXaxis()->SetTitle("#beta_{MC} at AMS center");
            axis->GetYaxis()->SetTitle(sigma ? "Core #sigma(1/#beta residual)"
                                              : "Core #mu(1/#beta residual)");
            axis->GetXaxis()->SetTitleSize(0.042);
            axis->GetYaxis()->SetTitleSize(0.042);
            axis->GetXaxis()->SetLabelSize(0.034);
            axis->GetYaxis()->SetLabelSize(0.034);
            axis->GetYaxis()->SetTitleOffset(1.55);
            if (!sigma)
            {
                TLine *zero = new TLine(result.betaLow, 0, 0.90, 0);
                zero->SetLineColor(TColor::GetColor("#95670B"));
                zero->SetLineStyle(2);
                zero->SetLineWidth(2);
                zero->Draw("SAME");
            }
            for (int method = 0; method < kNMethods; ++method)
                makeGraph(result.binned[method], sigma, colors[method], markers[method])->Draw("LP SAME");
            if (chargeIndex == 0 && row == 1)
            {
                TLegend *methodLegend = new TLegend(0.19, 0.68, 0.66, 0.91);
                methodLegend->SetBorderSize(0);
                methodLegend->SetFillStyle(0);
                methodLegend->SetTextSize(0.036);
                for (int method = 0; method < kNMethods; ++method)
                {
                    TGraphErrors *sample = makeGraph(result.binned[method], sigma,
                                                     colors[method], markers[method]);
                    methodLegend->AddEntry(sample, kMethodNames[method], "lp");
                }
                methodLegend->Draw();
            }
        }
    }

    const std::string pdf = std::string(outputDirectory) + "/center_zeta_refit_comparison.pdf";
    const std::string png = std::string(outputDirectory) + "/center_zeta_refit_comparison.png";
    canvas.Print(pdf.c_str());
    canvas.Print(png.c_str());
    std::cout << "Wrote comparison summary and plots" << std::endl;
}
