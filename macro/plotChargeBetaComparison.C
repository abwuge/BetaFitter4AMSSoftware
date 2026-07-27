#include <TCanvas.h>
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

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "beta_fitter_macro_utils.h"

namespace
{
constexpr int kNBetaBins = 28;
constexpr int kNMethods = 2;
constexpr double kResidualMin = -0.5;
constexpr double kResidualMax = 1.5;
constexpr double kIqrToSigma = 0.741301109252801;

struct FitSeries
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
    double betaMin = 0;
    double betaMax = 1;
    std::array<FitSeries, kNMethods> methods;
};

bool fitCore(TH1D *hist, double &mean, double &meanError,
             double &sigma, double &sigmaError)
{
    if (!hist || hist->GetEntries() < 300)
        return false;

    const double probabilities[3] = {0.25, 0.50, 0.75};
    double quantiles[3] = {};
    hist->GetQuantiles(3, quantiles, probabilities);

    const double robustSigma = kIqrToSigma * (quantiles[2] - quantiles[0]);
    if (!(robustSigma > 0) || !std::isfinite(robustSigma))
        return false;

    const double fitMin = std::max(kResidualMin, quantiles[1] - 3 * robustSigma);
    const double fitMax = std::min(kResidualMax, quantiles[1] + 3 * robustSigma);
    TF1 gaussian(Form("gaussian_%s", hist->GetName()), "gaus", fitMin, fitMax);
    gaussian.SetParameters(hist->GetMaximum(), quantiles[1], robustSigma);
    gaussian.SetLineWidth(3);
    gaussian.SetLineStyle(2);

    const TFitResultPtr result = hist->Fit(&gaussian, "SQNR");
    if (!result.Get() || result->Status() != 0)
        return false;

    mean = result->Parameter(1);
    meanError = result->ParError(1);
    sigma = std::abs(result->Parameter(2));
    sigmaError = result->ParError(2);
    return std::isfinite(mean) && std::isfinite(meanError) &&
           std::isfinite(sigma) && std::isfinite(sigmaError);
}

ChargeResult analyzeCharge(const char *fileName, int charge)
{
    ChargeResult output;
    output.charge = charge;

    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: cannot open " << fileName << std::endl;
        return output;
    }

    TTree *tree = static_cast<TTree *>(file->Get("betaTree"));
    if (!tree || tree->GetEntries() == 0)
    {
        std::cerr << "Error: missing or empty betaTree in " << fileName << std::endl;
        return output;
    }

    output.betaMin = std::max(0.25, tree->GetMinimum("mcBeta"));
    output.betaMax = std::min(1.0, tree->GetMaximum("mcBeta"));
    const double betaWidth = (output.betaMax - output.betaMin) / kNBetaBins;

    std::array<std::array<TH1D *, kNBetaBins>, kNMethods> histograms;
    for (int method = 0; method < kNMethods; ++method)
        for (int bin = 0; bin < kNBetaBins; ++bin)
            histograms[method][bin] = new TH1D(
                Form("hResidual_Z%d_M%d_B%d", charge, method, bin), "",
                800, kResidualMin, kResidualMax);

    double mcBeta = 0;
    double linearBeta = 0;
    double nonlinearBeta = 0;
    tree->SetBranchAddress("mcBeta", &mcBeta);
    tree->SetBranchAddress("linearBeta", &linearBeta);
    tree->SetBranchAddress("nonlinearBeta", &nonlinearBeta);

    const Long64_t entries = tree->GetEntries();
    for (Long64_t entry = 0; entry < entries; ++entry)
    {
        tree->GetEntry(entry);
        if (!(mcBeta > 0) || !std::isfinite(mcBeta) ||
            mcBeta < output.betaMin || mcBeta > output.betaMax)
            continue;

        int bin = static_cast<int>((mcBeta - output.betaMin) / betaWidth);
        bin = std::max(0, std::min(kNBetaBins - 1, bin));
        const double reconstructed[kNMethods] = {nonlinearBeta, linearBeta};
        for (int method = 0; method < kNMethods; ++method)
        {
            if (!(reconstructed[method] > 0) || !std::isfinite(reconstructed[method]))
                continue;
            const double residual = 1.0 / reconstructed[method] - 1.0 / mcBeta;
            if (std::isfinite(residual))
                histograms[method][bin]->Fill(residual);
        }
    }

    for (int method = 0; method < kNMethods; ++method)
    {
        FitSeries &series = output.methods[method];
        for (int bin = 0; bin < kNBetaBins; ++bin)
        {
            double mean = 0;
            double meanError = 0;
            double sigma = 0;
            double sigmaError = 0;
            if (!fitCore(histograms[method][bin], mean, meanError, sigma, sigmaError))
                continue;

            series.beta.push_back(output.betaMin + (bin + 0.5) * betaWidth);
            series.mean.push_back(mean);
            series.meanError.push_back(meanError);
            series.sigma.push_back(sigma);
            series.sigmaError.push_back(sigmaError);
        }
    }

    std::cout << "Z=" << charge << ": " << entries << " entries, "
              << output.methods[0].beta.size() << " nonlinear bins, "
              << output.methods[1].beta.size() << " linear bins" << std::endl;
    file->Close();
    return output;
}

TGraphErrors *makeGraph(const FitSeries &series, bool drawSigma,
                        int color, int markerStyle)
{
    const std::vector<double> &values = drawSigma ? series.sigma : series.mean;
    const std::vector<double> &errors = drawSigma ? series.sigmaError : series.meanError;
    TGraphErrors *graph = new TGraphErrors(
        series.beta.size(), series.beta.data(), values.data(), nullptr, errors.data());
    graph->SetLineColor(color);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerSize(0.85);
    graph->SetLineWidth(2);
    return graph;
}

std::pair<double, double> graphRange(const ChargeResult &result, bool sigma)
{
    double minimum = sigma ? 0.0 : std::numeric_limits<double>::max();
    double maximum = sigma ? 0.0 : -std::numeric_limits<double>::max();
    for (int method = 0; method < kNMethods; ++method)
    {
        const FitSeries &series = result.methods[method];
        const std::vector<double> &values = sigma ? series.sigma : series.mean;
        const std::vector<double> &errors = sigma ? series.sigmaError : series.meanError;
        for (size_t i = 0; i < values.size(); ++i)
        {
            minimum = std::min(minimum, values[i] - errors[i]);
            maximum = std::max(maximum, values[i] + errors[i]);
        }
    }

    if (sigma)
        return std::make_pair(0.0, std::max(0.01, maximum * 1.18));

    minimum = std::min(minimum, 0.0);
    maximum = std::max(maximum, 0.0);
    const double padding = std::max(0.005, 0.12 * (maximum - minimum));
    return std::make_pair(minimum - padding, maximum + padding);
}
} // namespace

void plotChargeBetaComparison(
    const char *z2File = "results/11.root",
    const char *z6File = "results/12.root",
    const char *z8File = "results/13.root",
    const char *outputName = nullptr)
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(4);

    const std::array<const char *, 3> files = {z2File, z6File, z8File};
    const std::array<int, 3> charges = {2, 6, 8};
    std::array<ChargeResult, 3> results;
    for (int index = 0; index < 3; ++index)
        results[index] = analyzeCharge(files[index], charges[index]);

    TString resolvedOutput = BetaFitterMacro::output_path(
        outputName, "beta_comparison", "charge_beta_center_comparison.pdf");
    TString pngOutput = resolvedOutput;
    if (pngOutput.EndsWith(".pdf"))
        pngOutput.ReplaceAll(".pdf", ".png");
    else
        pngOutput += ".png";

    const int nonlinearColor = TColor::GetColor("#45364B");
    const int linearColor = TColor::GetColor("#5F6077");
    TCanvas canvas("chargeBetaComparison", "", 1500, 860);
    canvas.Divide(3, 2, 0.002, 0.002);

    for (int chargeIndex = 0; chargeIndex < 3; ++chargeIndex)
    {
        const ChargeResult &result = results[chargeIndex];
        for (int row = 0; row < 2; ++row)
        {
            canvas.cd(row * 3 + chargeIndex + 1);
            gPad->SetLeftMargin(chargeIndex == 0 ? 0.16 : 0.12);
            gPad->SetRightMargin(0.035);
            gPad->SetBottomMargin(row == 1 ? 0.15 : 0.11);
            gPad->SetTopMargin(row == 0 ? 0.12 : 0.06);
            gPad->SetGridx();
            gPad->SetGridy();

            const bool sigma = row == 1;
            const std::pair<double, double> range = graphRange(result, sigma);
            TH1F *frame = gPad->DrawFrame(result.betaMin, range.first,
                                          result.betaMax, range.second);
            frame->SetTitle(row == 0 ? Form("Z = %d", result.charge) : "");
            frame->GetXaxis()->SetTitle("#beta_{MC} at AMS center");
            frame->GetYaxis()->SetTitle(sigma ? "Core #sigma(1/#beta residual)"
                                               : "Core #mu(1/#beta residual)");
            frame->GetXaxis()->SetTitleSize(0.052);
            frame->GetYaxis()->SetTitleSize(0.052);
            frame->GetXaxis()->SetLabelSize(0.043);
            frame->GetYaxis()->SetLabelSize(0.043);
            frame->GetYaxis()->SetTitleOffset(chargeIndex == 0 ? 1.45 : 1.20);

            TLine zero(result.betaMin, 0, result.betaMax, 0);
            if (!sigma)
            {
                zero.SetLineColor(TColor::GetColor("#788AA3"));
                zero.SetLineStyle(2);
                zero.SetLineWidth(2);
                zero.Draw("SAME");
            }

            TGraphErrors *nonlinear = makeGraph(result.methods[0], sigma,
                                                 nonlinearColor, 20);
            TGraphErrors *linear = makeGraph(result.methods[1], sigma,
                                              linearColor, 21);
            nonlinear->Draw("LP SAME");
            linear->Draw("LP SAME");

            if (row == 0)
            {
                TLatex chargeLabel;
                chargeLabel.SetNDC();
                chargeLabel.SetTextFont(42);
                chargeLabel.SetTextSize(0.060);
                chargeLabel.SetTextAlign(22);
                chargeLabel.DrawLatex(0.52, 0.94, Form("Z = %d", result.charge));
            }

            if (row == 0 && chargeIndex == 0)
            {
                TLegend *legend = new TLegend(0.58, 0.67, 0.91, 0.84);
                legend->SetBorderSize(0);
                legend->SetFillStyle(0);
                legend->SetTextSize(0.045);
                legend->AddEntry(nonlinear, "Nonlinear", "lp");
                legend->AddEntry(linear, "Linear", "lp");
                legend->Draw();
            }
        }
    }

    canvas.SaveAs(resolvedOutput);
    canvas.SaveAs(pngOutput);
    std::cout << "Saved comparison to " << resolvedOutput << " and " << pngOutput << std::endl;
}
