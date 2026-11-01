#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
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
constexpr int kNCharges = 4;
constexpr int kCharges[kNCharges] = {1, 2, 6, 8};
constexpr int kNBetaBins = 8;
constexpr int kNMethods = 2;
constexpr double kBetaMin = 0.60;
constexpr double kBetaWidth = 0.05;
constexpr double kResidualMin = -0.16;
constexpr double kResidualMax = 0.16;
constexpr double kIqrToSigma = 0.741301109252801;
constexpr int kColors[kNMethods] = {kP6Blue, kP6Yellow};
const char *kLabels[kNMethods] = {"Hybrid", "Data-trained"};

struct FitResult
{
    bool valid = false;
    double mean = 0;
    double meanError = 0;
    double sigma = 0;
    double sigmaError = 0;
    double amplitude = 0;
    double fitLow = 0;
    double fitHigh = 0;
};

struct BinResult
{
    Long64_t entries = 0;
    std::array<TH1D *, kNMethods> hist = {};
    std::array<FitResult, kNMethods> fit;
};

FitResult fitCore(TH1D *hist, const char *name)
{
    FitResult output;
    if (!hist || hist->GetEntries() < 100)
        return output;

    const double probabilities[3] = {0.25, 0.50, 0.75};
    double quantiles[3] = {};
    hist->GetQuantiles(3, quantiles, probabilities);
    const double robustSigma = kIqrToSigma * (quantiles[2] - quantiles[0]);
    if (!(robustSigma > 0) || !std::isfinite(robustSigma))
        return output;

    output.fitLow = std::max(kResidualMin, quantiles[1] - 3 * robustSigma);
    output.fitHigh = std::min(kResidualMax, quantiles[1] + 3 * robustSigma);
    TF1 gaussian(name, "gaus", output.fitLow, output.fitHigh);
    gaussian.SetParameters(hist->GetMaximum(), quantiles[1], robustSigma);
    const TFitResultPtr result = hist->Fit(&gaussian, "SQNR");
    if (!result.Get() || result->Status() != 0)
        return output;

    output.mean = gaussian.GetParameter(1);
    output.meanError = gaussian.GetParError(1);
    output.sigma = std::abs(gaussian.GetParameter(2));
    output.sigmaError = gaussian.GetParError(2);
    output.amplitude = gaussian.GetParameter(0);
    output.valid = std::isfinite(output.mean) && std::isfinite(output.meanError) &&
                   std::isfinite(output.sigma) && std::isfinite(output.sigmaError) &&
                   output.sigma > 0;
    return output;
}

bool inBetaBin(double beta, int bin)
{
    const double low = kBetaMin + bin * kBetaWidth;
    const double high = low + kBetaWidth;
    if (bin + 1 == kNBetaBins)
        return beta >= low && beta <= 1.0;
    return beta >= low && beta < high;
}

void drawHeader(int charge, int bin, Long64_t entries)
{
    const double low = kBetaMin + bin * kBetaWidth;
    const double high = low + kBetaWidth;
    TLatex text;
    text.SetNDC();
    text.SetTextFont(42);
    text.SetTextSize(0.042);
    const char *right = (bin + 1 == kNBetaBins) ? "#leq" : "<";
    text.DrawLatex(0.12, 0.94,
                   Form("Z = %d, %.2f #leq #beta_{MC} %s %.2f", charge, low, right, high));
    text.SetTextSize(0.031);
    text.DrawLatex(0.82, 0.94, Form("N = %lld", entries));
}

void printPage(TCanvas &canvas, const TString &pdf, int page, int pages)
{
    if (page == 0)
        canvas.Print(pdf + "(");
    else if (page + 1 == pages)
        canvas.Print(pdf + ")");
    else
        canvas.Print(pdf);
}

void styleGraph(TGraphErrors *graph, int method)
{
    graph->SetLineColor(kColors[method]);
    graph->SetLineWidth(4);
    graph->SetMarkerStyle(1);
    graph->SetFillColorAlpha(kColors[method], 0.18);
}

void styleAxes(TH1D &frame)
{
    for (TAxis *axis : {frame.GetXaxis(), frame.GetYaxis()})
    {
        axis->SetLabelSize(0.038);
        axis->SetTitleSize(0.042);
        axis->CenterTitle();
    }
    frame.GetXaxis()->SetTitleOffset(1.15);
    frame.GetYaxis()->SetTitleOffset(1.25);
}

void drawZero(double xmin, double xmax)
{
    TLine *zero = new TLine(xmin, 0, xmax, 0);
    zero->SetLineColor(kGray + 1);
    zero->SetLineStyle(2);
    zero->SetLineWidth(2);
    zero->Draw();
}

} // namespace

void plotHybridComparison(const char *inputDir = "results/data_trained_comparison",
                          const char *outputDir = "results/data_trained_comparison/plots",
                          const char *baselineLabel = "Hybrid",
                          const char *candidateLabel = "Data-trained",
                          const char *baselineBranch = "hybridBeta",
                          const char *candidateBranch = "newBeta")
{
    kLabels[0] = baselineLabel;
    kLabels[1] = candidateLabel;
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gSystem->mkdir(outputDir, true);

    std::array<std::array<BinResult, kNBetaBins>, kNCharges> results;

    for (int iz = 0; iz < kNCharges; ++iz)
    {
        const int charge = kCharges[iz];
        TFile input(Form("%s/comparison_Z%d.root", inputDir, charge), "READ");
        TTree *tree = static_cast<TTree *>(input.Get("comparison"));
        if (!tree)
        {
            std::cerr << "Missing comparison tree for Z=" << charge << std::endl;
            continue;
        }

        double mcBeta = 0;
        double hybridBeta = 0;
        double newBeta = 0;
        tree->SetBranchAddress("mcBeta", &mcBeta);
        tree->SetBranchAddress(baselineBranch, &hybridBeta);
        tree->SetBranchAddress(candidateBranch, &newBeta);

        for (int bin = 0; bin < kNBetaBins; ++bin)
        {
            for (int method = 0; method < kNMethods; ++method)
            {
                results[iz][bin].hist[method] = new TH1D(
                    Form("residual_Z%d_B%d_M%d", charge, bin, method), "",
                    320, kResidualMin, kResidualMax);
                results[iz][bin].hist[method]->SetDirectory(nullptr);
            }
        }

        for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
        {
            tree->GetEntry(entry);
            if (!(mcBeta >= kBetaMin && mcBeta <= 1.0) ||
                !(hybridBeta > 0) || !(newBeta > 0) ||
                !std::isfinite(mcBeta) || !std::isfinite(hybridBeta) || !std::isfinite(newBeta))
                continue;
            int bin = static_cast<int>((mcBeta - kBetaMin) / kBetaWidth);
            bin = std::max(0, std::min(kNBetaBins - 1, bin));
            if (!inBetaBin(mcBeta, bin))
                continue;
            const double truth = 1.0 / mcBeta;
            const double residual[kNMethods] = {1.0 / hybridBeta - truth, 1.0 / newBeta - truth};
            for (int method = 0; method < kNMethods; ++method)
                results[iz][bin].hist[method]->Fill(residual[method]);
            ++results[iz][bin].entries;
        }

        for (int bin = 0; bin < kNBetaBins; ++bin)
            for (int method = 0; method < kNMethods; ++method)
                results[iz][bin].fit[method] = fitCore(
                    results[iz][bin].hist[method],
                    Form("fit_Z%d_B%d_M%d", charge, bin, method));
    }

    std::ofstream csv(Form("%s/per_bin_summary.csv", inputDir));
    csv << "Z,beta_low,beta_high,N,hybrid_mean,hybrid_mean_error,hybrid_sigma,hybrid_sigma_error,"
           "new_mean,new_mean_error,new_sigma,new_sigma_error,improvement_percent,improvement_error\n";
    csv << std::setprecision(10);
    for (int iz = 0; iz < kNCharges; ++iz)
        for (int bin = 0; bin < kNBetaBins; ++bin)
        {
            const FitResult &hybrid = results[iz][bin].fit[0];
            const FitResult &candidate = results[iz][bin].fit[1];
            const double improvement = 100.0 * (1.0 - candidate.sigma / hybrid.sigma);
            const double improvementError = 100.0 * std::sqrt(
                std::pow(candidate.sigmaError / hybrid.sigma, 2) +
                std::pow(candidate.sigma * hybrid.sigmaError / (hybrid.sigma * hybrid.sigma), 2));
            csv << kCharges[iz] << "," << kBetaMin + bin * kBetaWidth << ","
                << kBetaMin + (bin + 1) * kBetaWidth << "," << results[iz][bin].entries << ","
                << hybrid.mean << "," << hybrid.meanError << "," << hybrid.sigma << "," << hybrid.sigmaError << ","
                << candidate.mean << "," << candidate.meanError << "," << candidate.sigma << "," << candidate.sigmaError << ","
                << improvement << "," << improvementError << "\n";
        }
    csv.close();

    TFile detailFile(Form("%s/per_bin_details.root", inputDir), "RECREATE");

    // One page per beta bin: overlaid residual distributions and their Gaussian core fits.
    for (int iz = 0; iz < kNCharges; ++iz)
    {
        const int charge = kCharges[iz];
        TCanvas canvas(Form("bin_canvas_Z%d", charge), "", 1280, 720);
        const TString pdf = Form("%s/residual_fits_Z%d.pdf", outputDir, charge);
        for (int bin = 0; bin < kNBetaBins; ++bin)
        {
            canvas.Clear();
            canvas.SetLeftMargin(0.14);
            canvas.SetRightMargin(0.04);
            canvas.SetBottomMargin(0.16);
            canvas.SetTopMargin(0.12);

            std::array<TH1D *, kNMethods> shown = {};
            double maximum = 0;
            for (int method = 0; method < kNMethods; ++method)
            {
                TH1D *raw = results[iz][bin].hist[method];
                raw->Write();
                shown[method] = static_cast<TH1D *>(raw->Clone(Form("shown_Z%d_B%d_M%d", charge, bin, method)));
                shown[method]->SetDirectory(nullptr);
                const double integral = shown[method]->Integral(1, shown[method]->GetNbinsX());
                if (integral > 0)
                    shown[method]->Scale(1.0 / integral);
                shown[method]->SetLineColor(kColors[method]);
                shown[method]->SetLineWidth(3);
                shown[method]->SetFillStyle(0);
                maximum = std::max(maximum, shown[method]->GetMaximum());
            }

            TH1D frame(Form("frame_Z%d_B%d", charge, bin), "", 320, kResidualMin, kResidualMax);
            frame.GetXaxis()->SetTitle("1/#beta_{reco} - 1/#beta_{MC}");
            frame.GetYaxis()->SetTitle("Normalized entries");
            styleAxes(frame);
            frame.SetMinimum(0);
            frame.SetMaximum(maximum * 1.30);
            frame.Draw("AXIS");

            for (int method = 0; method < kNMethods; ++method)
            {
                shown[method]->Draw("HIST SAME");
                const FitResult &fit = results[iz][bin].fit[method];
                const double integral = results[iz][bin].hist[method]->Integral(1, 320);
                TF1 *curve = new TF1(Form("draw_fit_Z%d_B%d_M%d", charge, bin, method),
                                     "gaus", fit.fitLow, fit.fitHigh);
                curve->SetParameters(integral > 0 ? fit.amplitude / integral : 0, fit.mean, fit.sigma);
                curve->SetLineColor(kColors[method]);
                curve->SetLineStyle(2);
                curve->SetLineWidth(3);
                curve->Draw("SAME");
                TLatex annotation;
                annotation.SetNDC();
                annotation.SetTextFont(42);
                annotation.SetTextSize(0.030);
                annotation.SetTextColor(kColors[method]);
                const double top = 0.82 - 0.22 * method;
                annotation.DrawLatex(0.67, top, kLabels[method]);
                annotation.DrawLatex(0.67, top - 0.055,
                                     Form("#mu = %.5f #pm %.5f", fit.mean, fit.meanError));
                annotation.DrawLatex(0.67, top - 0.11,
                                     Form("#sigma = %.5f #pm %.5f", fit.sigma, fit.sigmaError));
                curve->Write();
            }
            drawHeader(charge, bin, results[iz][bin].entries);
            canvas.RedrawAxis();
            printPage(canvas, pdf, bin, kNBetaBins);

            for (TH1D *hist : shown)
                delete hist;
        }
    }

    detailFile.Close();

    // Summary line plots: sigma(beta), mean(beta), and relative sigma improvement.
    const char *summaryNames[3] = {"core_sigma_vs_beta", "core_mean_vs_beta", "relative_improvement_vs_beta"};
    for (int summary = 0; summary < 3; ++summary)
    {
        TCanvas canvas(Form("summary_canvas_%d", summary), "", 1280, 720);
        const TString pdf = Form("%s/%s.pdf", outputDir, summaryNames[summary]);
        for (int iz = 0; iz < kNCharges; ++iz)
        {
            canvas.Clear();
            canvas.SetLeftMargin(0.14);
            canvas.SetRightMargin(0.04);
            canvas.SetBottomMargin(0.16);
            canvas.SetTopMargin(0.12);

            double x[kNBetaBins] = {};
            double ex[kNBetaBins] = {};
            double y[kNMethods][kNBetaBins] = {};
            double ey[kNMethods][kNBetaBins] = {};
            double improvement[kNBetaBins] = {};
            double improvementError[kNBetaBins] = {};
            double ymin = std::numeric_limits<double>::max();
            double ymax = -std::numeric_limits<double>::max();

            for (int bin = 0; bin < kNBetaBins; ++bin)
            {
                x[bin] = kBetaMin + (bin + 0.5) * kBetaWidth;
                for (int method = 0; method < kNMethods; ++method)
                {
                    const FitResult &fit = results[iz][bin].fit[method];
                    y[method][bin] = (summary == 0) ? fit.sigma : fit.mean;
                    ey[method][bin] = (summary == 0) ? fit.sigmaError : fit.meanError;
                    if (summary < 2)
                    {
                        ymin = std::min(ymin, y[method][bin] - ey[method][bin]);
                        ymax = std::max(ymax, y[method][bin] + ey[method][bin]);
                    }
                }
                const FitResult &hybrid = results[iz][bin].fit[0];
                const FitResult &candidate = results[iz][bin].fit[1];
                improvement[bin] = 100.0 * (1.0 - candidate.sigma / hybrid.sigma);
                improvementError[bin] = 100.0 * std::sqrt(
                    std::pow(candidate.sigmaError / hybrid.sigma, 2) +
                    std::pow(candidate.sigma * hybrid.sigmaError / (hybrid.sigma * hybrid.sigma), 2));
                if (summary == 2)
                {
                    ymin = std::min(ymin, improvement[bin] - improvementError[bin]);
                    ymax = std::max(ymax, improvement[bin] + improvementError[bin]);
                }
            }

            const double span = std::max(1e-6, ymax - ymin);
            double frameLow = ymin - 0.18 * span;
            double frameHigh = ymax + 0.22 * span;
            if (summary == 0)
                frameLow = std::max(0.0, frameLow);
            if (summary == 1)
            {
                const double extent = std::max(std::abs(frameLow), std::abs(frameHigh));
                frameLow = -extent;
                frameHigh = extent;
            }
            if (summary == 2)
            {
                frameLow = std::min(0.0, frameLow);
                frameHigh = std::max(0.0, frameHigh);
            }

            TH1D frame(Form("summary_frame_%d_Z%d", summary, kCharges[iz]), "", 100, kBetaMin, 1.0);
            frame.GetXaxis()->SetTitle("#beta_{MC}");
            frame.GetYaxis()->SetTitle(summary == 0 ? "Gaussian core #sigma[1/#beta]" :
                                       summary == 1 ? "Gaussian core #mu[1/#beta]" :
                                                      Form("#sigma improvement over %s [%%]", baselineLabel));
            styleAxes(frame);
            frame.SetMinimum(frameLow);
            frame.SetMaximum(frameHigh);
            frame.Draw("AXIS");
            if (summary != 0)
                drawZero(kBetaMin, 1.0);

            TLegend legend(0.72, 0.73, 0.94, 0.87);
            legend.SetBorderSize(0);
            legend.SetFillStyle(0);
            legend.SetTextSize(0.032);

            std::vector<TGraphErrors *> graphCleanup;
            if (summary < 2)
            {
                for (int method = 0; method < kNMethods; ++method)
                {
                    TGraphErrors *band = new TGraphErrors(kNBetaBins, x, y[method], ex, ey[method]);
                    TGraphErrors *line = new TGraphErrors(kNBetaBins, x, y[method], ex, ey[method]);
                    styleGraph(band, method);
                    styleGraph(line, method);
                    band->Draw("3 SAME");
                    line->Draw("L SAME");
                    legend.AddEntry(line, kLabels[method], "l");
                    graphCleanup.push_back(band);
                    graphCleanup.push_back(line);
                }
            }
            else
            {
                TGraphErrors *band = new TGraphErrors(kNBetaBins, x, improvement, ex, improvementError);
                TGraphErrors *line = new TGraphErrors(kNBetaBins, x, improvement, ex, improvementError);
                styleGraph(band, 0);
                styleGraph(line, 0);
                band->Draw("3 SAME");
                line->Draw("L SAME");
                legend.AddEntry(line, candidateLabel, "l");
                graphCleanup.push_back(band);
                graphCleanup.push_back(line);
            }
            legend.Draw();

            TLatex title;
            title.SetNDC();
            title.SetTextFont(42);
            title.SetTextSize(0.043);
            title.DrawLatex(0.12, 0.94, Form("Z = %d", kCharges[iz]));
            canvas.RedrawAxis();
            printPage(canvas, pdf, iz, kNCharges);
            for (TGraphErrors *graph : graphCleanup)
                delete graph;
        }
    }

    std::cout << "Wrote per-bin fits and Hybrid comparison summaries to " << outputDir << std::endl;
}
