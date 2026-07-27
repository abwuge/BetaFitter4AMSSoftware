#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>
#include <TLegend.h>
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
#include <string>
#include <vector>

namespace
{
constexpr double kBetaMin = 0.30;
constexpr double kBetaMax = 1.00;
constexpr double kBetaWidth = 0.05;
constexpr int kNBins = 14;

double quantile(std::vector<double> values, double probability)
{
    if (values.empty())
        return 0;
    std::sort(values.begin(), values.end());
    const double position = probability * (values.size() - 1);
    const size_t lower = static_cast<size_t>(position);
    const size_t upper = std::min(lower + 1, values.size() - 1);
    const double fraction = position - lower;
    return values[lower] * (1 - fraction) + values[upper] * fraction;
}

struct BinSummary
{
    int total = 0;
    int boundary = 0;
    std::vector<double> interior;
};

struct ChargeSummary
{
    int charge = 0;
    std::array<BinSummary, kNBins> bins;
};

ChargeSummary analyze(const char *path, int charge)
{
    ChargeSummary result;
    result.charge = charge;
    TFile file(path, "READ");
    TTree *tree = static_cast<TTree *>(file.Get("scaleTree"));
    if (!tree)
    {
        std::cerr << "Error: missing scaleTree in " << path << std::endl;
        return result;
    }

    float beta = 0;
    float zeta = 0;
    tree->SetBranchAddress("mcBeta", &beta);
    tree->SetBranchAddress("energyLossScale", &zeta);
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
    {
        tree->GetEntry(entry);
        if (!std::isfinite(beta) || !std::isfinite(zeta) ||
            beta < kBetaMin || beta > kBetaMax)
            continue;
        int bin = std::min(kNBins - 1, static_cast<int>((beta - kBetaMin) / kBetaWidth));
        BinSummary &summary = result.bins[bin];
        ++summary.total;
        const bool atBoundary = zeta <= -7.999 || zeta >= 11.999;
        summary.boundary += atBoundary;
        if (!atBoundary)
            summary.interior.push_back(zeta);
    }
    return result;
}
}

void analyzeCenterZeta(
    const char *z2File = "results/validation_center_zeta/zeta_Z2.root",
    const char *z6File = "results/validation_center_zeta/zeta_Z6.root",
    const char *z8File = "results/validation_center_zeta/zeta_Z8.root",
    const char *outputDirectory = "results/validation_center_zeta")
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    gSystem->mkdir(outputDirectory, true);

    const std::array<const char *, 3> paths = {z2File, z6File, z8File};
    const std::array<int, 3> charges = {2, 6, 8};
    std::array<ChargeSummary, 3> summaries;
    for (int i = 0; i < 3; ++i)
        summaries[i] = analyze(paths[i], charges[i]);

    std::ofstream table((std::string(outputDirectory) + "/zeta_by_beta.tsv").c_str());
    table << "Z\tbeta_low\tbeta_high\ttotal\tinterior\tboundary_fraction\tq16\tmedian\tq84\n";
    table << std::setprecision(7);

    const int color = TColor::GetColor("#45364B");
    const int boundaryColor = TColor::GetColor("#E55812");
    TCanvas canvas("centerZeta", "", 1600, 900);
    canvas.Divide(3, 2, 0.002, 0.002);
    for (int chargeIndex = 0; chargeIndex < 3; ++chargeIndex)
    {
        const ChargeSummary &summary = summaries[chargeIndex];
        std::vector<double> beta, median, lowError, highError, betaError;
        std::vector<double> boundaryBeta, boundaryFraction;
        for (int bin = 0; bin < kNBins; ++bin)
        {
            const BinSummary &item = summary.bins[bin];
            const double low = kBetaMin + bin * kBetaWidth;
            const double high = low + kBetaWidth;
            const double fraction = item.total > 0 ? 1.0 * item.boundary / item.total : 0;
            const double q16 = quantile(item.interior, 0.16);
            const double q50 = quantile(item.interior, 0.50);
            const double q84 = quantile(item.interior, 0.84);
            table << summary.charge << '\t' << low << '\t' << high << '\t'
                  << item.total << '\t' << item.interior.size() << '\t' << fraction << '\t'
                  << q16 << '\t' << q50 << '\t' << q84 << '\n';
            if (item.interior.size() >= 30)
            {
                beta.push_back((low + high) / 2);
                betaError.push_back(kBetaWidth / 2);
                median.push_back(q50);
                lowError.push_back(q50 - q16);
                highError.push_back(q84 - q50);
            }
            if (item.total > 0)
            {
                boundaryBeta.push_back((low + high) / 2);
                boundaryFraction.push_back(fraction);
            }
        }

        canvas.cd(chargeIndex + 1);
        gPad->SetLeftMargin(0.17);
        gPad->SetRightMargin(0.035);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.12);
        gPad->SetGridx();
        gPad->SetGridy();
        TH1F *frame = gPad->DrawFrame(kBetaMin, -2, kBetaMax, 7);
        frame->SetTitle(Form("Z = %d", summary.charge));
        frame->GetXaxis()->SetTitle("#beta_{MC} at AMS center");
        frame->GetYaxis()->SetTitle("Fitted #zeta (16/50/84%)");
        frame->GetXaxis()->SetTitleSize(0.042);
        frame->GetYaxis()->SetTitleSize(0.042);
        frame->GetXaxis()->SetLabelSize(0.034);
        frame->GetYaxis()->SetLabelSize(0.034);
        frame->GetYaxis()->SetTitleOffset(1.55);
        TGraphAsymmErrors *graph = new TGraphAsymmErrors(
            beta.size(), beta.data(), median.data(), betaError.data(), betaError.data(),
            lowError.data(), highError.data());
        graph->SetLineColor(color);
        graph->SetMarkerColor(color);
        graph->SetFillColorAlpha(color, 0.18);
        graph->SetMarkerStyle(20);
        graph->SetLineWidth(2);
        graph->Draw("3LP SAME");
        TLatex chargeLabel;
        chargeLabel.SetNDC(true);
        chargeLabel.SetTextSize(0.038);
        chargeLabel.DrawLatex(0.20, 0.91, Form("Z = %d", summary.charge));

        canvas.cd(chargeIndex + 4);
        gPad->SetLeftMargin(0.17);
        gPad->SetRightMargin(0.035);
        gPad->SetBottomMargin(0.15);
        gPad->SetTopMargin(0.06);
        gPad->SetGridx();
        gPad->SetGridy();
        TH1F *boundaryFrame = gPad->DrawFrame(kBetaMin, 0, kBetaMax, 1.03);
        boundaryFrame->GetXaxis()->SetTitle("#beta_{MC} at AMS center");
        boundaryFrame->GetYaxis()->SetTitle("Boundary fraction");
        boundaryFrame->GetXaxis()->SetTitleSize(0.042);
        boundaryFrame->GetYaxis()->SetTitleSize(0.042);
        boundaryFrame->GetXaxis()->SetLabelSize(0.034);
        boundaryFrame->GetYaxis()->SetLabelSize(0.034);
        boundaryFrame->GetYaxis()->SetTitleOffset(1.55);
        TGraph *boundaryGraph = new TGraph(boundaryBeta.size(), boundaryBeta.data(), boundaryFraction.data());
        boundaryGraph->SetLineColor(boundaryColor);
        boundaryGraph->SetMarkerColor(boundaryColor);
        boundaryGraph->SetMarkerStyle(21);
        boundaryGraph->SetLineWidth(2);
        boundaryGraph->Draw("LP SAME");
    }

    const std::string pdf = std::string(outputDirectory) + "/center_zeta_diagnostic.pdf";
    const std::string png = std::string(outputDirectory) + "/center_zeta_diagnostic.png";
    canvas.Print(pdf.c_str());
    canvas.Print(png.c_str());
    std::cout << "Wrote " << table.tellp() << " bytes to zeta_by_beta.tsv" << std::endl;
}
