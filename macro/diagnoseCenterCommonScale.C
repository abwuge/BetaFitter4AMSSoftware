#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Summary {
    long long count = 0;
    double mean = 0;
    double q16 = 0;
    double median = 0;
    double q84 = 0;
};

double quantile(const std::vector<double> &sorted, double probability)
{
    if (sorted.empty())
        return 0;
    const double position = probability * (sorted.size() - 1);
    const size_t low = static_cast<size_t>(position);
    const size_t high = std::min(low + 1, sorted.size() - 1);
    const double fraction = position - low;
    return sorted[low] * (1.0 - fraction) + sorted[high] * fraction;
}

Summary summarize(std::vector<double> values)
{
    Summary result;
    result.count = values.size();
    if (values.empty())
        return result;
    double sum = 0;
    for (double value : values)
        sum += value;
    result.mean = sum / values.size();
    std::sort(values.begin(), values.end());
    result.q16 = quantile(values, 0.16);
    result.median = quantile(values, 0.50);
    result.q84 = quantile(values, 0.84);
    return result;
}

void writeSummary(std::ofstream &output, const char *sample, int charge,
                  const char *selection, const char *metric,
                  const Summary &summary)
{
    output << sample << '\t' << charge << '\t' << selection << '\t'
           << metric << '\t' << summary.count << '\t' << summary.mean
           << '\t' << summary.q16 << '\t' << summary.median << '\t'
           << summary.q84 << '\t'
           << 0.5 * (summary.q84 - summary.q16) << '\n';
}

double generatorMass(int charge)
{
    if (charge == 8)
        return 14.899169;
    return 2.0 * charge * 0.9314941;
}

void analyzeFile(std::ofstream &output, const char *sample, int requestedCharge,
                 const char *fileName, double rigidityLow, double rigidityHigh,
                 double betaLow, double betaHigh, const char *selectionName,
                 bool isMC, bool selectOnMCBeta = true)
{
    TFile file(fileName, "READ");
    TTree *tree = static_cast<TTree *>(file.Get("betaTree"));
    if (!tree) {
        std::cerr << "Missing betaTree in " << fileName << std::endl;
        return;
    }

    double mcBeta = 0;
    double rigidity = 0;
    double linearBeta = 0;
    double nonlinearBeta = 0;
    int charge = 0;
    if (isMC)
        tree->SetBranchAddress("mcBeta", &mcBeta);
    tree->SetBranchAddress("innerRigidity", &rigidity);
    tree->SetBranchAddress("linearBeta", &linearBeta);
    tree->SetBranchAddress("nonlinearBeta", &nonlinearBeta);
    tree->SetBranchAddress("Z", &charge);

    std::vector<double> beta;
    std::vector<double> deltaBeta;
    std::vector<double> deltaInverseBeta;
    std::vector<double> relativeMassChange;
    std::vector<double> predictedRelativeMassChange;
    std::vector<double> linearMassBias;
    std::vector<double> nonlinearMassBias;
    std::vector<double> rigidityBias;
    std::vector<double> linearInverseBetaResidual;
    std::vector<double> nonlinearInverseBetaResidual;

    const double mass = generatorMass(requestedCharge);
    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) {
        tree->GetEntry(entry);
        if (charge != requestedCharge || !std::isfinite(rigidity) ||
            rigidity <= rigidityLow || rigidity >= rigidityHigh ||
            !std::isfinite(linearBeta) || !std::isfinite(nonlinearBeta) ||
            linearBeta <= 0 || linearBeta >= 1 ||
            nonlinearBeta <= 0 || nonlinearBeta >= 1)
            continue;
        if (isMC && (!std::isfinite(mcBeta) || mcBeta <= 0 || mcBeta >= 1))
            continue;
        const double selectionBeta =
            isMC && selectOnMCBeta ? mcBeta : linearBeta;
        if (selectionBeta <= betaLow || selectionBeta >= betaHigh)
            continue;

        const double momentum = std::abs(rigidity) * requestedCharge;
        const double linearMass = momentum *
            std::sqrt(1.0 / (linearBeta * linearBeta) - 1.0);
        const double nonlinearMass = momentum *
            std::sqrt(1.0 / (nonlinearBeta * nonlinearBeta) - 1.0);
        const double betaMid = 0.5 * (linearBeta + nonlinearBeta);
        const double betaDifference = nonlinearBeta - linearBeta;

        beta.push_back(linearBeta);
        deltaBeta.push_back(betaDifference);
        deltaInverseBeta.push_back(1.0 / nonlinearBeta - 1.0 / linearBeta);
        relativeMassChange.push_back(nonlinearMass / linearMass - 1.0);
        predictedRelativeMassChange.push_back(
            -betaDifference / (betaMid * (1.0 - betaMid * betaMid)));

        if (isMC) {
            const double trueRigidity =
                (mass / requestedCharge) * mcBeta /
                std::sqrt(1.0 - mcBeta * mcBeta);
            rigidityBias.push_back(rigidity / trueRigidity - 1.0);
            linearMassBias.push_back(linearMass / mass - 1.0);
            nonlinearMassBias.push_back(nonlinearMass / mass - 1.0);
            linearInverseBetaResidual.push_back(
                1.0 / linearBeta - 1.0 / mcBeta);
            nonlinearInverseBetaResidual.push_back(
                1.0 / nonlinearBeta - 1.0 / mcBeta);
        }
    }

    writeSummary(output, sample, requestedCharge, selectionName,
                 "linear_beta", summarize(beta));
    writeSummary(output, sample, requestedCharge, selectionName,
                 "nonlinear_minus_linear_beta", summarize(deltaBeta));
    writeSummary(output, sample, requestedCharge, selectionName,
                 "nonlinear_minus_linear_invbeta",
                 summarize(deltaInverseBeta));
    writeSummary(output, sample, requestedCharge, selectionName,
                 "relative_mass_change", summarize(relativeMassChange));
    writeSummary(output, sample, requestedCharge, selectionName,
                 "linearized_relative_mass_change",
                 summarize(predictedRelativeMassChange));
    if (isMC) {
        writeSummary(output, sample, requestedCharge, selectionName,
                     "inner_rigidity_bias", summarize(rigidityBias));
        writeSummary(output, sample, requestedCharge, selectionName,
                     "linear_mass_bias", summarize(linearMassBias));
        writeSummary(output, sample, requestedCharge, selectionName,
                     "nonlinear_mass_bias", summarize(nonlinearMassBias));
        writeSummary(output, sample, requestedCharge, selectionName,
                     "linear_invbeta_residual",
                     summarize(linearInverseBetaResidual));
        writeSummary(output, sample, requestedCharge, selectionName,
                     "nonlinear_invbeta_residual",
                     summarize(nonlinearInverseBetaResidual));
    }
}

int chargeIndex(int charge)
{
    return charge == 2 ? 0 : (charge == 6 ? 1 : (charge == 8 ? 2 : -1));
}

TGraph *makeGraph(const double *values, double xOffset, Color_t color,
                  Style_t marker)
{
    auto *graph = new TGraph(3);
    for (int index = 0; index < 3; ++index)
        graph->SetPoint(index, index + 1.0 + xOffset, values[index]);
    graph->SetLineColor(color);
    graph->SetLineWidth(3);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(marker);
    graph->SetMarkerSize(1.6);
    return graph;
}

void styleFrame(TH1D *frame, const char *title, const char *yTitle,
                double minimum, double maximum)
{
    frame->SetStats(false);
    frame->SetTitle(title);
    frame->SetMinimum(minimum);
    frame->SetMaximum(maximum);
    frame->GetXaxis()->SetBinLabel(1, "Z=2");
    frame->GetXaxis()->SetBinLabel(2, "Z=6");
    frame->GetXaxis()->SetBinLabel(3, "Z=8");
    frame->GetYaxis()->SetTitle(yTitle);
    frame->GetXaxis()->SetLabelSize(0.050);
    frame->GetYaxis()->SetLabelSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.048);
    frame->GetYaxis()->SetTitleOffset(1.25);
}

void drawPanelTitle(const char *title)
{
    TLatex text;
    text.SetNDC();
    text.SetTextAlign(22);
    text.SetTextSize(0.040);
    text.DrawLatex(0.56, 0.955, title);
}

void plotResponsibility(const char *diagnosticName,
                        const char *massSummaryName,
                        const char *outputDirectory)
{
    const double nan = std::numeric_limits<double>::quiet_NaN();
    double dataLinearMass[3] = {nan, nan, nan};
    double dataNonlinearMass[3] = {nan, nan, nan};
    double mcRigidity[3] = {nan, nan, nan};
    double mcLinearMass[3] = {nan, nan, nan};
    double mcNonlinearMass[3] = {nan, nan, nan};
    double mcLinearInvBeta[3] = {nan, nan, nan};
    double mcNonlinearInvBeta[3] = {nan, nan, nan};

    std::ifstream massInput(massSummaryName);
    std::string line;
    std::getline(massInput, line);
    while (std::getline(massInput, line)) {
        std::istringstream row(line);
        std::string species;
        std::string method;
        int charge = 0;
        double knownMass = 0;
        long long selected = 0;
        int valid = 0;
        double mean = 0;
        double meanError = 0;
        double sigma = 0;
        double sigmaError = 0;
        double bias = 0;
        double resolution = 0;
        row >> species >> charge >> knownMass >> method >> selected >> valid
            >> mean >> meanError >> sigma >> sigmaError >> bias >> resolution;
        const int index = chargeIndex(charge);
        if (index < 0 || !valid)
            continue;
        if (method == "linear")
            dataLinearMass[index] = bias;
        else if (method == "nonlinear")
            dataNonlinearMass[index] = bias;
    }

    std::ifstream diagnosticInput(diagnosticName);
    std::getline(diagnosticInput, line);
    while (std::getline(diagnosticInput, line)) {
        std::istringstream row(line);
        std::string sample;
        std::string selection;
        std::string metric;
        int charge = 0;
        long long count = 0;
        double mean = 0;
        double q16 = 0;
        double median = 0;
        double q84 = 0;
        double sigma68 = 0;
        row >> sample >> charge >> selection >> metric >> count >> mean >> q16
            >> median >> q84 >> sigma68;
        const int index = chargeIndex(charge);
        if (sample != "mc" || selection != "0.8<R<1.4" || index < 0)
            continue;
        if (metric == "inner_rigidity_bias")
            mcRigidity[index] = 100.0 * median;
        else if (metric == "linear_mass_bias")
            mcLinearMass[index] = 100.0 * median;
        else if (metric == "nonlinear_mass_bias")
            mcNonlinearMass[index] = 100.0 * median;
        else if (metric == "linear_invbeta_residual")
            mcLinearInvBeta[index] = median;
        else if (metric == "nonlinear_invbeta_residual")
            mcNonlinearInvBeta[index] = median;
    }

    const double *allValues[] = {
        dataLinearMass, dataNonlinearMass, mcRigidity, mcLinearMass,
        mcNonlinearMass, mcLinearInvBeta, mcNonlinearInvBeta};
    for (const double *values : allValues)
        for (int index = 0; index < 3; ++index)
            if (!std::isfinite(values[index])) {
                std::cerr << "Missing responsibility input for Z index "
                          << index << std::endl;
                return;
            }

    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    gStyle->SetOptStat(0);
    const Color_t linearColor = TColor::GetColor("#45364B");
    const Color_t nonlinearColor = TColor::GetColor("#788AA3");
    const Color_t rigidityColor = TColor::GetColor("#447604");

    TCanvas canvas("responsibilityCanvas", "", 1920, 640);
    canvas.Divide(3, 1, 0.001, 0.001);

    canvas.cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.14);
    gPad->SetGridy();
    TH1D dataFrame("dataFrame", "", 3, 0.5, 3.5);
    styleFrame(&dataFrame, "Flight data mass-peak fit, 0.8<R<1.4 GV",
               "Peak bias [%]", -9.0, 0.5);
    dataFrame.Draw();
    drawPanelTitle("Flight data mass peak, 0.8<R<1.4 GV");
    TLine dataZero(0.5, 0, 3.5, 0);
    dataZero.SetLineStyle(7);
    dataZero.SetLineWidth(2);
    dataZero.SetLineColor(rigidityColor);
    dataZero.Draw("SAME");
    TGraph *dataLinear = makeGraph(dataLinearMass, -0.08, linearColor, 20);
    TGraph *dataNonlinear =
        makeGraph(dataNonlinearMass, 0.08, nonlinearColor, 21);
    dataLinear->Draw("LP SAME");
    dataNonlinear->Draw("LP SAME");
    TLegend dataLegend(0.18, 0.75, 0.54, 0.88);
    dataLegend.SetBorderSize(0);
    dataLegend.SetFillStyle(0);
    dataLegend.AddEntry(dataLinear, "Linear #beta", "lp");
    dataLegend.AddEntry(dataNonlinear, "Nonlinear #beta", "lp");
    dataLegend.Draw();

    canvas.cd(2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.14);
    gPad->SetGridy();
    TH1D mcFrame("mcFrame", "", 3, 0.5, 3.5);
    styleFrame(&mcFrame, "MC event-common medians, 0.8<R<1.4 GV",
               "Relative bias [%]", -8.0, 1.0);
    mcFrame.Draw();
    drawPanelTitle("MC event-common medians, 0.8<R<1.4 GV");
    TLine mcZero(0.5, 0, 3.5, 0);
    mcZero.SetLineStyle(7);
    mcZero.SetLineWidth(2);
    mcZero.SetLineColor(rigidityColor);
    mcZero.Draw("SAME");
    TGraph *rigidity = makeGraph(mcRigidity, -0.14, rigidityColor, 22);
    TGraph *linearMass = makeGraph(mcLinearMass, 0.0, linearColor, 20);
    TGraph *nonlinearMass =
        makeGraph(mcNonlinearMass, 0.14, nonlinearColor, 21);
    rigidity->Draw("LP SAME");
    linearMass->Draw("LP SAME");
    nonlinearMass->Draw("LP SAME");
    TLegend mcLegend(0.18, 0.70, 0.67, 0.88);
    mcLegend.SetBorderSize(0);
    mcLegend.SetFillStyle(0);
    mcLegend.AddEntry(rigidity, "Inner rigidity", "lp");
    mcLegend.AddEntry(linearMass, "Mass with linear #beta", "lp");
    mcLegend.AddEntry(nonlinearMass, "Mass with nonlinear #beta", "lp");
    mcLegend.Draw();

    canvas.cd(3);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.04);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.14);
    gPad->SetGridy();
    TH1D betaFrame("betaFrame", "", 3, 0.5, 3.5);
    styleFrame(&betaFrame, "MC event-common beta residual, 0.8<R<1.4 GV",
               "#Delta(1/#beta)", -0.001, 0.010);
    betaFrame.Draw();
    drawPanelTitle("MC event-common beta residual, 0.8<R<1.4 GV");
    TLine betaZero(0.5, 0, 3.5, 0);
    betaZero.SetLineStyle(7);
    betaZero.SetLineWidth(2);
    betaZero.SetLineColor(rigidityColor);
    betaZero.Draw("SAME");
    TGraph *linearInv =
        makeGraph(mcLinearInvBeta, -0.08, linearColor, 20);
    TGraph *nonlinearInv =
        makeGraph(mcNonlinearInvBeta, 0.08, nonlinearColor, 21);
    linearInv->Draw("LP SAME");
    nonlinearInv->Draw("LP SAME");
    TLegend betaLegend(0.18, 0.18, 0.58, 0.32);
    betaLegend.SetBorderSize(0);
    betaLegend.SetFillStyle(0);
    betaLegend.AddEntry(linearInv, "Linear #beta", "lp");
    betaLegend.AddEntry(nonlinearInv, "Nonlinear #beta", "lp");
    betaLegend.Draw();

    const TString outputBase = TString(outputDirectory) +
                               "/mass_offset_responsibility";
    canvas.SaveAs(outputBase + ".pdf");
    canvas.SaveAs(outputBase + ".png");
    std::cout << "Wrote " << outputBase << ".{pdf,png}" << std::endl;
}

} // namespace

void diagnoseCenterCommonScale(
    const char *dataFile =
        "results/data_p000-009_center_common_20260805_z268.root",
    const char *mcDirectory =
        "results/energy_loss_scale_20260805/beta_at_common_scale",
    const char *outputName =
        "macro/results/data_p000-009_center_common_20260805_z268_lowR/"
        "center_common_mc_diagnostics.tsv",
    const char *massSummaryName =
        "macro/results/data_p000-009_center_common_20260805_z268_lowR/"
        "rigidity_mass_peaks.tsv",
    double rigidityLow = 0.8,
    double rigidityHigh = 1.4)
{
    std::ofstream output(outputName);
    output << std::setprecision(10);
    output << "sample\tZ\tselection\tmetric\tcount\tmean\tq16"
              "\tmedian\tq84\tsigma68\n";
    for (int charge : {2, 6, 8}) {
        const double betaLow = charge == 2 ? 0.35 : (charge == 6 ? 0.40 : 0.45);
        const double betaHigh = charge == 2 ? 0.50 : 0.60;
        analyzeFile(output, "data", charge, dataFile, rigidityLow,
                    rigidityHigh, 0.0, 1.0, "0.8<R<1.4", false);
        const std::string mcFile =
            std::string(mcDirectory) + Form("/beta_Z%d.root", charge);
        analyzeFile(output, "mc", charge, mcFile.c_str(), rigidityLow,
                    rigidityHigh, 0.0, 1.0, "0.8<R<1.4", true);
        const std::string betaSelection =
            Form("sensitive_beta_%.2f_%.2f", betaLow, betaHigh);
        analyzeFile(output, "data", charge, dataFile, 0.0, 1.0e6,
                    betaLow, betaHigh, betaSelection.c_str(), false);
        analyzeFile(output, "mc", charge, mcFile.c_str(), 0.0, 1.0e6,
                    betaLow, betaHigh, betaSelection.c_str(), true);
        const std::string matchedSelection =
            Form("matched_linear_beta_%.2f_%.2f", betaLow, betaHigh);
        analyzeFile(output, "mc", charge, mcFile.c_str(), 0.0, 1.0e6,
                    betaLow, betaHigh, matchedSelection.c_str(), true, false);
    }
    output.close();
    const TString outputDirectory = gSystem->DirName(outputName);
    plotResponsibility(outputName, massSummaryName, outputDirectory.Data());
    std::cout << "Wrote " << outputName << std::endl;
}
