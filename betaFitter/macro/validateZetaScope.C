#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace
{
struct Summary
{
    double mean = 0;
    double rms = 0;
    double median = 0;
    double q16 = 0;
    double q84 = 0;
    double q95Abs = 0;
    double q99Abs = 0;
    double maxAbs = 0;
};

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

Summary summarize(const std::vector<double> &values)
{
    Summary result;
    if (values.empty())
        return result;

    result.mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double sumSquared = 0;
    std::vector<double> absoluteValues;
    absoluteValues.reserve(values.size());
    for (double value : values)
    {
        sumSquared += (value - result.mean) * (value - result.mean);
        absoluteValues.push_back(std::abs(value));
        result.maxAbs = std::max(result.maxAbs, std::abs(value));
    }
    result.rms = std::sqrt(sumSquared / values.size());
    result.median = quantile(values, 0.5);
    result.q16 = quantile(values, 0.16);
    result.q84 = quantile(values, 0.84);
    result.q95Abs = quantile(absoluteValues, 0.95);
    result.q99Abs = quantile(absoluteValues, 0.99);
    return result;
}

double correlation(const std::vector<double> &x, const std::vector<double> &y)
{
    if (x.size() != y.size() || x.empty())
        return 0;
    const Summary sx = summarize(x);
    const Summary sy = summarize(y);
    double covariance = 0;
    for (size_t i = 0; i < x.size(); ++i)
        covariance += (x[i] - sx.mean) * (y[i] - sy.mean);
    covariance /= x.size();
    return covariance / (sx.rms * sy.rms);
}

void writeSummary(std::ofstream &output, const std::string &name, const Summary &summary)
{
    output << name << ".mean=" << summary.mean << '\n'
           << name << ".rms=" << summary.rms << '\n'
           << name << ".median=" << summary.median << '\n'
           << name << ".q16=" << summary.q16 << '\n'
           << name << ".q84=" << summary.q84 << '\n'
           << name << ".q95_abs=" << summary.q95Abs << '\n'
           << name << ".q99_abs=" << summary.q99Abs << '\n'
           << name << ".max_abs=" << summary.maxAbs << '\n';
}
}

void validateZetaScope(
    const char *betaAllPath = "results/validation_zeta_scope/beta_all.root",
    const char *betaS2Path = "results/validation_zeta_scope/beta_s2.root",
    const char *zetaAllPath = "results/validation_zeta_scope/zeta_all.root",
    const char *zetaS2Path = "results/validation_zeta_scope/zeta_s2.root",
    const char *outputDirectory = "results/validation_zeta_scope",
    const char *alternativeLabel = "#zeta on S2 only")
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gROOT->ForceStyle();
    TH1::AddDirectory(false);
    gSystem->mkdir(outputDirectory, true);

    TFile betaAllFile(betaAllPath, "READ");
    TFile betaS2File(betaS2Path, "READ");
    TFile zetaAllFile(zetaAllPath, "READ");
    TFile zetaS2File(zetaS2Path, "READ");
    TTree *betaAllTree = static_cast<TTree *>(betaAllFile.Get("betaTree"));
    TTree *betaS2Tree = static_cast<TTree *>(betaS2File.Get("betaTree"));
    TTree *zetaAllTree = static_cast<TTree *>(zetaAllFile.Get("scaleTree"));
    TTree *zetaS2Tree = static_cast<TTree *>(zetaS2File.Get("scaleTree"));
    if (!betaAllTree || !betaS2Tree || !zetaAllTree || !zetaS2Tree)
    {
        std::cerr << "Error: one or more validation trees are missing" << std::endl;
        return;
    }
    if (betaAllTree->GetEntries() != betaS2Tree->GetEntries() ||
        zetaAllTree->GetEntries() != zetaS2Tree->GetEntries())
    {
        std::cerr << "Error: paired files have different entry counts" << std::endl;
        return;
    }

    double mcBetaAll = 0, mcBetaS2 = 0, nonlinearBetaAll = 0, nonlinearBetaS2 = 0;
    betaAllTree->SetBranchAddress("mcBeta", &mcBetaAll);
    betaAllTree->SetBranchAddress("nonlinearBeta", &nonlinearBetaAll);
    betaS2Tree->SetBranchAddress("mcBeta", &mcBetaS2);
    betaS2Tree->SetBranchAddress("nonlinearBeta", &nonlinearBetaS2);

    std::vector<double> betaValues, residualAll, residualS2, betaDifference;
    for (Long64_t entry = 0; entry < betaAllTree->GetEntries(); ++entry)
    {
        betaAllTree->GetEntry(entry);
        betaS2Tree->GetEntry(entry);
        if (std::abs(mcBetaAll - mcBetaS2) > 1e-12)
        {
            std::cerr << "Error: beta files are not entry-aligned at " << entry << std::endl;
            return;
        }
        if (!std::isfinite(mcBetaAll) || !std::isfinite(nonlinearBetaAll) || !std::isfinite(nonlinearBetaS2) ||
            mcBetaAll <= 0 || nonlinearBetaAll <= 0 || nonlinearBetaS2 <= 0)
            continue;
        betaValues.push_back(mcBetaAll);
        residualAll.push_back(1.0 / nonlinearBetaAll - 1.0 / mcBetaAll);
        residualS2.push_back(1.0 / nonlinearBetaS2 - 1.0 / mcBetaAll);
        betaDifference.push_back(1.0 / nonlinearBetaS2 - 1.0 / nonlinearBetaAll);
    }

    float zetaAllValue = 0, zetaS2Value = 0, zetaMcBetaAll = 0, zetaMcBetaS2 = 0;
    zetaAllTree->SetBranchAddress("energyLossScale", &zetaAllValue);
    zetaAllTree->SetBranchAddress("mcBeta", &zetaMcBetaAll);
    zetaS2Tree->SetBranchAddress("energyLossScale", &zetaS2Value);
    zetaS2Tree->SetBranchAddress("mcBeta", &zetaMcBetaS2);

    std::vector<double> zetaAll, zetaS2, zetaDifference;
    for (Long64_t entry = 0; entry < zetaAllTree->GetEntries(); ++entry)
    {
        zetaAllTree->GetEntry(entry);
        zetaS2Tree->GetEntry(entry);
        if (std::abs(zetaMcBetaAll - zetaMcBetaS2) > 1e-6)
        {
            std::cerr << "Error: zeta files are not entry-aligned at " << entry << std::endl;
            return;
        }
        if (zetaMcBetaAll <= 0 || zetaMcBetaAll >= 0.9 ||
            !std::isfinite(zetaAllValue) || !std::isfinite(zetaS2Value))
            continue;
        zetaAll.push_back(zetaAllValue);
        zetaS2.push_back(zetaS2Value);
        zetaDifference.push_back(zetaS2Value - zetaAllValue);
    }
    if (betaValues.empty() || zetaAll.empty())
    {
        std::cerr << "Error: no valid MC entries are available for validation" << std::endl;
        return;
    }

    const Summary residualAllSummary = summarize(residualAll);
    const Summary residualS2Summary = summarize(residualS2);
    const Summary betaDifferenceSummary = summarize(betaDifference);
    const Summary zetaAllSummary = summarize(zetaAll);
    const Summary zetaS2Summary = summarize(zetaS2);
    const Summary zetaDifferenceSummary = summarize(zetaDifference);

    const int nBetaBins = 10;
    const double betaMin = *std::min_element(betaValues.begin(), betaValues.end());
    const double betaMax = *std::max_element(betaValues.begin(), betaValues.end());
    double maxBinnedMeanDifference = 0;
    double maxBinnedRmsRelativeDifference = 0;
    std::vector<int> binCounts(nBetaBins, 0);
    std::vector<double> binLows(nBetaBins, 0), binHighs(nBetaBins, 0);
    std::vector<double> binAllMeans(nBetaBins, 0), binS2Means(nBetaBins, 0);
    std::vector<double> binAllRms(nBetaBins, 0), binS2Rms(nBetaBins, 0);
    for (int bin = 0; bin < nBetaBins; ++bin)
    {
        const double low = betaMin + (betaMax - betaMin) * bin / nBetaBins;
        const double high = betaMin + (betaMax - betaMin) * (bin + 1) / nBetaBins;
        std::vector<double> binAll, binS2;
        for (size_t i = 0; i < betaValues.size(); ++i)
        {
            if (betaValues[i] >= low && (betaValues[i] < high || (bin == nBetaBins - 1 && betaValues[i] <= high)))
            {
                binAll.push_back(residualAll[i]);
                binS2.push_back(residualS2[i]);
            }
        }
        if (binAll.empty())
            continue;
        const Summary allSummary = summarize(binAll);
        const Summary s2Summary = summarize(binS2);
        binCounts[bin] = binAll.size();
        binLows[bin] = low;
        binHighs[bin] = high;
        binAllMeans[bin] = allSummary.mean;
        binS2Means[bin] = s2Summary.mean;
        binAllRms[bin] = allSummary.rms;
        binS2Rms[bin] = s2Summary.rms;
        maxBinnedMeanDifference = std::max(maxBinnedMeanDifference, std::abs(s2Summary.mean - allSummary.mean));
        if (allSummary.rms > 0)
            maxBinnedRmsRelativeDifference = std::max(maxBinnedRmsRelativeDifference,
                                                      std::abs(s2Summary.rms / allSummary.rms - 1.0));
    }

    int zetaAllBoundaryCount = 0;
    int zetaS2BoundaryCount = 0;
    std::vector<double> zetaDifferenceAwayFromBounds;
    std::vector<double> zetaAllCommonDisplayRange;
    std::vector<double> zetaS2CommonDisplayRange;
    std::vector<double> zetaDifferenceCommonDisplayRange;
    for (size_t i = 0; i < zetaAll.size(); ++i)
    {
        const bool allAtBoundary = zetaAll[i] <= -7.999 || zetaAll[i] >= 11.999;
        const bool s2AtBoundary = zetaS2[i] <= -7.999 || zetaS2[i] >= 11.999;
        zetaAllBoundaryCount += allAtBoundary;
        zetaS2BoundaryCount += s2AtBoundary;
        if (!allAtBoundary && !s2AtBoundary)
            zetaDifferenceAwayFromBounds.push_back(zetaDifference[i]);
        if (zetaAll[i] >= -2 && zetaAll[i] <= 6 && zetaS2[i] >= -2 && zetaS2[i] <= 6)
        {
            zetaAllCommonDisplayRange.push_back(zetaAll[i]);
            zetaS2CommonDisplayRange.push_back(zetaS2[i]);
            zetaDifferenceCommonDisplayRange.push_back(zetaDifference[i]);
        }
    }

    const std::string summaryPath = std::string(outputDirectory) + "/summary.txt";
    std::ofstream output(summaryPath.c_str());
    output << std::setprecision(10)
           << "beta.entries=" << betaValues.size() << '\n'
           << "zeta.entries_mc_beta_lt_0.9=" << zetaAll.size() << '\n';
    writeSummary(output, "residual_all", residualAllSummary);
    writeSummary(output, "residual_s2", residualS2Summary);
    writeSummary(output, "beta_s2_minus_all", betaDifferenceSummary);
    writeSummary(output, "zeta_all", zetaAllSummary);
    writeSummary(output, "zeta_s2", zetaS2Summary);
    writeSummary(output, "zeta_s2_minus_all", zetaDifferenceSummary);
    writeSummary(output, "zeta_s2_minus_all_away_from_bounds", summarize(zetaDifferenceAwayFromBounds));
    output << "zeta_common_display_range.entries=" << zetaAllCommonDisplayRange.size() << '\n';
    writeSummary(output, "zeta_all_common_display_range", summarize(zetaAllCommonDisplayRange));
    writeSummary(output, "zeta_s2_common_display_range", summarize(zetaS2CommonDisplayRange));
    writeSummary(output, "zeta_s2_minus_all_common_display_range", summarize(zetaDifferenceCommonDisplayRange));
    output << "zeta.correlation=" << correlation(zetaAll, zetaS2) << '\n'
           << "zeta_all.boundary_fraction=" << static_cast<double>(zetaAllBoundaryCount) / zetaAll.size() << '\n'
           << "zeta_s2.boundary_fraction=" << static_cast<double>(zetaS2BoundaryCount) / zetaS2.size() << '\n'
           << "beta.max_binned_mean_difference=" << maxBinnedMeanDifference << '\n'
           << "beta.max_binned_rms_relative_difference=" << maxBinnedRmsRelativeDifference << '\n';
    for (int bin = 0; bin < nBetaBins; ++bin)
    {
        output << "beta.bin" << bin << ".range=" << binLows[bin] << "," << binHighs[bin] << '\n'
               << "beta.bin" << bin << ".entries=" << binCounts[bin] << '\n'
               << "beta.bin" << bin << ".all_mean=" << binAllMeans[bin] << '\n'
               << "beta.bin" << bin << ".s2_mean=" << binS2Means[bin] << '\n'
               << "beta.bin" << bin << ".all_rms=" << binAllRms[bin] << '\n'
               << "beta.bin" << bin << ".s2_rms=" << binS2Rms[bin] << '\n';
    }
    output.close();

    const int colorAll = TColor::GetColor("#45364B");
    const int colorS2 = TColor::GetColor("#E55812");
    const double residualLow = std::min(quantile(residualAll, 0.005), quantile(residualS2, 0.005));
    const double residualHigh = std::max(quantile(residualAll, 0.995), quantile(residualS2, 0.995));
    const double betaDifferenceLow = quantile(betaDifference, 0.005);
    const double betaDifferenceHigh = quantile(betaDifference, 0.995);
    const double zetaLow = std::min(quantile(zetaAll, 0.005), quantile(zetaS2, 0.005));
    const double zetaHigh = std::max(quantile(zetaAll, 0.995), quantile(zetaS2, 0.995));

    TH1D hResidualAll("hResidualAll", ";1/#beta_{reco} - 1/#beta_{MC};normalized entries", 100, residualLow, residualHigh);
    TH1D hResidualS2("hResidualS2", ";1/#beta_{reco} - 1/#beta_{MC};normalized entries", 100, residualLow, residualHigh);
    TH1D hBetaDifference("hBetaDifference", ";1/#beta_{alternative} - 1/#beta_{all};entries", 100, betaDifferenceLow, betaDifferenceHigh);
    TH1D hZetaAll("hZetaAll", ";fitted #zeta (#beta_{MC} < 0.9);normalized entries", 100, zetaLow, zetaHigh);
    TH1D hZetaS2("hZetaS2", ";fitted #zeta (#beta_{MC} < 0.9);normalized entries", 100, zetaLow, zetaHigh);
    TH2D hZetaCorrelation("hZetaCorrelation", ";#zeta_{all};#zeta_{alternative}", 80, zetaLow, zetaHigh, 80, zetaLow, zetaHigh);
    for (size_t i = 0; i < residualAll.size(); ++i)
    {
        hResidualAll.Fill(residualAll[i]);
        hResidualS2.Fill(residualS2[i]);
        hBetaDifference.Fill(betaDifference[i]);
    }
    for (size_t i = 0; i < zetaAll.size(); ++i)
    {
        hZetaAll.Fill(zetaAll[i]);
        hZetaS2.Fill(zetaS2[i]);
        hZetaCorrelation.Fill(zetaAll[i], zetaS2[i]);
    }
    hResidualAll.Scale(1.0 / hResidualAll.Integral());
    hResidualS2.Scale(1.0 / hResidualS2.Integral());
    hZetaAll.Scale(1.0 / hZetaAll.Integral());
    hZetaS2.Scale(1.0 / hZetaS2.Integral());
    hResidualAll.SetLineColor(colorAll);
    hResidualS2.SetLineColor(colorS2);
    hZetaAll.SetLineColor(colorAll);
    hZetaS2.SetLineColor(colorS2);
    hResidualAll.SetLineWidth(3);
    hResidualS2.SetLineWidth(3);
    hZetaAll.SetLineWidth(3);
    hZetaS2.SetLineWidth(3);

    TCanvas canvas("canvas", "zeta scope validation", 1280, 900);
    canvas.Divide(2, 2);
    canvas.cd(1);
    hResidualAll.Draw("hist");
    hResidualS2.Draw("hist same");
    TLegend residualLegend(0.62, 0.72, 0.88, 0.88);
    residualLegend.SetBorderSize(0);
    residualLegend.AddEntry(&hResidualAll, "#zeta on S1-S3", "l");
    residualLegend.AddEntry(&hResidualS2, alternativeLabel, "l");
    residualLegend.Draw();
    canvas.cd(2);
    gPad->SetLogy();
    hBetaDifference.Draw("hist");
    canvas.cd(3);
    hZetaAll.Draw("hist");
    hZetaS2.Draw("hist same");
    TLegend zetaLegend(0.62, 0.72, 0.88, 0.88);
    zetaLegend.SetBorderSize(0);
    zetaLegend.AddEntry(&hZetaAll, "#zeta on S1-S3", "l");
    zetaLegend.AddEntry(&hZetaS2, alternativeLabel, "l");
    zetaLegend.Draw();
    canvas.cd(4);
    gPad->SetRightMargin(0.16);
    hZetaCorrelation.Draw("colz");
    const std::string pdfPath = std::string(outputDirectory) + "/comparison.pdf";
    const std::string pngPath = std::string(outputDirectory) + "/comparison.png";
    canvas.Print(pdfPath.c_str());
    canvas.Print(pngPath.c_str());

    std::cout << "Wrote " << summaryPath << std::endl;
    std::cout << "Wrote " << pdfPath << std::endl;
    std::cout << "Wrote " << pngPath << std::endl;
}
