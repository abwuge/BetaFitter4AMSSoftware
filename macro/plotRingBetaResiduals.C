#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TF1.h>
#include <TProfile.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TROOT.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <TFitResult.h>

/**
 * Draw plots with log(Ek/n) as x-axis (displayed as Ek/n) and means of
 * 1/linearBeta - 1/mcBeta[17] and 1/nonlinearBeta - 1/mcBeta[17] as y-axis,
 * separated by concentric ring regions
 *
 * @param fileName Input ROOT file path
 * @param nucleonNumber Number of nucleons per nucleus, default is 4
 * @param outputName Output PDF file name
 */
void plotRingBetaResiduals(
    std::string fileName = "test.root",
    int nucleonNumber = 4,
    const char *outputName = "test_ring_beta_residuals.pdf")
{
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetEndErrorSize(20);

    // Open file
    // ------------------------------------------------------------------------
    if (fileName.size() < 5 || fileName.substr(fileName.size() - 5) != ".root")
        fileName += ".root";

    TFile *file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        fileName = "results/" + fileName;
        std::cout << "Trying to open file: " << fileName << std::endl;
        file = TFile::Open(fileName.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "Error: Unable to find/open file " << fileName << std::endl;
            return;
        }
    }

    // Get tree
    // ------------------------------------------------------------------------
    TTree *tree = (TTree *)file->Get("betaDiff");
    if (!tree)
    {
        std::cerr << "Error: Could not find betaDiff tree in " << fileName << std::endl;
        file->Close();
        delete file;
        return;
    }
    if (tree->GetEntries() == 0)
    {
        std::cerr << "Error: Empty tree, no data to plot" << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Set branch addresses
    // ------------------------------------------------------------------------
    float mevmom1[21]{};
    float mevcoo1[21][3]{};
    float mcBeta[21]{};
    float linearBeta{};
    float nonlinearBeta{};

    tree->SetBranchAddress("mevmom1", mevmom1);
    tree->SetBranchAddress("mevcoo1", mevcoo1);
    tree->SetBranchAddress("mcBeta", mcBeta);
    tree->SetBranchAddress("linearBeta", &linearBeta);
    tree->SetBranchAddress("nonlinearBeta", &nonlinearBeta);

    // Set kinetic energy range
    // ------------------------------------------------------------------------
    double minLogEk = -0.5;
    double maxLogEk = 3;

    // Convert log range to actual Ek values for display
    double minEk = pow(10, minLogEk);
    double maxEk = pow(10, maxLogEk);

    // Define histograms
    // ------------------------------------------------------------------------
    const int nBinsX = 100; // Number of bins for log(Ek/n)
    const int nBinsY = 100; // Number of bins for residuals

    // Fixed residual range for better visualization
    const double yMinRes = -0.2;
    const double yMaxRes = 0.2;

    // Create histograms with logarithmic x-axis bins
    // Use log binning internally but display actual Ek values
    double logEkBins[nBinsX + 1];
    double logEkBinWidth = (maxLogEk - minLogEk) / nBinsX;

    for (int i = 0; i <= nBinsX; ++i)
        logEkBins[i] = minLogEk + i * logEkBinWidth;

    // Define ring regions (concentric rings, each 8cm wide)
    // Extended range -60 to 60 cm
    const int numRings = 8;      // 8 rings (0-8, 8-16, 16-24, 24-32, 32-40, 40-48, 48-56, 56-60)
    const float ringWidth = 6.0; // cm

    // Create histograms for each ring region (linear and nonlinear)
    TH2F *hLinearResVsLogEk_Rings[numRings];
    TH2F *hNonlinearResVsLogEk_Rings[numRings];

    for (int i = 0; i < numRings; ++i)
    {
        float innerRadius = i * ringWidth;
        float outerRadius = (i + 1) * ringWidth;

        // Linear beta histograms
        hLinearResVsLogEk_Rings[i] = new TH2F(
            Form("hLinearResVsLogEk_Ring%d", i),
            Form("Ring %d (%d - %d cm): Linear Beta Residuals vs. E_{k}/n;E_{k}/n [GeV/nucleon];1/#beta_{linear} - 1/#beta_{MC}",
                 i, (int)innerRadius, (int)outerRadius),
            nBinsX, logEkBins, nBinsY, yMinRes, yMaxRes);

        // Nonlinear beta histograms
        hNonlinearResVsLogEk_Rings[i] = new TH2F(
            Form("hNonlinearResVsLogEk_Ring%d", i),
            Form("Ring %d (%d - %d cm): Nonlinear Beta Residuals vs. E_{k}/n;E_{k}/n [GeV/nucleon];1/#beta_{nonlinear} - 1/#beta_{MC}",
                 i, (int)innerRadius, (int)outerRadius),
            nBinsX, logEkBins, nBinsY, yMinRes, yMaxRes);

        // Set minimum value to 1 for log scale visibility
        hLinearResVsLogEk_Rings[i]->SetMinimum(1);
        hNonlinearResVsLogEk_Rings[i]->SetMinimum(1);
    }

    // Fill histograms
    // ------------------------------------------------------------------------
    for (int i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);

        // Skip invalid data
        if (mevmom1[2] == -1000 || mcBeta[2] == -1000)
            continue;

        // Skip invalid coordinates
        if (mevcoo1[10][0] == -1000 || mevcoo1[10][1] == -1000 || mevcoo1[10][2] == -1000)
            continue;

        // Calculate kinetic energy Ek = mc^2 * (gamma - 1)
        double momentum = mevmom1[2];                                         // GeV/c
        double mass = momentum / mcBeta[2] * sqrt(1 - mcBeta[2] * mcBeta[2]); // GeV/c^2
        double gamma = 1.0 / sqrt(1 - mcBeta[2] * mcBeta[2]);
        double kineticEnergy = mass * (gamma - 1.0); // GeV

        // Calculate kinetic energy per nucleon and take log
        double ekPerNucleon = kineticEnergy / nucleonNumber; // GeV/nucleon
        double logEkPerNucleon = TMath::Log10(ekPerNucleon);

        // Calculate residuals
        double linearRes = 1.0 / linearBeta - 1.0 / mcBeta[17];
        double nonlinearRes = 1.0 / nonlinearBeta - 1.0 / mcBeta[17];

        // Calculate radius from center
        double x = mevcoo1[10][0];
        double y = mevcoo1[10][1];
        double r = TMath::Sqrt(x * x + y * y);

        // Determine which ring this point belongs to
        int ringIndex = TMath::Floor(r / ringWidth);

        // Fill appropriate histogram if within range
        if (ringIndex >= 0 && ringIndex < numRings)
        {
            hLinearResVsLogEk_Rings[ringIndex]->Fill(logEkPerNucleon, linearRes);
            hNonlinearResVsLogEk_Rings[ringIndex]->Fill(logEkPerNucleon, nonlinearRes);
        }
    }

    // Create canvas
    // ------------------------------------------------------------------------
    TCanvas *canvas = new TCanvas("canvas", "", 3508, 2480); // A4 size
    canvas->SetLeftMargin(0.16);
    canvas->SetRightMargin(0.12);
    canvas->SetGridx();
    canvas->SetGridy();
    canvas->SetLogz();

    // Create PDF file
    canvas->Print(Form("%s[", outputName));

    // Create info text
    TPaveText *infoText = new TPaveText(0.2, 0.92, 0.8, 0.98, "NDC");
    infoText->SetFillColor(0);
    infoText->SetBorderSize(0);
    infoText->AddText(Form("Nucleon number n = %d", nucleonNumber));

    // Set X-axis to display actual Ek values instead of log values
    for (int i = 0; i < numRings; ++i)
    {
        TAxis *xAxisLinear = hLinearResVsLogEk_Rings[i]->GetXaxis();
        TAxis *xAxisNonlinear = hNonlinearResVsLogEk_Rings[i]->GetXaxis();

        xAxisLinear->SetMoreLogLabels();
        xAxisLinear->SetNoExponent();
        xAxisNonlinear->SetMoreLogLabels();
        xAxisNonlinear->SetNoExponent();

        // Convert bin labels from log to actual values
        for (int j = 1; j <= nBinsX; j++)
        {
            double logValue = xAxisLinear->GetBinCenter(j);
            double actualValue = pow(10, logValue);

            // Only set labels for some bins to avoid overcrowding
            if (j % 5 == 0 || j == 1 || j == nBinsX)
            {
                if (actualValue < 0.1)
                {
                    xAxisLinear->SetBinLabel(j, Form("%.3f", actualValue));
                    xAxisNonlinear->SetBinLabel(j, Form("%.3f", actualValue));
                }
                else if (actualValue < 1)
                {
                    xAxisLinear->SetBinLabel(j, Form("%.2f", actualValue));
                    xAxisNonlinear->SetBinLabel(j, Form("%.2f", actualValue));
                }
                else
                {
                    xAxisLinear->SetBinLabel(j, Form("%.1f", actualValue));
                    xAxisNonlinear->SetBinLabel(j, Form("%.1f", actualValue));
                }
            }
        }
    }

    // Fit
    // ------------------------------------------------------------------------
    // Vectors to store fit results for each ring
    std::vector<double> logEkValues;
    std::vector<double> ekValues; // Actual Ek values for display
    std::vector<std::vector<double>> linearMeans(numRings);
    std::vector<std::vector<double>> linearErrors(numRings);
    std::vector<std::vector<double>> nonlinearMeans(numRings);
    std::vector<std::vector<double>> nonlinearErrors(numRings);

    logEkValues.reserve(nBinsX);
    ekValues.reserve(nBinsX);

    for (int i = 0; i < numRings; ++i)
    {
        linearMeans[i].reserve(nBinsX);
        linearErrors[i].reserve(nBinsX);
        nonlinearMeans[i].reserve(nBinsX);
        nonlinearErrors[i].reserve(nBinsX);
    }

    // Create fit function
    TF1 *fGaus = new TF1("fGaus", "gaus", yMinRes, yMaxRes);
    fGaus->SetParLimits(1, yMinRes, yMaxRes);
    fGaus->SetNpx(1000);

    // Fit each column for each histogram (both linear and nonlinear)
    for (int ringIdx = 0; ringIdx < numRings; ++ringIdx)
    {
        TH2F *linearHist = hLinearResVsLogEk_Rings[ringIdx];
        TH2F *nonlinearHist = hNonlinearResVsLogEk_Rings[ringIdx];

        for (int i = 0; i < nBinsX; ++i)
        {
            double logEkLow = logEkBins[i];
            double logEkHigh = logEkBins[i + 1];
            double binCenter = (logEkLow + logEkHigh) / 2.0;
            double ekValue = pow(10, binCenter); // Convert to actual Ek value

            // Process linear histograms
            TH1D *linearProj = linearHist->ProjectionY(Form("linear_proj_%d_%d", ringIdx, i), i + 1, i + 1);
            Int_t linearEntries = linearProj->GetEntries();

            if (linearEntries > 50)
            {
                fGaus->SetParameters(linearProj->GetMaximum(), linearProj->GetMean(), linearProj->GetRMS());
                TFitResultPtr fitResult = linearProj->Fit(fGaus, "SQNR");

                if (fitResult.Get() != nullptr && !fitResult->IsEmpty() && fitResult->Status() == 0)
                {
                    // Store fit results
                    if (ringIdx == 0)
                    { // Only store Ek values once
                        if (std::find(logEkValues.begin(), logEkValues.end(), binCenter) == logEkValues.end())
                        {
                            logEkValues.push_back(binCenter);
                            ekValues.push_back(ekValue);
                        }
                    }

                    linearMeans[ringIdx].push_back(fitResult->Parameter(1));
                    linearErrors[ringIdx].push_back(fitResult->Parameter(2));
                }
            }

            // Process nonlinear histograms
            TH1D *nonlinearProj = nonlinearHist->ProjectionY(Form("nonlinear_proj_%d_%d", ringIdx, i), i + 1, i + 1);
            Int_t nonlinearEntries = nonlinearProj->GetEntries();

            if (nonlinearEntries > 50)
            {
                fGaus->SetParameters(nonlinearProj->GetMaximum(), nonlinearProj->GetMean(), nonlinearProj->GetRMS());
                TFitResultPtr fitResult = nonlinearProj->Fit(fGaus, "SQNR");

                if (fitResult.Get() != nullptr && !fitResult->IsEmpty() && fitResult->Status() == 0)
                {
                    nonlinearMeans[ringIdx].push_back(fitResult->Parameter(1));
                    nonlinearErrors[ringIdx].push_back(fitResult->Parameter(2));
                }
            }

            delete linearProj;
            delete nonlinearProj;
        }
    }

    // Create mean and error graphs for each ring
    // ------------------------------------------------------------------------
    // Colors for different rings
    const int ringColors[numRings] = {kBlue, kRed, kGreen + 2, kViolet, kOrange + 7, kCyan + 1, kMagenta + 2, kBlack};

    // Graphs for linear means and errors
    TGraph *grLinearMean_Rings[numRings];
    TGraph *grLinearError_Rings[numRings];

    // Graphs for nonlinear means and errors
    TGraph *grNonlinearMean_Rings[numRings];
    TGraph *grNonlinearError_Rings[numRings];

    for (int i = 0; i < numRings; ++i)
    {
        // Process linear data if available
        if (!linearMeans[i].empty())
        {
            // Mean graphs
            grLinearMean_Rings[i] = new TGraph(ekValues.size(), ekValues.data(), linearMeans[i].data());
            grLinearMean_Rings[i]->SetMarkerStyle(20 + i);
            grLinearMean_Rings[i]->SetMarkerSize(1.5);
            grLinearMean_Rings[i]->SetMarkerColor(ringColors[i]);
            grLinearMean_Rings[i]->SetLineColor(ringColors[i]);
            grLinearMean_Rings[i]->SetLineWidth(2);

            // Error graphs
            grLinearError_Rings[i] = new TGraph(ekValues.size(), ekValues.data(), linearErrors[i].data());
            grLinearError_Rings[i]->SetMarkerStyle(24 + i);
            grLinearError_Rings[i]->SetMarkerSize(1.5);
            grLinearError_Rings[i]->SetMarkerColor(ringColors[i]);
            grLinearError_Rings[i]->SetLineColor(ringColors[i]);
            grLinearError_Rings[i]->SetLineStyle(2);
            grLinearError_Rings[i]->SetLineWidth(2);
        }

        // Process nonlinear data if available
        if (!nonlinearMeans[i].empty())
        {
            // Mean graphs
            grNonlinearMean_Rings[i] = new TGraph(ekValues.size(), ekValues.data(), nonlinearMeans[i].data());
            grNonlinearMean_Rings[i]->SetMarkerStyle(20 + i);
            grNonlinearMean_Rings[i]->SetMarkerSize(1.5);
            grNonlinearMean_Rings[i]->SetMarkerColor(ringColors[i]);
            grNonlinearMean_Rings[i]->SetLineColor(ringColors[i]);
            grNonlinearMean_Rings[i]->SetLineWidth(2);

            // Error graphs
            grNonlinearError_Rings[i] = new TGraph(ekValues.size(), ekValues.data(), nonlinearErrors[i].data());
            grNonlinearError_Rings[i]->SetMarkerStyle(24 + i);
            grNonlinearError_Rings[i]->SetMarkerSize(1.5);
            grNonlinearError_Rings[i]->SetMarkerColor(ringColors[i]);
            grNonlinearError_Rings[i]->SetLineColor(ringColors[i]);
            grNonlinearError_Rings[i]->SetLineStyle(2);
            grNonlinearError_Rings[i]->SetLineWidth(2);
        }
    }

    // Create zero residual reference line
    TF1 *perfectResidual = new TF1("perfectResidual", "0", minEk, maxEk);
    perfectResidual->SetLineColor(kBlack);
    perfectResidual->SetLineStyle(2);

    // Find optimal y-axis range for means and errors based on actual data
    double linearMeanMin = 0, linearMeanMax = 0;
    double linearErrorMin = 0, linearErrorMax = 0;
    double nonlinearMeanMin = 0, nonlinearMeanMax = 0;
    double nonlinearErrorMin = 0, nonlinearErrorMax = 0;

    // Helper function to find min/max values
    auto updateRange = [](const std::vector<double> &data, double &minVal, double &maxVal)
    {
        if (data.empty())
            return;

        for (double val : data)
        {
            if (std::isnan(val) || std::isinf(val))
                continue;

            if (val < minVal)
                minVal = val;
            if (val > maxVal)
                maxVal = val;
        }
    };

    // Find min/max for linear means across all rings
    linearMeanMin = std::numeric_limits<double>::max();
    linearMeanMax = std::numeric_limits<double>::lowest();

    for (int i = 0; i < numRings; ++i)
    {
        updateRange(linearMeans[i], linearMeanMin, linearMeanMax);
    }

    // Add some padding (20%)
    double linearMeanPadding = (linearMeanMax - linearMeanMin) * 0.2;
    linearMeanMin = std::max(linearMeanMin - linearMeanPadding, -0.1);
    linearMeanMax = std::min(linearMeanMax + linearMeanPadding, 0.1);

    // Ensure reasonable limits if data is very sparse
    if (linearMeanMin > -0.005)
        linearMeanMin = -0.005;
    if (linearMeanMax < 0.005)
        linearMeanMax = 0.005;

    // Find min/max for linear errors across all rings
    linearErrorMin = 0; // Errors are always positive
    linearErrorMax = std::numeric_limits<double>::lowest();

    for (int i = 0; i < numRings; ++i)
    {
        updateRange(linearErrors[i], linearErrorMin, linearErrorMax);
    }

    // Add some padding (20%)
    double linearErrorPadding = linearErrorMax * 0.2;
    linearErrorMax = std::min(linearErrorMax + linearErrorPadding, 0.1);

    // Ensure reasonable limits if data is very sparse
    if (linearErrorMax < 0.01)
        linearErrorMax = 0.01;

    // Find min/max for nonlinear means across all rings
    nonlinearMeanMin = std::numeric_limits<double>::max();
    nonlinearMeanMax = std::numeric_limits<double>::lowest();

    for (int i = 0; i < numRings; ++i)
    {
        updateRange(nonlinearMeans[i], nonlinearMeanMin, nonlinearMeanMax);
    }

    // Add some padding (20%)
    double nonlinearMeanPadding = (nonlinearMeanMax - nonlinearMeanMin) * 0.2;
    nonlinearMeanMin = std::max(nonlinearMeanMin - nonlinearMeanPadding, -0.1);
    nonlinearMeanMax = std::min(nonlinearMeanMax + nonlinearMeanPadding, 0.1);

    // Ensure reasonable limits if data is very sparse
    if (nonlinearMeanMin > -0.005)
        nonlinearMeanMin = -0.005;
    if (nonlinearMeanMax < 0.005)
        nonlinearMeanMax = 0.005;

    // Find min/max for nonlinear errors across all rings
    nonlinearErrorMin = 0; // Errors are always positive
    nonlinearErrorMax = std::numeric_limits<double>::lowest();

    for (int i = 0; i < numRings; ++i)
    {
        updateRange(nonlinearErrors[i], nonlinearErrorMin, nonlinearErrorMax);
    }

    // Add some padding (20%)
    double nonlinearErrorPadding = nonlinearErrorMax * 0.2;
    nonlinearErrorMax = std::min(nonlinearErrorMax + nonlinearErrorPadding, 0.1);

    // Ensure reasonable limits if data is very sparse
    if (nonlinearErrorMax < 0.01)
        nonlinearErrorMax = 0.01;

    // Set up dummy histograms with dynamic ranges for linear data
    TH2F *dummyHistLinearMean = new TH2F("dummyHistLinearMean",
                                         "Linear Beta Residual Means vs. E_{k}/n;E_{k}/n [GeV/nucleon];1/#beta_{linear} - 1/#beta_{MC}",
                                         100, minEk, maxEk, 100, linearMeanMin, linearMeanMax);
    dummyHistLinearMean->SetStats(0);

    TH2F *dummyHistLinearError = new TH2F("dummyHistLinearError",
                                          "Linear Beta Residual Errors vs. E_{k}/n;E_{k}/n [GeV/nucleon];#sigma(1/#beta_{linear} - 1/#beta_{MC})",
                                          100, minEk, maxEk, 100, linearErrorMin, linearErrorMax);
    dummyHistLinearError->SetStats(0);

    // Set up dummy histograms with dynamic ranges for nonlinear data
    TH2F *dummyHistNonlinearMean = new TH2F("dummyHistNonlinearMean",
                                            "Nonlinear Beta Residual Means vs. E_{k}/n;E_{k}/n [GeV/nucleon];1/#beta_{nonlinear} - 1/#beta_{MC}",
                                            100, minEk, maxEk, 100, nonlinearMeanMin, nonlinearMeanMax);
    dummyHistNonlinearMean->SetStats(0);

    TH2F *dummyHistNonlinearError = new TH2F("dummyHistNonlinearError",
                                             "Nonlinear Beta Residual Errors vs. E_{k}/n;E_{k}/n [GeV/nucleon];#sigma(1/#beta_{nonlinear} - 1/#beta_{MC})",
                                             100, minEk, maxEk, 100, nonlinearErrorMin, nonlinearErrorMax);
    dummyHistNonlinearError->SetStats(0);

    // Draw linear mean plots
    // ------------------------------------------------------------------------
    canvas->SetLogx(); // Set logarithmic x-axis

    // Linear means comparison for all rings
    dummyHistLinearMean->Draw();
    perfectResidual->Draw("SAME");

    TLegend *legendLinearMean = new TLegend(0.65, 0.15, 0.85, 0.35);
    legendLinearMean->SetBorderSize(0);

    for (int i = 0; i < numRings; ++i)
    {
        if (linearMeans[i].empty())
            continue;
        grLinearMean_Rings[i]->Draw("LP");
        legendLinearMean->AddEntry(grLinearMean_Rings[i],
                                   Form("Ring %d (%d-%d cm)", i, (int)(i * ringWidth), (int)((i + 1) * ringWidth)),
                                   "lp");
    }

    legendLinearMean->Draw();
    infoText->Draw();
    canvas->Print(outputName);

    // Draw linear error plots
    // ------------------------------------------------------------------------

    // Linear errors comparison for all rings
    dummyHistLinearError->Draw();

    TLegend *legendLinearError = new TLegend(0.65, 0.15, 0.85, 0.35);
    legendLinearError->SetBorderSize(0);

    for (int i = 0; i < numRings; ++i)
    {
        if (linearErrors[i].empty())
            continue;
        grLinearError_Rings[i]->Draw("LP");
        legendLinearError->AddEntry(grLinearError_Rings[i],
                                    Form("Ring %d (%d-%d cm)", i, (int)(i * ringWidth), (int)((i + 1) * ringWidth)),
                                    "lp");
    }

    legendLinearError->Draw();
    infoText->Draw();
    canvas->Print(outputName);

    // Draw nonlinear mean plots
    // ------------------------------------------------------------------------

    // Nonlinear means comparison for all rings
    dummyHistNonlinearMean->Draw();
    perfectResidual->Draw("SAME");

    TLegend *legendNonlinearMean = new TLegend(0.65, 0.15, 0.85, 0.35);
    legendNonlinearMean->SetBorderSize(0);

    for (int i = 0; i < numRings; ++i)
    {
        if (nonlinearMeans[i].empty())
            continue;
        grNonlinearMean_Rings[i]->Draw("LP");
        legendNonlinearMean->AddEntry(grNonlinearMean_Rings[i],
                                      Form("Ring %d (%d-%d cm)", i, (int)(i * ringWidth), (int)((i + 1) * ringWidth)),
                                      "lp");
    }

    legendNonlinearMean->Draw();
    infoText->Draw();
    canvas->Print(outputName);

    // Draw nonlinear error plots
    // ------------------------------------------------------------------------

    // Nonlinear errors comparison for all rings
    dummyHistNonlinearError->Draw();

    TLegend *legendNonlinearError = new TLegend(0.65, 0.15, 0.85, 0.35);
    legendNonlinearError->SetBorderSize(0);

    for (int i = 0; i < numRings; ++i)
    {
        if (nonlinearErrors[i].empty())
            continue;
        grNonlinearError_Rings[i]->Draw("LP");
        legendNonlinearError->AddEntry(grNonlinearError_Rings[i],
                                       Form("Ring %d (%d-%d cm)", i, (int)(i * ringWidth), (int)((i + 1) * ringWidth)),
                                       "lp");
    }

    legendNonlinearError->Draw();
    infoText->Draw();
    canvas->Print(outputName);

    // Close PDF file
    canvas->Print(Form("%s]", outputName));

    std::cout << "Ring-based E_k/n vs Beta residuals plot saved to: " << outputName << std::endl;

    // Cleanup
    delete dummyHistLinearMean;
    delete dummyHistLinearError;
    delete dummyHistNonlinearMean;
    delete dummyHistNonlinearError;

    for (int i = 0; i < numRings; ++i)
    {
        delete hLinearResVsLogEk_Rings[i];
        delete hNonlinearResVsLogEk_Rings[i];

        // Only delete if they were created (had data)
        if (!linearMeans[i].empty())
        {
            delete grLinearMean_Rings[i];
            delete grLinearError_Rings[i];
        }

        if (!nonlinearMeans[i].empty())
        {
            delete grNonlinearMean_Rings[i];
            delete grNonlinearError_Rings[i];
        }
    }

    file->Close();
    delete file;
}