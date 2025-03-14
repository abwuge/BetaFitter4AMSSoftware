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

/**
 * Get Z value and energyLossScale for a file from README.md
 * @param fileName Name of the ROOT file
 * @return pair of Z value and energyLossScale. Returns {-1, -1} if not found
 */
std::pair<int, double> getParamsFromReadme(const std::string &fileName)
{
    std::ifstream readme("results/README.md");
    if (!readme)
    {
        std::cerr << "Error: Could not open results/README.md" << std::endl;
        return {-1, -1};
    }

    std::string line;
    std::string targetFile = fileName;
    // If fileName starts with "results/", remove it
    if (targetFile.substr(0, 8) == "results/")
    {
        targetFile = targetFile.substr(8);
    }
    else
    {
        std::cerr << "Warning: File name does not start with 'results/'"
                  << "\nCannot find the file in README.md" << std::endl;
        return {-1, -1};
    }

    int latestZ = -1;
    double energyLossScale = -1.0;
    // Read file from bottom to top to get the latest entry
    std::vector<std::string> lines;
    while (std::getline(readme, line))
    {
        lines.push_back(line);
    }

    // Search from the end to find the latest matching entry
    for (auto it = lines.rbegin(); it != lines.rend(); ++it)
    {
        size_t filePos = it->find("FILE = " + targetFile);
        if (filePos != std::string::npos)
        {
            // Get Z value
            size_t zPos = it->find("Z = ", filePos);
            if (zPos != std::string::npos)
            {
                zPos += 4; // Skip "Z = "
                size_t commaPos = it->find(",", zPos);
                if (commaPos != std::string::npos)
                {
                    std::string zStr = it->substr(zPos, commaPos - zPos);
                    try
                    {
                        latestZ = std::stoi(zStr);
                    }
                    catch (...)
                    {
                        continue;
                    }
                }
            }

            // Get energyLossScale value
            size_t scalePos = it->find("energyLossScale = ", filePos);
            if (scalePos != std::string::npos)
            {
                scalePos += 17; // Skip "energyLossScale = "
                size_t commaPos = it->find(",", scalePos);
                if (commaPos != std::string::npos)
                {
                    std::string scaleStr = it->substr(scalePos, commaPos - scalePos);
                    try
                    {
                        energyLossScale = std::stod(scaleStr);
                    }
                    catch (...)
                    {
                        continue;
                    }
                }
            }

            if (latestZ > 0 && energyLossScale > 0)
            {
                break;
            }
        }
    }

    return {latestZ, energyLossScale};
}

/**
 * Draw beta comparison plots from a ROOT file containing betaTree
 *
 * @param fileName Path to the input ROOT file
 * @param outputName Output file name (default: auto-generated based on Z and ELS values)
 */
void plotBetaComparison(std::string fileName = "test.root",
                        const char *outputName = nullptr)
{
    // Set batch mode to avoid GUI related issues
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetEndErrorSize(20);

    // Open the ROOT file
    if (fileName.size() < 5 || fileName.substr(fileName.size() - 5) != ".root")
        fileName += ".root";

    TFile *file = TFile::Open(fileName.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        fileName = "results/" + fileName;
        std::cout << "Try to open file: " << fileName << std::endl;

        file = TFile::Open(fileName.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "Error: Unable to find/open file " << fileName << std::endl;
            return;
        }
    }

    // Auto-generate output name if not specified

    // Get Z value and energyLossScale from the file
    auto params = getParamsFromReadme(fileName);
    int zValue = params.first;
    double energyLossScale = params.second;
    std::string actualOutputName;

    if (!outputName)
    {
        if (zValue > 0 && energyLossScale > 0)
            actualOutputName = Form("test_Z%d_ELS%.1f.pdf", zValue, energyLossScale);
        else
            actualOutputName = "test_beta_comparison.pdf";
        outputName = actualOutputName.c_str();
    }

    // Get the beta reconstruction tree
    TTree *tree = (TTree *)file->Get("betaTree");
    if (!tree)
    {
        std::cerr << "Error: Could not find betaTree in " << fileName << std::endl;
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

    // Use a more efficient approach to get beta range - default values
    double mcBetaMin = tree->GetMinimum("mcBeta");
    double mcBetaMax = tree->GetMaximum("mcBeta");

    // Define variables for reading from tree
    double mcBeta, linearBeta, nonlinearBeta;

    // Set branch addresses
    tree->SetBranchAddress("mcBeta", &mcBeta);
    tree->SetBranchAddress("linearBeta", &linearBeta);
    tree->SetBranchAddress("nonlinearBeta", &nonlinearBeta);

    // Create histograms for the 2D plots
    // Using reasonable beta range from 0.1 to 1.0
    int nBinsX = 40;                           // Number of bins in x direction
    int nBinsY = 100;                          // Number of bins in y direction for the residual plot
    double xMin = mcBetaMin, xMax = mcBetaMax; // Range for the beta values

    // Helper function to get range using quantiles
    auto getQuantileRange = [&tree](const char *branchName, double lowQuantile = 0.01, double highQuantile = 0.99) -> std::pair<double, double>
    {
        // Create a high-precision temporary histogram with reasonable initial range
        TH1D hTemp(Form("hTemp_%s", branchName), "", 10000, -1, 1);

        // Use TFormula for efficient calculation
        TString formula = TString::Format("1/%s - 1/mcBeta", branchName);
        tree->Draw(Form("%s>>hTemp_%s", formula.Data(), branchName), "", "goff");

        // Calculate quantiles
        Double_t xq[2] = {lowQuantile, highQuantile};
        Double_t yq[2] = {0, 0};
        hTemp.GetQuantiles(2, yq, xq);

        std::cout << "Range for " << branchName << " residuals using quantiles: [" << yq[0] << ", " << yq[1] << "]" << std::endl;
        return std::make_pair(yq[0], yq[1]);
    };

    // Get ranges for residuals using quantiles
    auto nonlinearResRange = getQuantileRange("nonlinearBeta", 1e-4, 1 - 1e-4);
    auto linearResRange = getQuantileRange("linearBeta", 1e-4, 1 - 1e-6);

    // Use independent ranges for first two pages
    double yMinResNL = nonlinearResRange.first;
    double yMaxResNL = nonlinearResRange.second;

    double yMinResL = linearResRange.first;
    double yMaxResL = linearResRange.second;

    // Use the wider range for the comparison page
    double yMinRes = std::min(nonlinearResRange.first, linearResRange.first);
    double yMaxRes = std::max(nonlinearResRange.second, linearResRange.second);
    double margin = 0.05 * (yMaxRes - yMinRes);
    yMinRes -= margin;
    yMaxRes += margin;

    // Create histograms for the residuals plot with independent ranges
    TH2F *hNonlinearResVsMC = new TH2F("hNonlinearResVsMC",
                                       "Non-linear Reconstruction Residuals;#beta_{MC};1/#beta_{non-linear} - 1/#beta_{MC}",
                                       nBinsX, xMin, xMax, nBinsY, yMinResNL, yMaxResNL);
    TH2F *hLinearResVsMC = new TH2F("hLinearResVsMC",
                                    "Linear Reconstruction Residuals;#beta_{MC};1/#beta_{linear} - 1/#beta_{MC}",
                                    nBinsX, xMin, xMax, nBinsY, yMinResL, yMaxResL);

    hNonlinearResVsMC->SetMinimum(1);
    hLinearResVsMC->SetMinimum(1);

    // Fill 2D histograms efficiently using TFormula
    tree->Draw("1/nonlinearBeta - 1/mcBeta:mcBeta>>hNonlinearResVsMC", "", "goff");
    tree->Draw("1/linearBeta - 1/mcBeta:mcBeta>>hLinearResVsMC", "", "goff");

    // Arrays to store fit results and weights
    const int nProfiles = nBinsX;
    double mcBetaValues[nProfiles];
    double nonlinearMean[nProfiles], nonlinearError[nProfiles], nonlinearWeight[nProfiles];
    double linearMean[nProfiles], linearError[nProfiles], linearWeight[nProfiles];

    // Print table header before the loop
    std::cout << std::endl;
    std::cout << "Bin\tBeta\t\tProj_NL\tProj_L\tNL_mean\t\tL_mean" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    double binWidth = (xMax - xMin) / nBinsX;
    for (int bin = 0; bin < nBinsX; ++bin)
    {
        double binCenter = xMin + (bin + 0.5) * binWidth;
        int binIdx = bin + 1; // ROOT histograms are 1-indexed

        // Calculate bin range for projection with wider window
        int binWindow = 1; // Use ±1 bins for better statistics
        int binLow = TMath::Max(1, binIdx - binWindow);
        int binHigh = TMath::Min(nBinsX, binIdx + binWindow);

        mcBetaValues[bin] = binCenter;

        // Process nonlinear beta residual distribution with explicit options
        TH1D *projNonlinear = hNonlinearResVsMC->ProjectionY(Form("projNL_%d", bin), binLow, binHigh, "e"); // 'e' option to keep errors
        Int_t projNLEntries = projNonlinear->GetEntries();

        if (projNLEntries) // Only fit if we have enough statistics
        {
            // Add reasonable initial parameters for the fit
            projNonlinear->Fit("gaus", "QN", "", yMinRes, yMaxRes); // Add fit range
            TF1 *fitNL = projNonlinear->GetFunction("gaus");
            if (fitNL && fitNL->GetProb() > 0.01) // Check if fit is reasonable
            {
                nonlinearMean[bin] = fitNL->GetParameter(1);
                nonlinearError[bin] = fitNL->GetParameter(2);
            }
            else
            {
                nonlinearMean[bin] = projNonlinear->GetMean(); // Use histogram mean if fit fails
                nonlinearError[bin] = projNonlinear->GetRMS();
            }
        }
        else
        {
            nonlinearMean[bin] = 0;
            nonlinearError[bin] = 0;
        }
        delete projNonlinear;

        // Process linear beta residual distribution with same improvements
        TH1D *projLinear = hLinearResVsMC->ProjectionY(Form("projL_%d", bin), binLow, binHigh, "e");
        Int_t projLEntries = projLinear->GetEntries();

        if (projLEntries)
        {
            projLinear->Fit("gaus", "QN", "", yMinRes, yMaxRes);
            TF1 *fitL = projLinear->GetFunction("gaus");
            if (fitL && fitL->GetProb() > 0.01)
            {
                linearMean[bin] = fitL->GetParameter(1);
                linearError[bin] = fitL->GetParameter(2);
            }
            else
            {
                linearMean[bin] = projLinear->GetMean();
                linearError[bin] = projLinear->GetRMS();
            }
        }
        else
        {
            linearMean[bin] = 0;
            linearError[bin] = 0;
        }
        delete projLinear;

        // Print data in aligned columns using tabs and store weights
        nonlinearWeight[bin] = projNLEntries;
        linearWeight[bin] = projLEntries;
        printf("%d\t%.6f\t%d\t%d\t%.6f\t%.6f\n",
               bin, binCenter, projNLEntries, projLEntries,
               nonlinearMean[bin], linearMean[bin]);
    }

    // Create canvas with multiple pages
    TCanvas *canvas = new TCanvas("canvas", "Beta Residuals Comparison", 3508, 2480);
    canvas->Print(Form("%s[", outputName)); // Open PDF file

    // First page - Non-linear beta residuals
    TCanvas *c1 = new TCanvas("c1", "Non-linear Beta Residuals", 3508, 2480);
    c1->SetLeftMargin(0.16);
    c1->SetRightMargin(0.11);
    c1->SetGridx();
    c1->SetGridy();
    c1->SetLogz();

    hNonlinearResVsMC->Draw("COLZ");

    // Perfect residual reference line (zero residual)
    TF1 *perfectLine1 = new TF1("perfectLine1", "0", xMin, xMax);
    perfectLine1->SetLineColor(kRed);
    perfectLine1->SetLineStyle(2);
    perfectLine1->Draw("SAME");

    // Create and draw TGraphErrors for nonlinear residuals
    TGraphErrors *grNonlinear = new TGraphErrors(nProfiles, mcBetaValues, nonlinearMean, 0, nonlinearError);
    grNonlinear->SetMarkerStyle(20);
    grNonlinear->SetMarkerColor(kBlack);
    grNonlinear->SetMarkerSize(3.0);
    grNonlinear->Draw("P");

    c1->Print(outputName);

    // Second page - Linear beta residuals
    TCanvas *c2 = new TCanvas("c2", "Linear Beta Residuals", 3508, 2480);
    c2->SetLeftMargin(0.16);
    c2->SetRightMargin(0.11);
    c2->SetGridx();
    c2->SetGridy();
    c2->SetLogz();

    hLinearResVsMC->Draw("COLZ");

    // Perfect residual reference line (zero residual)
    TF1 *perfectLine2 = new TF1("perfectLine2", "0", xMin, xMax);
    perfectLine2->SetLineColor(kRed);
    perfectLine2->SetLineStyle(2);
    perfectLine2->Draw("SAME");

    // Create and draw TGraphErrors for linear residuals
    TGraphErrors *grLinear = new TGraphErrors(nProfiles, mcBetaValues, linearMean, 0, linearError);
    grLinear->SetMarkerStyle(20);
    grLinear->SetMarkerColor(kBlack);
    grLinear->SetMarkerSize(3.0);
    grLinear->Draw("P");

    c2->Print(outputName);

    // Third page - Comparison of TGraphErrors
    TCanvas *c3 = new TCanvas("c3", "Beta Reconstruction Methods Comparison", 3508, 2480);
    c3->SetLeftMargin(0.16);
    c3->SetGridx();
    c3->SetGridy();

    // Create a frame for the comparison plot
    TH2F *hFrame = new TH2F("hFrame", "Beta Reconstruction Methods Comparison;#beta_{MC};1/#beta_{rec} - 1/#beta_{MC}",
                            100, xMin, xMax, 100, yMinRes, yMaxRes);
    hFrame->Draw();

    // Define fit functions
    TF1 *fNonlinear = new TF1("fNonlinear", "[0] * exp(-[1] * x)", xMin, xMax);
    fNonlinear->SetParameters(0.1, 5.0, 0.0); // Initial parameters
    fNonlinear->SetLineColor(kBlue);
    fNonlinear->SetLineStyle(2);

    TF1 *fLinear = new TF1("fLinear", "[0] * exp(-[1] * x)", xMin, xMax);
    fLinear->SetParameters(0.1, 5.0, 0.0); // Initial parameters
    fLinear->SetLineColor(kRed);
    fLinear->SetLineStyle(2);

    // Create and setup graphs with weights
    grNonlinear->SetMarkerStyle(20);
    grNonlinear->SetMarkerColor(kBlue);
    grNonlinear->SetLineColor(kBlue);
    grNonlinear->SetMarkerSize(3.0);

    grNonlinear->Draw("LEP");
    grLinear->Draw("LEP");

    // Then update errors with weights for fitting only
    TGraphErrors *grNonlinearFit = (TGraphErrors *)grNonlinear->Clone("grNonlinearFit");
    for (int i = 0; i < nProfiles; i++)
        if (nonlinearWeight[i] > 0)
            grNonlinearFit->SetPointError(i, 0, nonlinearError[i] / sqrt(nonlinearWeight[i]));

    grNonlinearFit->Fit(fNonlinear, "NQR"); // Q for quiet, R for using range
    fNonlinear->Draw("SAME");

    grLinear->SetMarkerStyle(21);
    grLinear->SetMarkerColor(kRed);
    grLinear->SetLineColor(kRed);
    grLinear->SetMarkerSize(3.0);

    // First draw points with original error bars
    TGraphErrors *grLinearFit = (TGraphErrors *)grLinear->Clone("grLinearFit");
    for (int i = 0; i < nProfiles; i++)
        if (linearWeight[i] > 0)
            grLinearFit->SetPointError(i, 0, linearError[i] / sqrt(linearWeight[i]));

    grLinearFit->Fit(fLinear, "NQR");
    fLinear->Draw("SAME");

    // Add legend with fit equations
    TLegend *legend = new TLegend(0.45, 0.65, 0.85, 0.85);
    legend->AddEntry(grNonlinear, "Non-linear Method", "lp");
    if (fNonlinear->GetParameter(1) > 0)
        legend->AddEntry(fNonlinear, Form("Fit: %.2fe^{-%.2fx}", fNonlinear->GetParameter(0), fNonlinear->GetParameter(1)), "l");
    else
        legend->AddEntry(fNonlinear, Form("Fit: %.2fe^{%.2fx}", fNonlinear->GetParameter(0), -fNonlinear->GetParameter(1)), "l");
    legend->AddEntry(grLinear, "Linear Method", "lp");
    if (fLinear->GetParameter(1) > 0)
        legend->AddEntry(fLinear, Form("Fit: %.2fe^{-%.2fx}", fLinear->GetParameter(0), fLinear->GetParameter(1)), "l");
    else
        legend->AddEntry(fLinear, Form("Fit: %.2fe^{%.2fx}", fLinear->GetParameter(0), -fLinear->GetParameter(1)), "l");
    legend->SetBorderSize(0);
    legend->Draw();

    // Draw zero line
    TF1 *zeroLine = new TF1("zeroLine", "0", xMin, xMax);
    zeroLine->SetLineStyle(2);
    zeroLine->SetLineColor(kGray + 2);
    zeroLine->Draw("SAME");

    // Add Z value and energyLossScale at the top center if available
    TPaveText *infoText = nullptr;
    if (zValue > 0)
    {
        infoText = new TPaveText(0.2, 0.92, 0.8, 0.98, "NDC");
        infoText->SetFillColor(0);
        infoText->SetBorderSize(0);
        infoText->AddText(Form("Z = %d, Energy Loss Scale = %.1f", zValue, energyLossScale));
        infoText->Draw();
    }

    c3->Print(outputName);
    c3->Print(Form("%s]", outputName)); // Close PDF file

    // Clean up
    if (zValue > 0)
    {
        delete infoText;
    }
    delete zeroLine;

    // Close and delete the file after all references are gone
    file->Close();

    std::cout << "Beta residuals comparison plot saved to: " << outputName << std::endl;
}
