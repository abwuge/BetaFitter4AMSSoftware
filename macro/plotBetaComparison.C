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
        std::cout << "Try to open file: " << fileName << std::endl;

        file = TFile::Open(fileName.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "Error: Unable to find/open file " << fileName << std::endl;
            return;
        }
    }

    // Generate output name
    // ------------------------------------------------------------------------

    // Get Z value and energyLossScale from the file
    auto params = getParamsFromReadme(fileName);
    int zValue = params.first;
    double energyLossScale = params.second;
    std::string actualOutputName;

    // Set output name
    if (!outputName)
    {
        if (zValue > 0 && energyLossScale > 0)
            actualOutputName = Form("test_beta_Z%d_zeta%.3f.pdf", zValue, energyLossScale);
        else
            actualOutputName = "test_beta_comparison.pdf";
        outputName = actualOutputName.c_str();
    }

    // Get tree
    // ------------------------------------------------------------------------

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

    const double xMin = tree->GetMinimum("mcBeta");
    const double xMax = tree->GetMaximum("mcBeta");

    double mcBeta, linearBeta, nonlinearBeta;

    tree->SetBranchAddress("mcBeta", &mcBeta);
    tree->SetBranchAddress("linearBeta", &linearBeta);
    tree->SetBranchAddress("nonlinearBeta", &nonlinearBeta);

    // Bin
    // ------------------------------------------------------------------------

    // Set bin numbers
    const int nBinsX = 40;
    const int nBinsY = 100;

    // Bin width
    const double binWidth = (xMax - xMin) / nBinsX;

    // Get Quantile Range
    // ------------------------------------------------------------------------

    // Helper function to get range using quantiles
    auto getQuantileRange = [&tree](const char *branchName, double lowQuantile = 0.01, double highQuantile = 0.99, bool diff = true) -> std::pair<double, double>
    {
        TH1D hTemp(Form("hTemp_%s", branchName), "", 10000, 0, 0);

        if (diff)
            tree->Draw(Form("1/%s - 1/mcBeta>>hTemp_%s", branchName, branchName), "", "goff");
        else
            tree->Draw(Form("%s>>hTemp_%s", branchName, branchName), "", "goff");

        Double_t xq[2] = {lowQuantile, highQuantile};
        Double_t yq[2] = {0, 0};
        hTemp.GetQuantiles(2, yq, xq);

        std::cout << "Range for " << branchName << " residuals using quantiles: [" << yq[0] << ", " << yq[1] << "]" << std::endl;
        return std::make_pair(yq[0], yq[1]);
    };

    // Get beta range
    const auto betaRange = getQuantileRange("nonlinearBeta", 1e-4, 1 - 1e-4, false);
    const double yMin = betaRange.first;
    const double yMax = betaRange.second;

    // Get residuals range
    const auto resRange = getQuantileRange("nonlinearBeta", 1e-4, 1 - 1e-4);
    const double yMinRes = resRange.first;
    const double yMaxRes = resRange.second;

    // Draw
    // ------------------------------------------------------------------------

    // Create beta histograms
    TH2F *hNonlinearVsMC = new TH2F("hNonlinearVsMC",
                                    "Non-linear Beta vs MC Beta;#beta_{MC};#beta_{non-linear}",
                                    100, xMin, xMax, 100, yMin, yMax);
    TH2F *hLinearVsMC = new TH2F("hLinearVsMC",
                                 "Linear Beta vs MC Beta;#beta_{MC};#beta_{linear}",
                                 100, xMin, xMax, 100, yMin, yMax);

    // Create residuals histograms
    TH2F *hNonlinearResVsMC = new TH2F("hNonlinearResVsMC",
                                       "Non-linear Reconstruction Residuals;#beta_{MC};1/#beta_{non-linear} - 1/#beta_{MC}",
                                       nBinsX, xMin, xMax, nBinsY, yMinRes, yMaxRes);
    TH2F *hLinearResVsMC = new TH2F("hLinearResVsMC",
                                    "Linear Reconstruction Residuals;#beta_{MC};1/#beta_{linear} - 1/#beta_{MC}",
                                    nBinsX, xMin, xMax, nBinsY, yMinRes, yMaxRes);

    // Create comparison plot
    TH2F *hComparison = new TH2F("hComparison",
                                 "Beta Reconstruction Methods Comparison;#beta_{MC};1/#beta_{rec} - 1/#beta_{MC}",
                                 nBinsX, xMin, xMax, nBinsY, yMinRes, yMaxRes);

    // Set minimum value to 1
    hNonlinearVsMC->SetMinimum(1);
    hLinearVsMC->SetMinimum(1);
    hNonlinearResVsMC->SetMinimum(1);
    hLinearResVsMC->SetMinimum(1);

    // Fill histograms
    tree->Draw("nonlinearBeta:mcBeta>>hNonlinearVsMC",
               Form("nonlinearBeta > %f && nonlinearBeta < %f", yMin, yMax),
               "goff");
    tree->Draw("linearBeta:mcBeta>>hLinearVsMC",
               Form("linearBeta > %f && linearBeta < %f", yMin, yMax),
               "goff");
    tree->Draw("1/nonlinearBeta - 1/mcBeta:mcBeta>>hNonlinearResVsMC",
               Form("1/nonlinearBeta - 1/mcBeta > %f && 1/nonlinearBeta - 1/mcBeta < %f", yMinRes, yMaxRes),
               "goff");
    tree->Draw("1/linearBeta - 1/mcBeta:mcBeta>>hLinearResVsMC",
               Form("1/linearBeta - 1/mcBeta > %f && 1/linearBeta - 1/mcBeta < %f", yMinRes, yMaxRes),
               "goff");

    // Print beta
    // ------------------------------------------------------------------------

    // Create canvas
    TCanvas *canvas = new TCanvas("canvas", "", 3508, 2480); // A4 landscape
    canvas->SetLeftMargin(0.16);
    canvas->SetRightMargin(0.12);
    canvas->SetGridx();
    canvas->SetGridy();
    canvas->SetLogz();

    // Open PDF file
    canvas->Print(Form("%s[", outputName));

    // Create perfect correlation line
    TF1 *perfectCorrelation = new TF1("perfectCorrelation", "x", xMin, xMax);
    perfectCorrelation->SetLineColor(kRed);
    perfectCorrelation->SetLineStyle(2);

    // Create info text
    TPaveText *infoText = nullptr;
    if (zValue > 0)
    {
        infoText = new TPaveText(0.2, 0.92, 0.8, 0.98, "NDC");
        infoText->SetFillColor(0);
        infoText->SetBorderSize(0);
        infoText->AddText(Form("Z = %d, #zeta = %.3f", zValue, energyLossScale));
    }

    // Draw nonlinear beta vs MC beta
    hNonlinearVsMC->Draw("COLZ");
    perfectCorrelation->Draw("SAME");
    if (infoText)
        infoText->Draw();
    canvas->Print(outputName);

    // Draw linear beta vs MC beta
    hLinearVsMC->Draw("COLZ");
    perfectCorrelation->Draw("SAME");
    if (infoText)
        infoText->Draw();
    canvas->Print(outputName);

    // Fit
    // ------------------------------------------------------------------------

    // Vectors to store fit results
    std::vector<double> nonlinearMCBetaValues, linearMCBetaValues;
    std::vector<double> nonlinearMeans, nonlinearErrors;
    std::vector<double> linearMeans, linearErrors;

    // Prepare vectors
    nonlinearMCBetaValues.reserve(nBinsX);
    nonlinearMeans.reserve(nBinsX);
    nonlinearErrors.reserve(nBinsX);
    linearMCBetaValues.reserve(nBinsX);
    linearMeans.reserve(nBinsX);
    linearErrors.reserve(nBinsX);

    // Alias pointers for easier access
    TH2F *hRes[2] = {hNonlinearResVsMC, hLinearResVsMC};
    std::vector<double> *mcBetaValues[2] = {&nonlinearMCBetaValues, &linearMCBetaValues};
    std::vector<double> *means[2] = {&nonlinearMeans, &linearMeans};
    std::vector<double> *errors[2] = {&nonlinearErrors, &linearErrors};

    // Create fit info text
    TPaveText *fitInfoText = new TPaveText(0.3, 0.92, 0.7, 0.98, "NDC");
    fitInfoText->SetFillColor(0);
    fitInfoText->SetBorderSize(0);

    // Create fit function
    TF1 *fGaus = new TF1("fGaus", "gaus", yMinRes, yMaxRes);
    fGaus->SetParLimits(1, yMinRes, yMaxRes);
    fGaus->SetNpx(1000);

    // Enable log y
    canvas->SetLogy(1);

    // Fit each column
    for (int i = 0; i < nBinsX; ++i)
    {
        double binCenter = xMin + (i + 0.5) * binWidth;
        bool isValid = false;

        // Process each residual histogram
        for (int j = 0; j < 2; ++j)
        {
            TH1D *proj = hRes[j]->ProjectionY(Form("proj_%d_%d", j, i), i + 1, i + 1);
            Int_t projEntries = proj->GetEntries();

            fitInfoText->Clear();
            proj->Draw();

            if (projEntries > 50)
            {
                fGaus->SetParameters(proj->GetMaximum(), proj->GetMean(), proj->GetRMS());
                TFitResultPtr fitResult = proj->Fit(fGaus, "SQNR");

                if (fitResult->Status() == 0)
                {
                    isValid = true;
                    mcBetaValues[j]->push_back(binCenter);
                    means[j]->push_back(fitResult->Parameter(1));
                    errors[j]->push_back(fitResult->Parameter(2));

                    fGaus->Draw("SAME");
                    fitInfoText->AddText(
                        Form("%s #beta: %.3f #mu: %.3f #sigma: %.3f",
                             j ? "Linear" : "Nonlinear", binCenter,
                             fitResult->Parameter(1), fitResult->Parameter(2)));
                }
                else
                    fitInfoText->AddText(
                        Form("%s #beta: %.3f FIT FAILED",
                             j ? "Linear" : "Nonlinear", binCenter));
            }
            else
                fitInfoText->AddText(
                    Form("%s #beta: %.3f Entries (%d) too low",
                         j ? "Linear" : "Nonlinear", binCenter, projEntries));

            fitInfoText->Draw();
            canvas->Print(outputName);

            delete proj;
        }
    }

    // Disable log y
    canvas->SetLogy(0);

    // Print residuals
    // ------------------------------------------------------------------------

    // Create TGraphErrors for residuals
    TGraphErrors *grNonlinear = new TGraphErrors(nonlinearMCBetaValues.size(),
                                                 nonlinearMCBetaValues.data(),
                                                 nonlinearMeans.data(),
                                                 nullptr,
                                                 nonlinearErrors.data());
    grNonlinear->SetMarkerStyle(20);
    grNonlinear->SetMarkerSize(3.0);

    TGraphErrors *grLinear = new TGraphErrors(linearMCBetaValues.size(),
                                              linearMCBetaValues.data(),
                                              linearMeans.data(),
                                              nullptr,
                                              linearErrors.data());
    grLinear->SetMarkerStyle(20);
    grLinear->SetMarkerSize(3.0);

    // Create fit functions & Fit
    TF1 *fNonlinear = new TF1("fNonlinear", "[0] * exp(-[1] * x)", xMin, xMax);
    fNonlinear->SetParameters(0.1, 5.0, 0.0);
    fNonlinear->SetLineColor(kRed);
    fNonlinear->SetLineStyle(2);
    grNonlinear->Fit(fNonlinear, "QNR");

    TF1 *fLinear = new TF1("fLinear", "[0] * exp(-[1] * x)", xMin, xMax);
    fLinear->SetParameters(0.1, 5.0, 0.0);
    fLinear->SetLineColor(kBlue);
    fLinear->SetLineStyle(2);
    grLinear->Fit(fLinear, "QNR");

    // Create perfect residual reference line (zero residual)
    TF1 *perfectResidual = new TF1("perfectResidual", "0", xMin, xMax);
    perfectResidual->SetLineColor(kRed);
    perfectResidual->SetLineStyle(2);

    // Draw nonlinear residuals vs MC beta
    hNonlinearResVsMC->Draw("COLZ");
    perfectResidual->Draw("SAME");
    grNonlinear->Draw("P");
    if (infoText)
        infoText->Draw();
    canvas->Print(outputName);

    // Draw linear residuals vs MC beta
    hLinearResVsMC->Draw("COLZ");
    perfectResidual->Draw("SAME");
    grLinear->Draw("P");
    if (infoText)
        infoText->Draw();
    canvas->Print(outputName);

    // Draw comparison plot
    hComparison->Draw("COLZ");
    grNonlinear->SetMarkerColor(kRed);
    grNonlinear->SetLineColor(kRed);
    grNonlinear->Draw("LEP");
    fNonlinear->Draw("SAME");

    grLinear->SetMarkerColor(kBlue);
    grLinear->SetLineColor(kBlue);
    grLinear->Draw("LEP");
    perfectResidual->Draw("SAME");
    fLinear->Draw("SAME");

    if (infoText)
        infoText->Draw();

    // Add legend
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

    canvas->Print(outputName);

    // Close PDF file
    canvas->Print(Form("%s]", outputName));

    std::cout << "Beta residuals comparison plot saved to: " << outputName << std::endl;
}
