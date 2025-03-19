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

void plotRealDataBetaComparison(std::string fileName = "test.root",
                                const char *restrict = "",
                                const char *outputName = "test_real_beta_comparison.pdf")
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

    const double rigidityMin = tree->GetMinimum("innerRigidity");
    const double rigidityMax = tree->GetMaximum("innerRigidity");
    const double betaMin = tree->GetMinimum("betaRigidity");
    const double betaMax = tree->GetMaximum("betaRigidity");

    const char *baseBetaBranch = "betaRigidity";
    const char *baseBetaTitle = "#beta_{rigidity}";
    const char *rigidityAxisBranch = "innerRigidity";
    const char *rigidityAxisTitle = "Rigidity [GV]";
    const char *betaAxisBranch = "betaRigidity";
    const char *betaAxisTitle = "#beta_{rigidity}";

    double baseBeta, rigidityAxis, betaAxis;
    double linearBeta, nonlinearBeta;

    tree->SetBranchAddress(baseBetaBranch, &baseBeta);
    tree->SetBranchAddress(rigidityAxisBranch, &rigidityAxis);
    tree->SetBranchAddress(betaAxisBranch, &betaAxis);
    tree->SetBranchAddress("linearBeta", &linearBeta);
    tree->SetBranchAddress("nonlinearBeta", &nonlinearBeta);

    // Bin
    // ------------------------------------------------------------------------

    // Set bin numbers
    const int nBinsX = 40;
    const int nBinsY = 100;

    // Linear bins in x-axis (beta)
    const double betaBinWidth = (betaMax - betaMin) / nBinsX;

    // Logarithmic bins in x-axis (rigidity)
    double rigidityBins[nBinsX + 1]{};
    const double logRigidityMin = TMath::Log10(rigidityMin);
    const double logRigidityMax = TMath::Log10(rigidityMax);
    const double rigidityBinWidth = (logRigidityMax - logRigidityMin) / nBinsX;
    for (int i = 0; i <= nBinsX; i++)
        rigidityBins[i] = TMath::Power(10, logRigidityMin + i * rigidityBinWidth);

    // Get Quantile Range
    // ------------------------------------------------------------------------

    // Helper function to get range using quantiles
    auto getQuantileRange = [&](const char *branchName, double lowQuantile = 0.01, double highQuantile = 0.99, bool diff = true) -> std::pair<double, double>
    {
        TH1D hTemp(Form("hTemp_%s", branchName), "", 10000, 0, 0);

        if (diff)
            tree->Draw(Form("1/%s - 1/%s>>hTemp_%s", branchName, baseBetaBranch, branchName), restrict, "goff");
        else
            tree->Draw(Form("%s>>hTemp_%s", branchName, branchName), restrict, "goff");

        Double_t xq[2] = {lowQuantile, highQuantile};
        Double_t yq[2] = {0, 0};
        hTemp.GetQuantiles(2, yq, xq);

        if (diff)
            std::cout << "Range for " << branchName << " residuals using quantiles: [" << yq[0] << ", " << yq[1] << "]" << std::endl;
        else
            std::cout << "Range for " << branchName << " using quantiles: [" << yq[0] << ", " << yq[1] << "]" << std::endl;

        return std::make_pair(yq[0], yq[1]);
    };

    // Get beta range
    const auto nonlinearBetaRange = getQuantileRange("nonlinearBeta", 1e-4, 1 - 1e-4, false);
    const double yMinNL = nonlinearBetaRange.first;
    const double yMaxNL = nonlinearBetaRange.second;

    const auto linearBetaRange = getQuantileRange("linearBeta", 1e-4, 0.99, false);
    const double yMinL = linearBetaRange.first;
    const double yMaxL = linearBetaRange.second;

    // Get residuals range
    const auto nonlinearResRange = getQuantileRange("nonlinearBeta", 1e-4, 1 - 1e-4);
    const double yMinResNL = nonlinearResRange.first;
    const double yMaxResNL = nonlinearResRange.second;

    const auto linearResRange = getQuantileRange("linearBeta", 0.001, 1 - 1e-6);
    const double yMinResL = linearResRange.first;
    const double yMaxResL = linearResRange.second;

    const double yMinRes = std::min(nonlinearResRange.first, linearResRange.first);
    const double yMaxRes = std::max(nonlinearResRange.second, linearResRange.second);

    // Draw
    // ------------------------------------------------------------------------

    // Create beta histograms
    TH2F *hNonlinearVsRigidity = new TH2F("hNonlinearVsRigidity",
                                          Form(";%s;#beta_{non-linear}", rigidityAxisTitle),
                                          nBinsX, rigidityBins, 100, yMinNL, yMaxNL);
    TH2F *hLinearVsRigidity = new TH2F("hLinearVsRigidity",
                                       Form(";%s;#beta_{linear}", rigidityAxisTitle),
                                       nBinsX, rigidityBins, 100, yMinL, yMaxL);
    TH2F *hNonlinearVsBeta = new TH2F("hNonlinearVsBeta",
                                      Form(";%s;#beta_{non-linear}", betaAxisTitle),
                                      nBinsX, betaMin, betaMax, 100, yMinNL, yMaxNL);
    TH2F *hLinearVsBeta = new TH2F("hLinearVsBeta",
                                   Form(";%s;#beta_{linear}", betaAxisTitle),
                                   nBinsX, betaMin, betaMax, 100, yMinL, yMaxL);

    // Create Residuals histograms
    TH2F *hNonlinearResVsRigidity = new TH2F("hNonlinearResVsRigidity",
                                             Form(";%s;1/#beta_{non-linear} - 1/%s",
                                                  rigidityAxisTitle, baseBetaTitle),
                                             nBinsX, rigidityBins, nBinsY, yMinResNL, yMaxResNL);
    TH2F *hLinearResVsRigidity = new TH2F("hLinearResVsRigidity",
                                          Form(";%s;1/#beta_{linear} - 1/%s",
                                               rigidityAxisTitle, baseBetaTitle),
                                          nBinsX, rigidityBins, nBinsY, yMinResL, yMaxResL);
    TH2F *hNonlinearResVsBeta = new TH2F("hNonlinearResVsBeta",
                                         Form(";%s;1/#beta_{non-linear} - 1/%s",
                                              betaAxisTitle, baseBetaTitle),
                                         nBinsX, betaMin, betaMax, nBinsY, yMinResNL, yMaxResNL);
    TH2F *hLinearResVsBeta = new TH2F("hLinearResVsBeta",
                                      Form(";%s;1/#beta_{linear} - 1/%s",
                                           betaAxisTitle, baseBetaTitle),
                                      nBinsX, betaMin, betaMax, nBinsY, yMinResNL, yMaxResNL);

    // Create comparison plot
    TH2F *hComparisonRigidity = new TH2F("hComparisonRigidity",
                                         Form(";%s;#beta_{rec}", rigidityAxisTitle),
                                         nBinsX, rigidityBins, nBinsY, yMinResNL, yMaxResNL);
    TH2F *hComparisonBeta = new TH2F("hComparisonBeta",
                                     Form(";%s;#beta_{rec}", betaAxisTitle),
                                     nBinsX, betaMin, betaMax, nBinsY, yMinResNL, yMaxResNL);

    // Set minimum value to 1
    hNonlinearVsRigidity->SetMinimum(1);
    hLinearVsRigidity->SetMinimum(1);
    hNonlinearVsBeta->SetMinimum(1);
    hLinearVsBeta->SetMinimum(1);
    hNonlinearResVsRigidity->SetMinimum(1);
    hLinearResVsRigidity->SetMinimum(1);
    hNonlinearResVsBeta->SetMinimum(1);
    hLinearResVsBeta->SetMinimum(1);

    // Fill histograms
    tree->Draw(Form("nonlinearBeta:%s>>hNonlinearVsRigidity", rigidityAxisBranch), restrict, "goff");
    tree->Draw(Form("linearBeta:%s>>hLinearVsRigidity", rigidityAxisBranch), restrict, "goff");
    tree->Draw(Form("nonlinearBeta:%s>>hNonlinearVsBeta", betaAxisBranch), restrict, "goff");
    tree->Draw(Form("linearBeta:%s>>hLinearVsBeta", betaAxisBranch), restrict, "goff");
    tree->Draw(Form("1/nonlinearBeta - 1/%s:%s>>hNonlinearResVsRigidity", baseBetaBranch, rigidityAxisBranch), restrict, "goff");
    tree->Draw(Form("1/linearBeta - 1/%s:%s>>hLinearResVsRigidity", baseBetaBranch, rigidityAxisBranch), restrict, "goff");
    tree->Draw(Form("1/nonlinearBeta - 1/%s:%s>>hNonlinearResVsBeta", baseBetaBranch, betaAxisBranch), restrict, "goff");
    tree->Draw(Form("1/linearBeta - 1/%s:%s>>hLinearResVsBeta", baseBetaBranch, betaAxisBranch), restrict, "goff");

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

    // Open log(x)
    canvas->SetLogx(1);

    // Draw nonlinear beta vs rigidity
    hNonlinearVsRigidity->Draw("COLZ");
    canvas->Print(outputName);

    // Draw linear beta vs rigidity
    hLinearVsRigidity->Draw("COLZ");
    canvas->Print(outputName);

    // Close log(x)
    canvas->SetLogx(0);

    // Draw nonlinear beta vs beta
    hNonlinearVsBeta->Draw("COLZ");
    canvas->Print(outputName);

    // Draw linear beta vs beta
    hLinearVsBeta->Draw("COLZ");
    canvas->Print(outputName);

    // Fit
    // ------------------------------------------------------------------------

    // Arrays to store fit results
    double betaValues[nBinsX]{}, rigidityValues[nBinsX]{};
    double rigidityNLMean[nBinsX]{}, rigidityNLError[nBinsX]{};
    double rigidityLMean[nBinsX]{}, rigidityLError[nBinsX]{};
    double betaNLMean[nBinsX]{}, betaNLError[nBinsX]{};
    double betaLMean[nBinsX]{}, betaLError[nBinsX]{};

    // Alias pointers for easier access
    TH2F *hRes[4] = {hNonlinearResVsRigidity, hLinearResVsRigidity, hNonlinearResVsBeta, hLinearResVsBeta};
    double *means[4] = {rigidityNLMean, rigidityLMean, betaNLMean, betaLMean};
    double *errors[4] = {rigidityNLError, rigidityLError, betaNLError, betaLError};

    // Fit each column
    TF1 *fGaus = new TF1("fGaus", "gaus", yMinRes, yMaxRes);
    fGaus->SetParLimits(1, yMinRes, yMaxRes);
    for (int i = 0; i < nBinsX; ++i)
    {
        double rigidityLower = TMath::Log10(rigidityBins[i]);
        double rigidityUpper = TMath::Log10(rigidityBins[i + 1]);
        double rigidityCenter = (rigidityLower + rigidityUpper) / 2;
        double rigidityBinCenter = TMath::Power(10, rigidityCenter);
        double betaBinCenter = betaMin + (i + 0.5) * betaBinWidth;

        rigidityValues[i] = rigidityBinCenter;
        betaValues[i] = betaBinCenter;

        // Process each residual histogram
        for (int j = 0; j < 4; ++j)
        {
            TH1D *proj = hRes[j]->ProjectionY(Form("proj_%d_%d", j, i), i + 1, i + 1);
            Int_t projEntries = proj->GetEntries();

            if (projEntries > 10)
            {
                fGaus->SetParameters(proj->GetMaximum(), proj->GetMean(), proj->GetRMS());
                TFitResultPtr fitResult = proj->Fit(fGaus, "QS");

                if (fitResult->Status() == 0)
                {
                    means[j][i] = fitResult->Parameter(1);
                    errors[j][i] = fitResult->Parameter(2);
                }
                else
                {
                    means[j][i] = proj->GetMean();
                    errors[j][i] = proj->GetRMS();
                }
            }

            // proj->Draw();
            // canvas->Print(outputName);

            delete proj;
        }
    }

    // Print residuals
    // ------------------------------------------------------------------------

    // Create TGraphErrors for residuals
    TGraphErrors *grNLRigidity = new TGraphErrors(nBinsX, rigidityValues, rigidityNLMean, 0, rigidityNLError);
    grNLRigidity->SetMarkerStyle(20);
    grNLRigidity->SetMarkerSize(3.0);

    TGraphErrors *grLRigidity = new TGraphErrors(nBinsX, rigidityValues, rigidityLMean, 0, rigidityLError);
    grLRigidity->SetMarkerStyle(20);
    grLRigidity->SetMarkerSize(3.0);

    TGraphErrors *grNLBeta = new TGraphErrors(nBinsX, betaValues, betaNLMean, 0, betaNLError);
    grNLBeta->SetMarkerStyle(20);
    grNLBeta->SetMarkerSize(3.0);

    TGraphErrors *grLBeta = new TGraphErrors(nBinsX, betaValues, betaLMean, 0, betaLError);
    grLBeta->SetMarkerStyle(20);
    grLBeta->SetMarkerSize(3.0);

    // Create perfect residual reference line (zero residual)
    TF1 *perfectRigidityResidualLine = new TF1("perfectRigidityResidualLine", "0", rigidityMin, rigidityMax);
    perfectRigidityResidualLine->SetLineColor(kRed);
    perfectRigidityResidualLine->SetLineStyle(2);

    TF1 *perfectBetaResidualLine = new TF1("perfectBetaResidualLine", "0", betaMin, betaMax);
    perfectBetaResidualLine->SetLineColor(kRed);
    perfectBetaResidualLine->SetLineStyle(2);

    // Open log(x)
    canvas->SetLogx(1);

    // Draw nonlinear residuals vs rigidity
    hNonlinearResVsRigidity->Draw("COLZ");
    perfectRigidityResidualLine->Draw("SAME");
    grNLRigidity->Draw("P");
    canvas->Print(outputName);

    // Draw linear residuals vs rigidity
    hLinearResVsRigidity->Draw("COLZ");
    perfectRigidityResidualLine->Draw("SAME");
    grLRigidity->Draw("P");
    canvas->Print(outputName);

    // Close log(x)
    canvas->SetLogx(0);

    // Draw nonlinear residuals vs beta
    hNonlinearResVsBeta->Draw("COLZ");
    perfectBetaResidualLine->Draw("SAME");
    grNLBeta->Draw("P");
    canvas->Print(outputName);

    // Draw linear residuals vs beta
    hLinearResVsBeta->Draw("COLZ");
    perfectBetaResidualLine->Draw("SAME");
    grLBeta->Draw("P");
    canvas->Print(outputName);

    // Draw comparison plot
    canvas->SetLogx(1);
    hComparisonRigidity->Draw("COLZ");
    grNLRigidity->SetMarkerColor(kRed);
    grNLRigidity->SetLineColor(kRed);
    grNLRigidity->Draw("LEP");
    grLRigidity->SetMarkerColor(kBlue);
    grLRigidity->SetLineColor(kBlue);
    grLRigidity->Draw("LEP");
    perfectRigidityResidualLine->Draw("SAME");
    canvas->Print(outputName);

    canvas->SetLogx(0);
    hComparisonBeta->Draw("COLZ");
    grNLBeta->SetMarkerColor(kRed);
    grNLBeta->SetLineColor(kRed);
    grNLBeta->Draw("LEP");
    grLBeta->SetMarkerColor(kBlue);
    grLBeta->SetLineColor(kBlue);
    grLBeta->Draw("LEP");
    perfectBetaResidualLine->Draw("SAME");
    canvas->Print(outputName);

    // Close PDF file
    canvas->Print(Form("%s]", outputName));

    std::cout << "Real data beta residuals comparison plot saved to: " << outputName << std::endl;
}
