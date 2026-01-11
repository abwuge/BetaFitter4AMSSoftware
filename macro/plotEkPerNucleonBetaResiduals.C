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
 * separated by NaF and AgL regions
 * 
 * @param fileName Input ROOT file path
 * @param nucleonNumber Number of nucleons per nucleus, default is 4
 * @param outputName Output PDF file name
 */
void plotEkPerNucleonBetaResiduals(
    std::string fileName = "test.root",
    int nucleonNumber = 4,
    const char *outputName = "test_ek_per_nucleon_beta_residuals.pdf")
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
    
    // NaF histograms (inner region)
    TH2F *hLinearResVsLogEk_NaF = new TH2F("hLinearResVsLogEk_NaF", 
                                          "NaF: Linear Beta Residuals vs. E_{k}/n;E_{k}/n [GeV/nucleon];1/#beta_{linear} - 1/#beta_{MC}",
                                          nBinsX, logEkBins, nBinsY, yMinRes, yMaxRes);
    
    TH2F *hNonlinearResVsLogEk_NaF = new TH2F("hNonlinearResVsLogEk_NaF", 
                                             "NaF: Non-linear Beta Residuals vs. E_{k}/n;E_{k}/n [GeV/nucleon];1/#beta_{non-linear} - 1/#beta_{MC}",
                                             nBinsX, logEkBins, nBinsY, yMinRes, yMaxRes);
    
    // AgL histograms (outer region)
    TH2F *hLinearResVsLogEk_AgL = new TH2F("hLinearResVsLogEk_AgL", 
                                          "AgL: Linear Beta Residuals vs. E_{k}/n;E_{k}/n [GeV/nucleon];1/#beta_{linear} - 1/#beta_{MC}",
                                          nBinsX, logEkBins, nBinsY, yMinRes, yMaxRes);
    
    TH2F *hNonlinearResVsLogEk_AgL = new TH2F("hNonlinearResVsLogEk_AgL", 
                                             "AgL: Non-linear Beta Residuals vs. E_{k}/n;E_{k}/n [GeV/nucleon];1/#beta_{non-linear} - 1/#beta_{MC}",
                                             nBinsX, logEkBins, nBinsY, yMinRes, yMaxRes);
    
    TH2F *hComparisonVsLogEk = new TH2F("hComparisonVsLogEk", 
                                       "Beta Residuals Comparison vs. E_{k}/n;E_{k}/n [GeV/nucleon];1/#beta_{rec} - 1/#beta_{MC}",
                                       nBinsX, logEkBins, nBinsY, yMinRes, yMaxRes);
    
    // Set minimum value to 1 for log scale visibility
    hLinearResVsLogEk_NaF->SetMinimum(1);
    hNonlinearResVsLogEk_NaF->SetMinimum(1);
    hLinearResVsLogEk_AgL->SetMinimum(1);
    hNonlinearResVsLogEk_AgL->SetMinimum(1);
    
    // Fill histograms
    // ------------------------------------------------------------------------

    for (int i = 0; i < tree->GetEntries(); ++i)
    {
        tree->GetEntry(i);
        
        // Skip invalid data
        if (mevmom1[2] == -1000 || mcBeta[2] == -1000)
            continue;
            
        // Skip invalid coordinates
        if (mevcoo1[17][0] == -1000 || mevcoo1[17][1] == -1000 || mevcoo1[17][2] == -1000)
            continue;
        
        // Calculate kinetic energy Ek = mc^2 * (gamma - 1)
        double momentum = mevmom1[2]; // GeV/c
        double mass = momentum / mcBeta[2] * sqrt(1 - mcBeta[2] * mcBeta[2]); // GeV/c^2
        double gamma = 1.0 / sqrt(1 - mcBeta[2] * mcBeta[2]);
        double kineticEnergy = mass * (gamma - 1.0); // GeV
        
        // Calculate kinetic energy per nucleon and take log
        double ekPerNucleon = kineticEnergy / nucleonNumber; // GeV/nucleon
        double logEkPerNucleon = TMath::Log10(ekPerNucleon);
        
        // Calculate residuals
        double linearRes = 1.0/linearBeta - 1.0/mcBeta[17];
        double nonlinearRes = 1.0/nonlinearBeta - 1.0/mcBeta[17];
        
        // Determine if the hit is in AgL (outer region) or NaF (inner region)
        bool isAgL = (TMath::Abs(mevcoo1[17][0]) > 34 || TMath::Abs(mevcoo1[17][1]) > 34);
        
        // Fill appropriate histograms based on region
        if (isAgL) {
            hLinearResVsLogEk_AgL->Fill(logEkPerNucleon, linearRes);
            hNonlinearResVsLogEk_AgL->Fill(logEkPerNucleon, nonlinearRes);
        } else {
            hLinearResVsLogEk_NaF->Fill(logEkPerNucleon, linearRes);
            hNonlinearResVsLogEk_NaF->Fill(logEkPerNucleon, nonlinearRes);
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
    TH2F* histsToUpdate[] = {
        hLinearResVsLogEk_NaF, hNonlinearResVsLogEk_NaF,
        hLinearResVsLogEk_AgL, hNonlinearResVsLogEk_AgL,
        hComparisonVsLogEk
    };
    
    for (auto hist : histsToUpdate) {
        TAxis* xAxis = hist->GetXaxis();
        xAxis->SetMoreLogLabels();
        xAxis->SetNoExponent();
        
        // Convert bin labels from log to actual values
        for (int i = 1; i <= nBinsX; i++) {
            double logValue = xAxis->GetBinCenter(i);
            double actualValue = pow(10, logValue);
            
            // Only set labels for some bins to avoid overcrowding
            if (i % 5 == 0 || i == 1 || i == nBinsX) {
                if (actualValue < 0.1) {
                    xAxis->SetBinLabel(i, Form("%.3f", actualValue));
                } else if (actualValue < 1) {
                    xAxis->SetBinLabel(i, Form("%.2f", actualValue));
                } else {
                    xAxis->SetBinLabel(i, Form("%.1f", actualValue));
                }
            }
        }
    }
    
    // Draw residual plots
    // ------------------------------------------------------------------------

    // NaF plots
    hLinearResVsLogEk_NaF->Draw("COLZ");
    infoText->Draw();
    canvas->Print(outputName);
    
    hNonlinearResVsLogEk_NaF->Draw("COLZ");
    infoText->Draw();
    canvas->Print(outputName);
    
    // AgL plots
    hLinearResVsLogEk_AgL->Draw("COLZ");
    infoText->Draw();
    canvas->Print(outputName);
    
    hNonlinearResVsLogEk_AgL->Draw("COLZ");
    infoText->Draw();
    canvas->Print(outputName);
    
    // Fit
    // ------------------------------------------------------------------------

    // Vectors to store fit results
    std::vector<double> logEkValues;
    std::vector<double> ekValues; // Actual Ek values for display
    std::vector<double> linearMeans_NaF, linearErrors_NaF;
    std::vector<double> nonlinearMeans_NaF, nonlinearErrors_NaF;
    std::vector<double> linearMeans_AgL, linearErrors_AgL;
    std::vector<double> nonlinearMeans_AgL, nonlinearErrors_AgL;
    
    logEkValues.reserve(nBinsX);
    ekValues.reserve(nBinsX);
    linearMeans_NaF.reserve(nBinsX);
    linearErrors_NaF.reserve(nBinsX);
    nonlinearMeans_NaF.reserve(nBinsX);
    nonlinearErrors_NaF.reserve(nBinsX);
    linearMeans_AgL.reserve(nBinsX);
    linearErrors_AgL.reserve(nBinsX);
    nonlinearMeans_AgL.reserve(nBinsX);
    nonlinearErrors_AgL.reserve(nBinsX);
    
    // Create fit function
    TF1 *fGaus = new TF1("fGaus", "gaus", yMinRes, yMaxRes);
    fGaus->SetParLimits(1, yMinRes, yMaxRes);
    fGaus->SetNpx(1000);
    
    // Enable log y
    canvas->SetLogy(1);
    
    // Create fit info text
    TPaveText *fitInfoText = new TPaveText(0.3, 0.92, 0.7, 0.98, "NDC");
    fitInfoText->SetFillColor(0);
    fitInfoText->SetBorderSize(0);
    
    // Arrays of histograms to process
    TH2F* histograms[4] = {
        hLinearResVsLogEk_NaF, hNonlinearResVsLogEk_NaF,
        hLinearResVsLogEk_AgL, hNonlinearResVsLogEk_AgL
    };
    
    const char* regionNames[4] = {
        "NaF Linear", "NaF Nonlinear",
        "AgL Linear", "AgL Nonlinear"
    };
    
    // Fit each column for each histogram
    for (int histIdx = 0; histIdx < 4; ++histIdx) {
        TH2F* currentHist = histograms[histIdx];
        const char* regionName = regionNames[histIdx];
        
        for (int i = 0; i < nBinsX; ++i) {
            double logEkLow = logEkBins[i];
            double logEkHigh = logEkBins[i+1];
            double binCenter = (logEkLow + logEkHigh) / 2.0;
            double ekValue = pow(10, binCenter); // Convert to actual Ek value
            
            TH1D *proj = currentHist->ProjectionY(Form("proj_%d_%d", histIdx, i), i + 1, i + 1);
            Int_t projEntries = proj->GetEntries();
            
            fitInfoText->Clear();
            proj->SetTitle(Form("%s: E_{k}/n = %.3f GeV", regionName, ekValue));
            proj->Draw();
            
            if (projEntries > 50) {
                fGaus->SetParameters(proj->GetMaximum(), proj->GetMean(), proj->GetRMS());
                TFitResultPtr fitResult = proj->Fit(fGaus, "SQNR");
                
                if (fitResult->Status() == 0 || fitResult->IsEmpty()) {
                    // Store fit results in appropriate vectors
                    if (histIdx == 0) { // NaF Linear
                        if (std::find(logEkValues.begin(), logEkValues.end(), binCenter) == logEkValues.end()) {
                            logEkValues.push_back(binCenter);
                            ekValues.push_back(ekValue);
                        }
                        linearMeans_NaF.push_back(fitResult->Parameter(1));
                        linearErrors_NaF.push_back(fitResult->Parameter(2));
                    } else if (histIdx == 1) { // NaF Non-linear
                        nonlinearMeans_NaF.push_back(fitResult->Parameter(1));
                        nonlinearErrors_NaF.push_back(fitResult->Parameter(2));
                    } else if (histIdx == 2) { // AgL Linear
                        linearMeans_AgL.push_back(fitResult->Parameter(1));
                        linearErrors_AgL.push_back(fitResult->Parameter(2));
                    } else { // AgL Non-linear
                        nonlinearMeans_AgL.push_back(fitResult->Parameter(1));
                        nonlinearErrors_AgL.push_back(fitResult->Parameter(2));
                    }
                    
                    fGaus->Draw("SAME");
                    fitInfoText->AddText(
                        Form("#mu: %.3f #sigma: %.3f", 
                             fitResult->Parameter(1), 
                             fitResult->Parameter(2)));
                } else {
                    fitInfoText->AddText("Fit failed");
                }
            } else {
                fitInfoText->AddText(
                    Form("Entries (%d) too few", projEntries));
            }
            
            fitInfoText->Draw();
            canvas->Print(outputName);
            
            delete proj;
        }
    }
    
    // Disable log y
    canvas->SetLogy(0);
    
    // Create separate mean and error graphs
    // ------------------------------------------------------------------------
    // Means graphs
    TGraph *grLinearMean_NaF = new TGraph(ekValues.size(), ekValues.data(), linearMeans_NaF.data());
    grLinearMean_NaF->SetMarkerStyle(20);
    grLinearMean_NaF->SetMarkerSize(2.0);
    grLinearMean_NaF->SetMarkerColor(kBlue);
    grLinearMean_NaF->SetLineColor(kBlue);
    grLinearMean_NaF->SetLineWidth(2);
    
    TGraph *grNonlinearMean_NaF = new TGraph(ekValues.size(), ekValues.data(), nonlinearMeans_NaF.data());
    grNonlinearMean_NaF->SetMarkerStyle(20);
    grNonlinearMean_NaF->SetMarkerSize(2.0);
    grNonlinearMean_NaF->SetMarkerColor(kGreen+2);
    grNonlinearMean_NaF->SetLineColor(kGreen+2);
    grNonlinearMean_NaF->SetLineWidth(2);
    
    TGraph *grLinearMean_AgL = new TGraph(ekValues.size(), ekValues.data(), linearMeans_AgL.data());
    grLinearMean_AgL->SetMarkerStyle(21);
    grLinearMean_AgL->SetMarkerSize(2.0);
    grLinearMean_AgL->SetMarkerColor(kRed);
    grLinearMean_AgL->SetLineColor(kRed);
    grLinearMean_AgL->SetLineWidth(2);
    
    TGraph *grNonlinearMean_AgL = new TGraph(ekValues.size(), ekValues.data(), nonlinearMeans_AgL.data());
    grNonlinearMean_AgL->SetMarkerStyle(21);
    grNonlinearMean_AgL->SetMarkerSize(2.0);
    grNonlinearMean_AgL->SetMarkerColor(kViolet);
    grNonlinearMean_AgL->SetLineColor(kViolet);
    grNonlinearMean_AgL->SetLineWidth(2);
    
    // Error graphs
    TGraph *grLinearError_NaF = new TGraph(ekValues.size(), ekValues.data(), linearErrors_NaF.data());
    grLinearError_NaF->SetMarkerStyle(24);
    grLinearError_NaF->SetMarkerSize(1.5);
    grLinearError_NaF->SetMarkerColor(kBlue);
    grLinearError_NaF->SetLineColor(kBlue);
    grLinearError_NaF->SetLineStyle(2);
    grLinearError_NaF->SetLineWidth(2);
    
    TGraph *grNonlinearError_NaF = new TGraph(ekValues.size(), ekValues.data(), nonlinearErrors_NaF.data());
    grNonlinearError_NaF->SetMarkerStyle(24);
    grNonlinearError_NaF->SetMarkerSize(1.5);
    grNonlinearError_NaF->SetMarkerColor(kGreen+2);
    grNonlinearError_NaF->SetLineColor(kGreen+2);
    grNonlinearError_NaF->SetLineStyle(2);
    grNonlinearError_NaF->SetLineWidth(2);
    
    TGraph *grLinearError_AgL = new TGraph(ekValues.size(), ekValues.data(), linearErrors_AgL.data());
    grLinearError_AgL->SetMarkerStyle(25);
    grLinearError_AgL->SetMarkerSize(1.5);
    grLinearError_AgL->SetMarkerColor(kRed);
    grLinearError_AgL->SetLineColor(kRed);
    grLinearError_AgL->SetLineStyle(2);
    grLinearError_AgL->SetLineWidth(2);
    
    TGraph *grNonlinearError_AgL = new TGraph(ekValues.size(), ekValues.data(), nonlinearErrors_AgL.data());
    grNonlinearError_AgL->SetMarkerStyle(25);
    grNonlinearError_AgL->SetMarkerSize(1.5);
    grNonlinearError_AgL->SetMarkerColor(kViolet);
    grNonlinearError_AgL->SetLineColor(kViolet);
    grNonlinearError_AgL->SetLineStyle(2);
    grNonlinearError_AgL->SetLineWidth(2);
    
    // Create zero residual reference line
    TF1 *perfectResidual = new TF1("perfectResidual", "0", minEk, maxEk);
    perfectResidual->SetLineColor(kBlack);
    perfectResidual->SetLineStyle(2);
    
    // Draw comparison plots
    // ------------------------------------------------------------------------

    // Find optimal y-axis range for means based on actual data
    double meanMin = 0, meanMax = 0;
    double errorMin = 0, errorMax = 0;
    
    // Helper function to find min/max values
    auto updateRange = [](const std::vector<double>& data, double& minVal, double& maxVal) {
        if (data.empty()) return;
        
        for (double val : data) {
            if (std::isnan(val) || std::isinf(val)) continue;
            
            if (val < minVal) minVal = val;
            if (val > maxVal) maxVal = val;
        }
    };
    
    // Find min/max for means
    meanMin = std::numeric_limits<double>::max();
    meanMax = std::numeric_limits<double>::lowest();
    
    updateRange(linearMeans_NaF, meanMin, meanMax);
    updateRange(nonlinearMeans_NaF, meanMin, meanMax);
    updateRange(linearMeans_AgL, meanMin, meanMax);
    updateRange(nonlinearMeans_AgL, meanMin, meanMax);
    
    // Add some padding (20%)
    double meanPadding = (meanMax - meanMin) * 0.2;
    meanMin = std::max(meanMin - meanPadding, -0.1);
    meanMax = std::min(meanMax + meanPadding, 0.1);
    
    // Ensure reasonable limits if data is very sparse
    if (meanMin > -0.005) meanMin = -0.005;
    if (meanMax < 0.005) meanMax = 0.005;
    
    // Find min/max for errors
    errorMin = 0; // Errors are always positive
    errorMax = std::numeric_limits<double>::lowest();
    
    updateRange(linearErrors_NaF, errorMin, errorMax);
    updateRange(nonlinearErrors_NaF, errorMin, errorMax);
    updateRange(linearErrors_AgL, errorMin, errorMax);
    updateRange(nonlinearErrors_AgL, errorMin, errorMax);
    
    // Add some padding (20%)
    double errorPadding = errorMax * 0.2;
    errorMax = std::min(errorMax + errorPadding, 0.1);
    
    // Ensure reasonable limits if data is very sparse
    if (errorMax < 0.01) errorMax = 0.01;
    
    // Set up dummy histograms with dynamic ranges
    TH2F *dummyHistMean = new TH2F("dummyHistMean", "Beta Residual Means vs. E_{k}/n;E_{k}/n [GeV/nucleon];1/#beta_{rec} - 1/#beta_{MC}",
                                  100, minEk, maxEk, 100, meanMin, meanMax);
    dummyHistMean->SetStats(0);
    
    TH2F *dummyHistError = new TH2F("dummyHistError", "Beta Residual Errors vs. E_{k}/n;E_{k}/n [GeV/nucleon];#sigma(1/#beta_{rec} - 1/#beta_{MC})",
                                   100, minEk, maxEk, 100, errorMin, errorMax);
    dummyHistError->SetStats(0);

    // Draw mean plots
    // ------------------------------------------------------------------------
    
    canvas->SetLogx(); // Set logarithmic x-axis

    // All means comparison
    dummyHistMean->Draw();
    perfectResidual->Draw("SAME");
    grLinearMean_NaF->Draw("LP");
    grNonlinearMean_NaF->Draw("LP");
    grLinearMean_AgL->Draw("LP");
    grNonlinearMean_AgL->Draw("LP");
    
    TLegend *legendMean = new TLegend(0.65, 0.15, 0.85, 0.35);
    legendMean->SetBorderSize(0);
    legendMean->AddEntry(grLinearMean_NaF, "NaF Linear", "lp");
    legendMean->AddEntry(grNonlinearMean_NaF, "NaF Non-linear", "lp");
    legendMean->AddEntry(grLinearMean_AgL, "AgL Linear", "lp");
    legendMean->AddEntry(grNonlinearMean_AgL, "AgL Non-linear", "lp");
    legendMean->Draw();
    
    infoText->Draw();
    canvas->Print(outputName);
    
    // Linear means only
    dummyHistMean->Draw();
    perfectResidual->Draw("SAME");
    grLinearMean_NaF->Draw("LP");
    grLinearMean_AgL->Draw("LP");
    
    TLegend *legendLinearMean = new TLegend(0.65, 0.15, 0.85, 0.25);
    legendLinearMean->SetBorderSize(0);
    legendLinearMean->AddEntry(grLinearMean_NaF, "NaF Linear", "lp");
    legendLinearMean->AddEntry(grLinearMean_AgL, "AgL Linear", "lp");
    legendLinearMean->Draw();
    
    infoText->Draw();
    canvas->Print(outputName);
    
    // Non-linear means only
    dummyHistMean->Draw();
    perfectResidual->Draw("SAME");
    grNonlinearMean_NaF->Draw("LP");
    grNonlinearMean_AgL->Draw("LP");
    
    TLegend *legendNonlinearMean = new TLegend(0.65, 0.15, 0.85, 0.25);
    legendNonlinearMean->SetBorderSize(0);
    legendNonlinearMean->AddEntry(grNonlinearMean_NaF, "NaF Non-linear", "lp");
    legendNonlinearMean->AddEntry(grNonlinearMean_AgL, "AgL Non-linear", "lp");
    legendNonlinearMean->Draw();
    
    infoText->Draw();
    canvas->Print(outputName);

    // Draw error plots
    // ------------------------------------------------------------------------
    
    // All errors comparison
    dummyHistError->Draw();
    grLinearError_NaF->Draw("LP");
    grNonlinearError_NaF->Draw("LP");
    grLinearError_AgL->Draw("LP");
    grNonlinearError_AgL->Draw("LP");
    
    TLegend *legendError = new TLegend(0.65, 0.15, 0.85, 0.35);
    legendError->SetBorderSize(0);
    legendError->AddEntry(grLinearError_NaF, "NaF Linear", "lp");
    legendError->AddEntry(grNonlinearError_NaF, "NaF Non-linear", "lp");
    legendError->AddEntry(grLinearError_AgL, "AgL Linear", "lp");
    legendError->AddEntry(grNonlinearError_AgL, "AgL Non-linear", "lp");
    legendError->Draw();
    
    infoText->Draw();
    canvas->Print(outputName);
    
    // Linear errors only
    dummyHistError->Draw();
    grLinearError_NaF->Draw("LP");
    grLinearError_AgL->Draw("LP");
    
    TLegend *legendLinearError = new TLegend(0.65, 0.15, 0.85, 0.25);
    legendLinearError->SetBorderSize(0);
    legendLinearError->AddEntry(grLinearError_NaF, "NaF Linear", "lp");
    legendLinearError->AddEntry(grLinearError_AgL, "AgL Linear", "lp");
    legendLinearError->Draw();
    
    infoText->Draw();
    canvas->Print(outputName);
    
    // Non-linear errors only
    dummyHistError->Draw();
    grNonlinearError_NaF->Draw("LP");
    grNonlinearError_AgL->Draw("LP");
    
    TLegend *legendNonlinearError = new TLegend(0.65, 0.15, 0.85, 0.25);
    legendNonlinearError->SetBorderSize(0);
    legendNonlinearError->AddEntry(grNonlinearError_NaF, "NaF Non-linear", "lp");
    legendNonlinearError->AddEntry(grNonlinearError_AgL, "AgL Non-linear", "lp");
    legendNonlinearError->Draw();
    
    infoText->Draw();
    canvas->Print(outputName);
    
    // Close PDF file
    canvas->Print(Form("%s]", outputName));
    
    std::cout << "E_k/n vs Beta residuals plot (NaF vs AgL comparison) saved to: " << outputName << std::endl;
    std::cout << "Dynamic Y-axis ranges: Means [" << meanMin << ", " << meanMax << "], Errors [" << errorMin << ", " << errorMax << "]" << std::endl;
    
    // Cleanup
    delete dummyHistMean;
    delete dummyHistError;
    file->Close();
    delete file;
}