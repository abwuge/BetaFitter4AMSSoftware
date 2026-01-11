#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLegend.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <TMath.h>
#include <TGraphErrors.h>

/**
 * Draw energy loss distributions from a ROOT file containing energyLoss tree
 *
 * @param fileName Path to the input ROOT file
 * @param outputName Output file name
 */
void plotEnergyLoss(std::string fileName = "test.root",
                    const char *outputName = "test_energy_loss_plots.pdf")
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
        std::string resultsPath = "results/" + fileName;
        std::cout << "Try to open file: " << resultsPath << std::endl;

        file = TFile::Open(resultsPath.c_str(), "READ");
        if (!file || file->IsZombie())
        {
            std::cerr << "Error: Unable to open file " << fileName << " or " << resultsPath << std::endl;
            return;
        }
    }

    // Get the energy loss tree
    TTree *tree = (TTree *)file->Get("energyLoss");
    if (!tree)
    {
        std::cerr << "Error: Could not find energyLoss tree in " << fileName << std::endl;
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

    // Define variables for reading from tree
    float energyDepositedS1S2 = 0.0;  // energy deposited from before S1 to after S2
    float energyDepositedTotal = 0.0; // total energy deposited
    float energyLoss_S1S2_ = 0.0;     // energy loss from before S1 to after S2
    float energyLoss_S1S4_ = 0.0;     // energy loss from before S1 to after S4
    float energyLossScaleS1S2 = 0.0;  // energy loss scale factor from before S1 to after S2
    float energyLossScaleTotal = 0.0; // energy loss scale factor from before S1 to after S4
    float energyLossS2__S3 = 0.0;     // energy loss from after S2 to before S3
    float energyLossS2S3_Total = 0.0; // energy loss from after S2 to before S3, normalized to total energy loss
    float mcBeta = 0.0;               // Monte Carlo beta
    int tof_qs = 0;                   // Q Status (1111: all unoverlapped, 0000: all overlapped)

    // Set branch addresses
    tree->SetBranchAddress("energyDepositedS1S2", &energyDepositedS1S2);
    tree->SetBranchAddress("energyDepositedTotal", &energyDepositedTotal);
    tree->SetBranchAddress("energyLoss_S1S2_", &energyLoss_S1S2_);
    tree->SetBranchAddress("energyLoss_S1S4_", &energyLoss_S1S4_);
    tree->SetBranchAddress("energyLossScaleS1S2", &energyLossScaleS1S2);
    tree->SetBranchAddress("energyLossScaleTotal", &energyLossScaleTotal);
    tree->SetBranchAddress("energyLossS2__S3", &energyLossS2__S3);
    tree->SetBranchAddress("energyLossS2S3_Total", &energyLossS2S3_Total);
    tree->SetBranchAddress("mcBeta", &mcBeta);
    tree->SetBranchAddress("tof_qs", &tof_qs);

    // Helper function to get range using quantiles
    auto getQuantileRange = [&tree](const char* branchName, double lowQuantile = 0.01, double highQuantile = 0.99) {
        // Create a high-precision temporary histogram
        TH1D hTemp(Form("hTemp_%s", branchName), "", 10000, tree->GetMinimum(branchName), tree->GetMaximum(branchName));
        tree->Draw(Form("%s>>hTemp_%s", branchName, branchName), "", "goff");
        
        // Calculate quantiles
        Double_t xq[2] = {lowQuantile, highQuantile};
        Double_t yq[2] = {0, 0};
        hTemp.GetQuantiles(2, yq, xq);
        
        std::cout << "Range for " << branchName << " using quantiles: [" << yq[0] << ", " << yq[1] << "]" << std::endl;
        return std::make_pair(yq[0], yq[1]);
    };

    // Get ranges for all variables using quantiles
    auto energyDepositedS1S2Range = getQuantileRange("energyDepositedS1S2");
    auto energyLoss_S1S2_Range = getQuantileRange("energyLoss_S1S2_");
    auto energyLossScaleS1S2Range = getQuantileRange("energyLossScaleS1S2");
    auto energyLossScaleTotalRange = getQuantileRange("energyLossScaleTotal");
    auto energyLossS2__S3Range = getQuantileRange("energyLossS2__S3");
    auto energyLossS2S3_TotalRange = getQuantileRange("energyLossS2S3_Total");

    // Extract values for readability
    double mcBetaMin = tree->GetMinimum("mcBeta");
    double mcBetaMax = tree->GetMaximum("mcBeta");
    double energyDepositedS1S2Min = energyDepositedS1S2Range.first;
    double energyDepositedS1S2Max = energyDepositedS1S2Range.second;
    double energyLoss_S1S2_Min = energyLoss_S1S2_Range.first;
    double energyLoss_S1S2_Max = energyLoss_S1S2_Range.second;
    double energyLossScaleS1S2Min = energyLossScaleS1S2Range.first;
    double energyLossScaleS1S2Max = energyLossScaleS1S2Range.second;
    double energyLossScaleTotalMin = energyLossScaleTotalRange.first;
    double energyLossScaleTotalMax = energyLossScaleTotalRange.second;
    double energyLossS2__S3Min = energyLossS2__S3Range.first;
    double energyLossS2__S3Max = energyLossS2__S3Range.second;
    double energyLossS2S3_TotalMin = energyLossS2S3_TotalRange.first;
    double energyLossS2S3_TotalMax = energyLossS2S3_TotalRange.second;

    // Number of bins for histograms
    int nBins = 100;  // Number of bins for 1D histograms
    int nBinsX = 40;  // Number of bins in X direction (beta) for 2D histograms
    int nBinsY = 100; // Number of bins in Y direction for 2D histograms

    // Create 1D histograms with ranges from quantiles
    TH1F *hEnergyLossScaleS1S2 = new TH1F("hEnergyLossScaleS1S2",
                                          ";Energy Loss Scale Factor (S1-S2);",
                                          nBins, energyLossScaleS1S2Min, energyLossScaleS1S2Max);

    TH1F *hEnergyLossScaleTotal = new TH1F("hEnergyLossScaleTotal",
                                           ";Total Energy Loss Scale Factor (S1-S4);",
                                           nBins, energyLossScaleTotalMin, energyLossScaleTotalMax);

    TH1F *hEnergyLossS2__S3 = new TH1F("hEnergyLossS2__S3",
                                       ";Energy Loss Between S2 and S3 [GeV];",
                                       nBins, energyLossS2__S3Min, energyLossS2__S3Max);

    TH1F *hEnergyLossS2S3_Total = new TH1F("hEnergyLossS2S3_Total",
                                           ";Energy Loss Between S2 and S3 (Total);",
                                           nBins, energyLossS2S3_TotalMin, energyLossS2S3_TotalMax);

    TH2F *hEnergyLossScaleS1S2VsBeta = new TH2F("hEnergyLossScaleS1S2VsBeta",
                                                ";#beta_{MC};Energy Loss Scale Factor (S1-S2)",
                                                nBinsX, mcBetaMin, mcBetaMax, 
                                                nBinsY, energyLossScaleS1S2Min, energyLossScaleS1S2Max);

    TH2F *hEnergyLossScaleTotalVsBeta = new TH2F("hEnergyLossScaleTotalVsBeta",
                                                 ";#beta_{MC};Total Energy Loss Scale Factor (S1-S4)",
                                                 nBinsX, mcBetaMin, mcBetaMax, 
                                                 nBinsY, energyLossScaleTotalMin, energyLossScaleTotalMax);

    TH2F *hEnergyLossS2__S3VsBeta = new TH2F("hEnergyLossS2__S3VsBeta",
                                             ";#beta_{MC};Energy Loss Between S2 and S3 [GeV]",
                                             nBinsX, mcBetaMin, mcBetaMax, 
                                             nBinsY, energyLossS2__S3Min, energyLossS2__S3Max);

    TH2F *hEnergyLossS2S3_TotalVsBeta = new TH2F("hEnergyLossS2S3_TotalVsBeta",
                                                 ";#beta_{MC};Energy Loss Between S2 and S3 (Total)",
                                                 nBinsX, mcBetaMin, mcBetaMax, 
                                                 nBinsY, energyLossS2S3_TotalMin, energyLossS2S3_TotalMax);

    TH2F *hEnergyLossS1S2 = new TH2F("hEnergyLossS1S2",
                                     ";Energy Deposited (S1-S2) [GeV];Energy Loss (S1-S2) [GeV];",
                                     nBins, energyDepositedS1S2Min, energyDepositedS1S2Max, 
                                     nBins, energyLoss_S1S2_Min, energyLoss_S1S2_Max);

    TH2F *hEnergyLossS1S2_Unoverlapped = new TH2F("hEnergyLossS1S2_Unoverlapped",
                                                  "Unoverlapped S1&S2;Energy Deposited (S1-S2) [GeV];Energy Loss (S1-S2) [GeV];",
                                                  nBins, energyDepositedS1S2Min, energyDepositedS1S2Max, 
                                                  nBins, energyLoss_S1S2_Min, energyLoss_S1S2_Max);

    TH2F *hEnergyLossS1S2_Overlapped = new TH2F("hEnergyLossS1S2_Overlapped",
                                                "Overlapped S1&S2;Energy Deposited (S1-S2) [GeV];Energy Loss (S1-S2) [GeV];",
                                                nBins, energyDepositedS1S2Min, energyDepositedS1S2Max, 
                                                nBins, energyLoss_S1S2_Min, energyLoss_S1S2_Max);

    // Fill histograms
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);

        hEnergyLossScaleS1S2->Fill(energyLossScaleS1S2);
        hEnergyLossScaleTotal->Fill(energyLossScaleTotal);
        hEnergyLossS2__S3->Fill(energyLossS2__S3);
        hEnergyLossS2S3_Total->Fill(energyLossS2S3_Total);
        hEnergyLossS1S2->Fill(energyDepositedS1S2, energyLoss_S1S2_);

        // Check if S1 and S2 are unoverlapped (11xx)
        bool isS1S2Unoverlapped = tof_qs / 100 == 11;
        if (isS1S2Unoverlapped)
        {
            hEnergyLossS1S2_Unoverlapped->Fill(energyDepositedS1S2, energyLoss_S1S2_);
        }
        else
        {
            hEnergyLossS1S2_Overlapped->Fill(energyDepositedS1S2, energyLoss_S1S2_);
        }

        // Fill 2D histograms vs Beta
        hEnergyLossScaleS1S2VsBeta->Fill(mcBeta, energyLossScaleS1S2);
        hEnergyLossScaleTotalVsBeta->Fill(mcBeta, energyLossScaleTotal);
        hEnergyLossS2__S3VsBeta->Fill(mcBeta, energyLossS2__S3);
        hEnergyLossS2S3_TotalVsBeta->Fill(mcBeta, energyLossS2S3_Total);
    }

    // Arrays to store fit results for each 2D histogram
    const int nProfiles = nBinsX;
    double betaValues[nProfiles];
    double scaleS1S2Mean[nProfiles], scaleS1S2Error[nProfiles];
    double scaleTotalMean[nProfiles], scaleTotalError[nProfiles];
    double lossS2S3Mean[nProfiles], lossS2S3Error[nProfiles];
    double lossS2S3TotalMean[nProfiles], lossS2S3TotalError[nProfiles];

    // Print table header for the fit results
    std::cout << std::endl;
    std::cout << "Bin\tBeta\t\tS1S2_mean\tTotal_mean\tS2S3_mean\tS2S3Total_mean" << std::endl;
    std::cout << std::string(85, '-') << std::endl;

    // Process each bin to get fit results for all 2D histograms
    double binWidth = (mcBetaMax - mcBetaMin) / nBinsX;
    for (int bin = 0; bin < nBinsX; ++bin)
    {
        double binCenter = mcBetaMin + (bin + 0.5) * binWidth;
        int binIdx = bin + 1; // ROOT histograms are 1-indexed

        // Calculate bin range for projection with wider window
        int binWindow = 1; // Use ±1 bins for better statistics
        int binLow = TMath::Max(1, binIdx - binWindow);
        int binHigh = TMath::Min(nBinsX, binIdx + binWindow);

        betaValues[bin] = binCenter;

        // Process S1S2 scale distribution
        TH1D *projS1S2 = hEnergyLossScaleS1S2VsBeta->ProjectionY(Form("projS1S2_%d", bin), binLow, binHigh, "e");
        Int_t projS1S2Entries = projS1S2->GetEntries();

        if (projS1S2Entries)
        {
            projS1S2->Fit("gaus", "QN", "", 0, 3);
            TF1 *fitS1S2 = projS1S2->GetFunction("gaus");
            if (fitS1S2 && fitS1S2->GetProb() > 0.01)
            {
                scaleS1S2Mean[bin] = fitS1S2->GetParameter(1);
                scaleS1S2Error[bin] = fitS1S2->GetParameter(2);
            }
            else
            {
                scaleS1S2Mean[bin] = projS1S2->GetMean();
                scaleS1S2Error[bin] = projS1S2->GetRMS();
            }
        }
        else
        {
            scaleS1S2Mean[bin] = 0;
            scaleS1S2Error[bin] = 0;
        }
        delete projS1S2;

        // Process Total scale distribution
        TH1D *projTotal = hEnergyLossScaleTotalVsBeta->ProjectionY(Form("projTotal_%d", bin), binLow, binHigh, "e");
        Int_t projTotalEntries = projTotal->GetEntries();

        if (projTotalEntries)
        {
            projTotal->Fit("gaus", "QN", "", 0, 4);
            TF1 *fitTotal = projTotal->GetFunction("gaus");
            if (fitTotal && fitTotal->GetProb() > 0.01)
            {
                scaleTotalMean[bin] = fitTotal->GetParameter(1);
                scaleTotalError[bin] = fitTotal->GetParameter(2);
            }
            else
            {
                scaleTotalMean[bin] = projTotal->GetMean();
                scaleTotalError[bin] = projTotal->GetRMS();
            }
        }
        else
        {
            scaleTotalMean[bin] = 0;
            scaleTotalError[bin] = 0;
        }
        delete projTotal;

        // Process S2S3 energy loss distribution
        TH1D *projS2S3 = hEnergyLossS2__S3VsBeta->ProjectionY(Form("projS2S3_%d", bin), binLow, binHigh, "e");
        Int_t projS2S3Entries = projS2S3->GetEntries();

        if (projS2S3Entries)
        {
            projS2S3->Fit("gaus", "QN", "", 0, 1);
            TF1 *fitS2S3 = projS2S3->GetFunction("gaus");
            if (fitS2S3 && fitS2S3->GetProb() > 0.01)
            {
                lossS2S3Mean[bin] = fitS2S3->GetParameter(1);
                lossS2S3Error[bin] = fitS2S3->GetParameter(2);
            }
            else
            {
                lossS2S3Mean[bin] = projS2S3->GetMean();
                lossS2S3Error[bin] = projS2S3->GetRMS();
            }
        }
        else
        {
            lossS2S3Mean[bin] = 0;
            lossS2S3Error[bin] = 0;
        }
        delete projS2S3;

        // Process S2S3 Total distribution
        TH1D *projS2S3Tot = hEnergyLossS2S3_TotalVsBeta->ProjectionY(Form("projS2S3Tot_%d", bin), binLow, binHigh, "e");
        Int_t projS2S3TotEntries = projS2S3Tot->GetEntries();

        if (projS2S3TotEntries)
        {
            projS2S3Tot->Fit("gaus", "QN", "", 0, 1);
            TF1 *fitS2S3Tot = projS2S3Tot->GetFunction("gaus");
            if (fitS2S3Tot && fitS2S3Tot->GetProb() > 0.01)
            {
                lossS2S3TotalMean[bin] = fitS2S3Tot->GetParameter(1);
                lossS2S3TotalError[bin] = fitS2S3Tot->GetParameter(2);
            }
            else
            {
                lossS2S3TotalMean[bin] = projS2S3Tot->GetMean();
                lossS2S3TotalError[bin] = projS2S3Tot->GetRMS();
            }
        }
        else
        {
            lossS2S3TotalMean[bin] = 0;
            lossS2S3TotalError[bin] = 0;
        }
        delete projS2S3Tot;

        // Print fit results for this bin
        printf("%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
               bin, binCenter, scaleS1S2Mean[bin], scaleTotalMean[bin],
               lossS2S3Mean[bin], lossS2S3TotalMean[bin]);
    }

    // Create canvas with multiple pages
    TCanvas *canvas = new TCanvas("canvas", "Energy Loss Distributions", 3508, 2480);
    canvas->Print(Form("%s[", outputName)); // Open PDF file

    // Original 1D histograms
    canvas->SetLeftMargin(0.16);
    canvas->SetRightMargin(0.11);
    canvas->SetGridx();
    canvas->SetGridy();
    canvas->SetLogy();

    hEnergyLossScaleS1S2->Draw();
    canvas->Print(outputName);

    hEnergyLossScaleTotal->Draw();
    canvas->Print(outputName);

    hEnergyLossS2__S3->Draw();
    canvas->Print(outputName);

    hEnergyLossS2S3_Total->Draw();
    canvas->Print(outputName);

    hEnergyLossS1S2->Draw("colz");
    canvas->Print(outputName);

    hEnergyLossS1S2_Unoverlapped->Draw("colz");
    canvas->Print(outputName);

    hEnergyLossS1S2_Overlapped->Draw("colz");
    canvas->Print(outputName);

    // 2D histograms with beta correlation
    TCanvas *c1 = new TCanvas("c1", "Energy Loss Scale S1S2 vs Beta", 3508, 2480);
    c1->SetLeftMargin(0.16);
    c1->SetRightMargin(0.11);
    c1->SetGridx();
    c1->SetGridy();
    c1->SetLogz();

    hEnergyLossScaleS1S2VsBeta->Draw("COLZ");

    // Create and draw TGraphErrors for S1S2 scale
    TGraphErrors *grS1S2 = new TGraphErrors(nProfiles, betaValues, scaleS1S2Mean, 0, scaleS1S2Error);
    grS1S2->SetMarkerStyle(20);
    grS1S2->SetMarkerColor(kBlack);
    grS1S2->SetMarkerSize(3.0);
    grS1S2->Draw("P");

    c1->Print(outputName);

    // Energy Loss Scale Total vs Beta
    TCanvas *c2 = new TCanvas("c2", "Energy Loss Scale Total vs Beta", 3508, 2480);
    c2->SetLeftMargin(0.16);
    c2->SetRightMargin(0.11);
    c2->SetGridx();
    c2->SetGridy();
    c2->SetLogz();

    hEnergyLossScaleTotalVsBeta->Draw("COLZ");

    // Create and draw TGraphErrors for Total scale
    TGraphErrors *grTotal = new TGraphErrors(nProfiles, betaValues, scaleTotalMean, 0, scaleTotalError);
    grTotal->SetMarkerStyle(20);
    grTotal->SetMarkerColor(kBlack);
    grTotal->SetMarkerSize(3.0);
    grTotal->Draw("P");

    c2->Print(outputName);

    // Energy Loss S2S3 vs Beta
    TCanvas *c3 = new TCanvas("c3", "Energy Loss S2S3 vs Beta", 3508, 2480);
    c3->SetLeftMargin(0.16);
    c3->SetRightMargin(0.11);
    c3->SetGridx();
    c3->SetGridy();
    c3->SetLogz();

    hEnergyLossS2__S3VsBeta->Draw("COLZ");

    // Create and draw TGraphErrors for S2S3 loss
    TGraphErrors *grS2S3 = new TGraphErrors(nProfiles, betaValues, lossS2S3Mean, 0, lossS2S3Error);
    grS2S3->SetMarkerStyle(20);
    grS2S3->SetMarkerColor(kBlack);
    grS2S3->SetMarkerSize(3.0);
    grS2S3->Draw("P");

    c3->Print(outputName);

    // Energy Loss S2S3 Total vs Beta
    TCanvas *c4 = new TCanvas("c4", "Energy Loss S2S3 Total vs Beta", 3508, 2480);
    c4->SetLeftMargin(0.16);
    c4->SetRightMargin(0.11);
    c4->SetGridx();
    c4->SetGridy();
    c4->SetLogz();

    hEnergyLossS2S3_TotalVsBeta->Draw("COLZ");

    // Create and draw TGraphErrors for S2S3 total loss
    TGraphErrors *grS2S3Tot = new TGraphErrors(nProfiles, betaValues, lossS2S3TotalMean, 0, lossS2S3TotalError);
    grS2S3Tot->SetMarkerStyle(20);
    grS2S3Tot->SetMarkerColor(kBlack);
    grS2S3Tot->SetMarkerSize(3.0);
    grS2S3Tot->Draw("P");

    c4->Print(outputName);

    // Compare all energy loss trends in one plot
    TCanvas *c5 = new TCanvas("c5", "Energy Loss Trends vs Beta", 3508, 2480);
    c5->SetLeftMargin(0.16);
    c5->SetGridx();
    c5->SetGridy();

    // Create a frame for the comparison plot
    TH2F *hFrame = new TH2F("hFrame", "Energy Loss Trends vs Beta;#beta_{MC};Energy Loss Scale/Value",
                            100, mcBetaMin, mcBetaMax, 100, 0, 4);
    hFrame->Draw();

    // Draw all graphs with different styles
    grS1S2->SetMarkerStyle(20);
    grS1S2->SetMarkerColor(kBlue);
    grS1S2->SetLineColor(kBlue);
    grS1S2->SetMarkerSize(2.0);
    grS1S2->Draw("LP");

    grTotal->SetMarkerStyle(21);
    grTotal->SetMarkerColor(kRed);
    grTotal->SetLineColor(kRed);
    grTotal->SetMarkerSize(2.0);
    grTotal->Draw("LP");

    grS2S3->SetMarkerStyle(22);
    grS2S3->SetMarkerColor(kGreen + 2);
    grS2S3->SetLineColor(kGreen + 2);
    grS2S3->SetMarkerSize(2.0);
    grS2S3->Draw("LP");

    grS2S3Tot->SetMarkerStyle(23);
    grS2S3Tot->SetMarkerColor(kMagenta);
    grS2S3Tot->SetLineColor(kMagenta);
    grS2S3Tot->SetMarkerSize(2.0);
    grS2S3Tot->Draw("LP");

    // Add legend
    TLegend *legend = new TLegend(0.55, 0.7, 0.85, 0.85);
    legend->AddEntry(grS1S2, "Scale S1-S2", "lp");
    legend->AddEntry(grTotal, "Scale Total (S1-S4)", "lp");
    legend->AddEntry(grS2S3, "Loss S2-S3 [GeV]", "lp");
    legend->AddEntry(grS2S3Tot, "Loss S2-S3 (Total)", "lp");
    legend->SetBorderSize(0);
    legend->Draw();

    c5->Print(outputName);
    c5->Print(Form("%s]", outputName)); // Close PDF file

    // Close and delete the file after all references are gone
    file->Close();

    std::cout << "Energy loss plots saved to: " << outputName << std::endl;

    // Since this is the end of this macro
    // We do NOT need to delete the pointers
}