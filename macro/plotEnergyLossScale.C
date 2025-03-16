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
#include <TPaveText.h>
#include <TGraphErrors.h>

/**
 * Draw energy loss scale distribution from a ROOT file containing scaleTree tree
 *
 * @param fileName Path to the input ROOT file
 * @param outputName Output file name
 */
void plotEnergyLossScale(std::string fileName = "test.root",
                         const char *outputName = "test_energy_loss_scale.pdf")
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

    // Get the energy loss scale tree
    TTree *tree = (TTree *)file->Get("scaleTree");
    if (!tree)
    {
        std::cerr << "Error: Could not find scaleTree in " << fileName << std::endl;
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

    // Helper function to get range using quantiles
    auto getQuantileRange = [&tree](const char *branchName, double lowQuantile = 0.01, double highQuantile = 0.99)
    {
        TH1D hTemp(Form("hTemp_%s", branchName), "", 10000, tree->GetMinimum(branchName), tree->GetMaximum(branchName));
        tree->Draw(Form("%s>>hTemp_%s", branchName, branchName), "", "goff");

        Double_t xq[2] = {lowQuantile, highQuantile};
        Double_t yq[2] = {0, 0};
        hTemp.GetQuantiles(2, yq, xq);

        std::cout << "Range for " << branchName << " using quantiles: [" << yq[0] << ", " << yq[1] << "]" << std::endl;
        return std::make_pair(yq[0], yq[1]);
    };

    // Get ranges for variables using quantiles
    auto energyLossScaleRange = getQuantileRange("energyLossScale");

    // Extract values for readability
    double energyLossScaleMin = energyLossScaleRange.first;
    double energyLossScaleMax = energyLossScaleRange.second;
    double mcBetaMin = tree->GetMinimum("mcBeta");
    double mcBetaMax = tree->GetMaximum("mcBeta");

    // Number of bins for histograms
    int nBins = 100;
    int nBinsX = 40;  // Number of bins in beta direction
    int nBinsY = 100; // Number of bins in scale direction

    // Create 1D histogram
    TH1F *hEnergyLossScale = new TH1F("hEnergyLossScale",
                                      ";Energy Loss Scale Factor;Entries",
                                      nBins, energyLossScaleMin, energyLossScaleMax);

    // Create 2D histogram
    TH2F *hScaleVsBeta = new TH2F("hScaleVsBeta",
                                  ";#beta_{MC};Energy Loss Scale Factor",
                                  nBinsX, mcBetaMin, mcBetaMax,
                                  nBinsY, energyLossScaleMin, energyLossScaleMax);

#if false
    // Define variables for reading from tree
    float energyLossScale;
    float mcBeta;
    float position[3];
    float direction[3];

    // Set branch addresses
    tree->SetBranchAddress("energyLossScale", &energyLossScale);
    tree->SetBranchAddress("mcBeta", &mcBeta);
    tree->SetBranchAddress("position", position);
    tree->SetBranchAddress("direction", direction);
#endif

    // Create canvas and draw 1D histogram
    TCanvas *canvas1 = new TCanvas("canvas1", "Energy Loss Scale Distribution", 800, 600);
    canvas1->SetLeftMargin(0.16);
    canvas1->SetRightMargin(0.11);
    canvas1->SetGridx();
    canvas1->SetGridy();

    tree->Draw("energyLossScale>>hEnergyLossScale");
    hEnergyLossScale->SetLineColor(kBlue);
    hEnergyLossScale->SetLineWidth(2);
    hEnergyLossScale->SetFillColor(kBlue - 10);
    hEnergyLossScale->SetFillStyle(3004);
    hEnergyLossScale->Draw();

    canvas1->Print(outputName);

    // Create canvas for 2D plot
    TCanvas *canvas2 = new TCanvas("canvas2", "Energy Loss Scale vs Beta", 800, 600);
    canvas2->SetLeftMargin(0.16);
    canvas2->SetRightMargin(0.11);
    canvas2->SetGridx();
    canvas2->SetGridy();
    canvas2->SetLogz();

    tree->Draw("energyLossScale:mcBeta>>hScaleVsBeta");
    hScaleVsBeta->Draw("COLZ");

    // Prepare for profile analysis
    const int nProfiles = nBinsX;
    double betaValues[nProfiles];
    double scaleMean[nProfiles], scaleError[nProfiles];
    double binWidth = (mcBetaMax - mcBetaMin) / nBinsX;

    // Print table header before the loop
    std::cout << "\nBin\tBeta\tEntries\tMean\n"
              << std::string(65, '-') << std::endl;

    // Process each bin
    for (int bin = 0; bin < nBinsX; ++bin)
    {
        double binCenter = mcBetaMin + (bin + 0.5) * binWidth;
        int binIdx = bin + 1; // ROOT histograms are 1-indexed

        // Calculate bin range for projection
        int binWindow = 1; // Use ±1 bins for better statistics
        int binLow = TMath::Max(1, binIdx - binWindow);
        int binHigh = TMath::Min(nBinsX, binIdx + binWindow);

        betaValues[bin] = binCenter;

        // Get projection and fit
        TH1D *proj = hScaleVsBeta->ProjectionY(Form("proj_%d", bin), binLow, binHigh, "e");
        if (proj->GetEntries() > 0)
        {
            proj->Fit("gaus", "QN", "", energyLossScaleMin, energyLossScaleMax);
            TF1 *fit = proj->GetFunction("gaus");
            if (fit && fit->GetProb() > 0.01)
            {
                scaleMean[bin] = fit->GetParameter(1);
                scaleError[bin] = fit->GetParameter(2);
            }
            else
            {
                scaleMean[bin] = proj->GetMean();
                scaleError[bin] = proj->GetRMS();
            }
        }
        else
        {
            scaleMean[bin] = 0;
            scaleError[bin] = 0;
        }

        printf("%d\t%.5f\t%d\t%.5f\n",
               bin, binCenter, (Int_t)proj->GetEntries(), scaleMean[bin]);

        delete proj;
    }

    // Create and draw TGraphErrors
    TGraphErrors *grScale = new TGraphErrors(nProfiles, betaValues, scaleMean, 0, scaleError);
    grScale->SetMarkerStyle(20);
    grScale->SetMarkerColor(kBlack);
    grScale->SetMarkerSize(1.5);
    grScale->Draw("P SAME");

    // Update canvas
    canvas2->Modified();
    canvas2->Update();
    canvas2->Print(Form("%s)", outputName)); // Append to PDF

    // Close file
    file->Close();
    std::cout << "Energy loss scale plots saved to: " << outputName << std::endl;
}