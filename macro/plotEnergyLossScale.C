#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TF1.h>
#include <TLegend.h>
#include <TROOT.h>
#include <iostream>
#include <string>
#include <TPaveText.h>

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

    // Get ranges for variables using quantiles
    auto energyLossScaleRange = getQuantileRange("energyLossScale");
    
    // Extract values for readability
    double energyLossScaleMin = energyLossScaleRange.first;
    double energyLossScaleMax = energyLossScaleRange.second;
    
    // Add some margin to ranges (5%)
    auto addMargin = [](double min, double max) {
        double margin = 0.05 * (max - min);
        return std::make_pair(min - margin, max + margin);
    };
    
    // Apply margins to ranges
    std::tie(energyLossScaleMin, energyLossScaleMax) = addMargin(energyLossScaleMin, energyLossScaleMax);

    // Number of bins for histograms
    int nBins = 100;

    // Create 1D histogram with range from quantiles
    TH1F *hEnergyLossScale = new TH1F("hEnergyLossScale",
                                    ";Energy Loss Scale Factor;Entries",
                                    nBins, energyLossScaleMin, energyLossScaleMax);

    // Fill histograms
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);
        hEnergyLossScale->Fill(energyLossScale);
    }

    // Create canvas
    TCanvas *canvas = new TCanvas("canvas", "Energy Loss Scale Distribution", 800, 600);
    canvas->SetLeftMargin(0.16);
    canvas->SetRightMargin(0.11);
    canvas->SetGridx();
    canvas->SetGridy();
    
    // Draw histogram
    hEnergyLossScale->SetLineColor(kBlue);
    hEnergyLossScale->SetLineWidth(2);
    hEnergyLossScale->SetFillColor(kBlue-10);
    hEnergyLossScale->SetFillStyle(3004);
    hEnergyLossScale->Draw();
    
    // Print to PDF
    canvas->Print(outputName);
    
    // Close file
    file->Close();
    std::cout << "Energy loss scale plot saved to: " << outputName << std::endl;
}