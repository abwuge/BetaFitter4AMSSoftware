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
#include <TPaveText.h>

/**
 * Draw energy loss distributions from a ROOT file containing energyLoss tree
 *
 * @param fileName Path to the input ROOT file
 * @param outputName Output file name (default: energy_loss_plots.pdf)
 */
void plotEnergyLoss(const char *fileName = "test.root",
                    const char *outputName = "test_energy_loss_plots.pdf")
{
    // Set batch mode to avoid GUI related issues
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);

    // Open the ROOT file
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return;
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

    // Set branch addresses
    tree->SetBranchAddress("energyDepositedS1S2", &energyDepositedS1S2);
    tree->SetBranchAddress("energyDepositedTotal", &energyDepositedTotal);
    tree->SetBranchAddress("energyLoss_S1S2_", &energyLoss_S1S2_);
    tree->SetBranchAddress("energyLoss_S1S4_", &energyLoss_S1S4_);
    tree->SetBranchAddress("energyLossScaleS1S2", &energyLossScaleS1S2);
    tree->SetBranchAddress("energyLossScaleTotal", &energyLossScaleTotal);
    tree->SetBranchAddress("energyLossS2__S3", &energyLossS2__S3);
    tree->SetBranchAddress("energyLossS2S3_Total", &energyLossS2S3_Total);

    // Create histograms
    int nBins = 100; // Number of bins for histograms

    TH1F *hEnergyLossScaleS1S2 = new TH1F("hEnergyLossScaleS1S2",
                                          ";Energy Loss Scale Factor (S1-S2);",
                                          nBins, 0, 3);

    TH1F *hEnergyLossScaleTotal = new TH1F("hEnergyLossScaleTotal",
                                           ";Total Energy Loss Scale Factor (S1-S4);",
                                           nBins, 0, 4);

    TH1F *hEnergyLossS2__S3 = new TH1F("hEnergyLossS2__S3",
                                       ";Energy Loss Between S2 and S3 [GeV];",
                                       nBins, 0, 1);

    TH1F *hEnergyLossS2S3_Total = new TH1F("hEnergyLossS2S3_Total",
                                      ";Energy Loss Between S2 and S3 (Total);",
                                      nBins, 0, 1);
    
    TH2F *hEnergyLossS1S2 = new TH2F("hEnergyLossS1S2",
                                     ";Energy Deposited (S1-S2) [GeV];Energy Loss (S1-S2) [GeV];",
                                     nBins, 0.1, 0.7, nBins, 0, 1);

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
    }
    
    // Create canvas with multiple pages
    TCanvas *canvas = new TCanvas("canvas", "Energy Loss Distributions", 3508, 2480);
    canvas->SetLeftMargin(0.16);
    canvas->SetRightMargin(0.11);
    canvas->SetGridx();
    canvas->SetGridy();
        
    hEnergyLossScaleS1S2->Draw();
    canvas->Print(Form("%s(", outputName));

    hEnergyLossScaleTotal->Draw();
    canvas->Print(outputName);
    
    hEnergyLossS2__S3->Draw();
    canvas->Print(outputName);
    
    hEnergyLossS2S3_Total->Draw();
    canvas->Print(outputName);
    
    hEnergyLossS1S2->Draw("colz");
    canvas->Print(Form("%s)", outputName));
    
    // Clean up
    delete hEnergyLossS2__S3;
    delete hEnergyLossScaleTotal;
    delete hEnergyLossScaleS1S2;
    delete hEnergyLossS2S3_Total;
    delete hEnergyLossS1S2;
    delete canvas;
    
    // Reset tree branches to avoid dangling pointers
    tree->ResetBranchAddresses();
    
    // Close and delete the file after all references are gone
    file->Close();
    delete file;
    
    std::cout << "Energy loss plots saved to: " << outputName << std::endl;
}