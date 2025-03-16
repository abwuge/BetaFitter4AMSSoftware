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
#include <TMath.h>

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

    // Extract values for readability
    double energyLossScaleMin = -2;
    double energyLossScaleMax = 6;
    double mcBetaMin = tree->GetMinimum("mcBeta");
    double mcBetaMax = tree->GetMaximum("mcBeta");

    // Number of bins for histograms
    int nBins = 200;
    int nBinsX = 40;  // Number of bins in beta direction
    int nBinsY = 100; // Number of bins in scale direction

    // Create 1D histogram
    TH1F *hEnergyLossScale = new TH1F("hEnergyLossScale",
                                      ";Energy Loss Scale Factor;Entries",
                                      nBins, energyLossScaleMin, energyLossScaleMax);

    // Create 2D histograms for angles
    TH2F *hScaleVsXZAngle = new TH2F("hScaleVsXZAngle",
                                     ";XZ Angle;Energy Loss Scale Factor",
                                     nBinsX, -TMath::Pi() * 0.3, TMath::Pi() * 0.3,
                                     nBinsY, energyLossScaleMin, energyLossScaleMax);

    TH2F *hScaleVsYZAngle = new TH2F("hScaleVsYZAngle",
                                     ";YZ Angle;Energy Loss Scale Factor",
                                     nBinsX, -TMath::Pi() * 0.3, TMath::Pi() * 0.3,
                                     nBinsY, energyLossScaleMin, energyLossScaleMax);

    // Create 2D histogram for beta
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

    tree->Draw("energyLossScale>>hEnergyLossScale", "mcBeta < 0.9");
    hEnergyLossScale->SetLineColor(kBlue);
    hEnergyLossScale->SetLineWidth(2);
    hEnergyLossScale->SetFillColor(kBlue - 10);
    hEnergyLossScale->SetFillStyle(3004);
    hEnergyLossScale->Draw();

    TPaveText *infoText = new TPaveText(0.2, 0.92, 0.8, 0.98, "NDC");
    infoText->SetFillColor(0);
    infoText->SetBorderSize(0);
    infoText->AddText(Form("mcBeta < 0.9"));
    infoText->Draw();

    canvas1->Print(Form("%s(", outputName));

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
    canvas2->Print(Form("%s", outputName));

    // Create canvas for XZ angle
    TCanvas *canvas3 = new TCanvas("canvas3", "Energy Loss Scale vs XZ Angle", 800, 600);
    canvas3->SetLeftMargin(0.16);
    canvas3->SetRightMargin(0.11);
    canvas3->SetGridx();
    canvas3->SetGridy();
    canvas3->SetLogz();

    tree->Draw("energyLossScale:TMath::ATan2(direction[0], direction[2])>>hScaleVsXZAngle", "mcBeta < 0.9");
    hScaleVsXZAngle->Draw("COLZ");

    // Prepare for profile analysis - XZ Angle
    const int nProfilesXZ = nBinsX;
    double xzAngleValues[nProfilesXZ];
    double xzScaleMean[nProfilesXZ], xzScaleError[nProfilesXZ];
    double xzBinWidth = 360.0 / nBinsX;

    // Print table header for XZ angle
    std::cout << "\nXZ Angle Bin\tAngle\tEntries\tMean\n"
              << std::string(65, '-') << std::endl;

    // Process each bin for XZ angle
    for (int bin = 0; bin < nBinsX; ++bin)
    {
        double binCenter = -180 + (bin + 0.5) * xzBinWidth;
        int binIdx = bin + 1; // ROOT histograms are 1-indexed

        // Calculate bin range for projection
        int binWindow = 1; // Use ±1 bins for better statistics
        int binLow = TMath::Max(1, binIdx - binWindow);
        int binHigh = TMath::Min(nBinsX, binIdx + binWindow);

        xzAngleValues[bin] = binCenter;

        // Get projection and fit
        TH1D *proj = hScaleVsXZAngle->ProjectionY(Form("projXZ_%d", bin), binLow, binHigh, "e");
        if (proj->GetEntries() > 10)
        {
            proj->Fit("gaus", "QN", "", energyLossScaleMin, energyLossScaleMax);
            TF1 *fit = proj->GetFunction("gaus");
            if (fit && fit->GetProb() > 0.01)
            {
                xzScaleMean[bin] = fit->GetParameter(1);
                xzScaleError[bin] = fit->GetParameter(2);
            }
            else
            {
                xzScaleMean[bin] = proj->GetMean();
                xzScaleError[bin] = proj->GetRMS();
            }
        }
        else
        {
            xzScaleMean[bin] = 0;
            xzScaleError[bin] = 0;
        }

        printf("%d\t%.5f\t%d\t%.5f\n",
               bin, binCenter, (Int_t)proj->GetEntries(), xzScaleMean[bin]);

        delete proj;
    }

    // Create and draw TGraphErrors for XZ angle
    TGraphErrors *grXZScale = new TGraphErrors(nProfilesXZ, xzAngleValues, xzScaleMean, 0, xzScaleError);
    grXZScale->SetMarkerStyle(20);
    grXZScale->SetMarkerColor(kBlack);
    grXZScale->SetMarkerSize(1.5);
    grXZScale->Draw("P SAME");

    // Update canvas for XZ angle
    canvas3->Modified();
    canvas3->Update();
    canvas3->Print(Form("%s", outputName));

    // Create canvas for YZ angle
    TCanvas *canvas4 = new TCanvas("canvas4", "Energy Loss Scale vs YZ Angle", 800, 600);
    canvas4->SetLeftMargin(0.16);
    canvas4->SetRightMargin(0.11);
    canvas4->SetGridx();
    canvas4->SetGridy();
    canvas4->SetLogz();

    tree->Draw("energyLossScale:TMath::ATan2(direction[1], direction[2])>>hScaleVsYZAngle", "mcBeta < 0.9");
    hScaleVsYZAngle->Draw("COLZ");

    // Prepare for profile analysis - YZ Angle
    const int nProfilesYZ = nBinsX;
    double yzAngleValues[nProfilesYZ];
    double yzScaleMean[nProfilesYZ], yzScaleError[nProfilesYZ];
    double yzBinWidth = 360.0 / nBinsX;

    // Print table header for YZ angle
    std::cout << "\nYZ Angle Bin\tAngle\tEntries\tMean\n"
              << std::string(65, '-') << std::endl;

    // Process each bin for YZ angle
    for (int bin = 0; bin < nBinsX; ++bin)
    {
        double binCenter = -180 + (bin + 0.5) * yzBinWidth;
        int binIdx = bin + 1; // ROOT histograms are 1-indexed

        // Calculate bin range for projection
        int binWindow = 1; // Use ±1 bins for better statistics
        int binLow = TMath::Max(1, binIdx - binWindow);
        int binHigh = TMath::Min(nBinsX, binIdx + binWindow);

        yzAngleValues[bin] = binCenter;

        // Get projection and fit
        TH1D *proj = hScaleVsYZAngle->ProjectionY(Form("projYZ_%d", bin), binLow, binHigh, "e");
        if (proj->GetEntries() > 10)
        {
            proj->Fit("gaus", "QN", "", energyLossScaleMin, energyLossScaleMax);
            TF1 *fit = proj->GetFunction("gaus");
            if (fit && fit->GetProb() > 0.01)
            {
                yzScaleMean[bin] = fit->GetParameter(1);
                yzScaleError[bin] = fit->GetParameter(2);
            }
            else
            {
                yzScaleMean[bin] = proj->GetMean();
                yzScaleError[bin] = proj->GetRMS();
            }
        }
        else
        {
            yzScaleMean[bin] = 0;
            yzScaleError[bin] = 0;
        }

        printf("%d\t%.5f\t%d\t%.5f\n",
               bin, binCenter, (Int_t)proj->GetEntries(), yzScaleMean[bin]);

        delete proj;
    }

    // Create and draw TGraphErrors for YZ angle
    TGraphErrors *grYZScale = new TGraphErrors(nProfilesYZ, yzAngleValues, yzScaleMean, 0, yzScaleError);
    grYZScale->SetMarkerStyle(20);
    grYZScale->SetMarkerColor(kBlack);
    grYZScale->SetMarkerSize(1.5);
    grYZScale->Draw("P SAME");

    // Update canvas for YZ angle
    canvas4->Modified();
    canvas4->Update();
    canvas4->Print(Form("%s)", outputName));

    // Close file
    file->Close();
    std::cout << "Energy loss scale plots saved to: " << outputName << std::endl;
}