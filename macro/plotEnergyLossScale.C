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
#include <fstream>
#include <vector>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TFitResult.h>

/**
 * Get Z value for a file from README.md
 * @param fileName Name of the ROOT file
 * @return Z value. Returns -1 if not found
 */
int getParamsFromReadme(const std::string &fileName)
{
    std::ifstream readme("results/README.md");
    if (!readme)
    {
        std::cerr << "Error: Could not open results/README.md" << std::endl;
        return -1;
    }
    std::string line;
    std::string targetFile = fileName;

    // If fileName starts with "results/", remove it
    if (targetFile.substr(0, 8) == "results/")
        targetFile = targetFile.substr(8);
    else
    {
        std::cerr << "Warning: File name does not start with 'results/'"
                  << "\nCannot find the file in README.md" << std::endl;
        return -1;
    }

    int latestZ = -1;
    // Read file from bottom to top to get the latest entry
    std::vector<std::string> lines;
    while (std::getline(readme, line))
        lines.push_back(line);

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
            if (latestZ > 0)
                break;
        }
    }
    return latestZ;
}

/**
 * Draw energy loss scale distribution from a ROOT file containing scaleTree tree
 *
 * @param fileName Path to the input ROOT file
 * @param outputName Output file name
 */
void plotEnergyLossScale(std::string fileName = "test.root",
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

    // Get Z value and energyLossScale from the file
    int zValue = getParamsFromReadme(fileName);

    std::string actualOutputName;
    if (!outputName)
    {
        if (zValue > 0)
            actualOutputName = Form("test_zeta_Z%d.pdf", zValue);
        else
            actualOutputName = "test_zeta.pdf";
        outputName = actualOutputName.c_str();
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
    double energyLossScaleMin = -6;
    double energyLossScaleMax = 10;
    double mcBetaMin = tree->GetMinimum("mcBeta");
    double mcBetaMax = tree->GetMaximum("mcBeta");

    // Number of bins for histograms
    int nBins = 100;
    int nBinsX = 40;   // Number of bins in beta direction
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
    canvas1->SetLogy();

    tree->Draw("energyLossScale>>hEnergyLossScale", "mcBeta < 0.95");
    hEnergyLossScale->SetLineColor(kBlue);
    hEnergyLossScale->SetLineWidth(2);
    hEnergyLossScale->SetFillColor(kBlue - 10);
    hEnergyLossScale->SetFillStyle(3004);
    hEnergyLossScale->Draw();

    TPaveText *infoText = new TPaveText(0.2, 0.92, 0.8, 0.98, "NDC");
    infoText->SetFillColor(0);
    infoText->SetBorderSize(0);

    if (zValue > 0)
        infoText->AddText(Form("Z = %d, mcBeta < 0.9", zValue));
    else
        infoText->AddText(Form("mcBeta < 0.9"));

    infoText->Draw();

    TLegend *legend = new TLegend(0.62, 0.67, 0.87, 0.87);
    legend->SetBorderSize(kNone);
    legend->AddEntry("", Form("#mu = %.4g", hEnergyLossScale->GetMean()), "");
    legend->AddEntry("", Form("#sigma = %.4g", hEnergyLossScale->GetRMS()), "");
    legend->AddEntry("", Form("Entries: %d", (int)hEnergyLossScale->GetEntries()), "");
    legend->Draw();

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

    std::vector<double> betaValues, scaleMean, scaleError;
    betaValues.reserve(nBinsX);
    scaleMean.reserve(nBinsX);
    scaleError.reserve(nBinsX);
    double binWidth = (mcBetaMax - mcBetaMin) / nBinsX;

    // Process each bin
    TF1 *fGaus = new TF1("fGaus", "gaus", energyLossScaleMin, energyLossScaleMax);
    fGaus->SetParLimits(1, energyLossScaleMin, energyLossScaleMax);
    for (int i = 0; i < nBinsX; ++i)
    {
        double binCenter = mcBetaMin + (i + 0.5) * binWidth;

        // Get projection and fit
        TH1D *proj = hScaleVsBeta->ProjectionY(Form("proj_%d", i), i + 1, i + 1);
        if (proj->GetEntries() > 50)
        {
            TFitResultPtr fitResult = proj->Fit(fGaus, "SQNR");
            if (fitResult->Status() == 0)
            {
                betaValues.push_back(binCenter);
                scaleMean.push_back(fitResult->Parameter(1));
                scaleError.push_back(fitResult->Parameter(2));
            }

            delete proj;
        }
    }

    // Create and draw TGraphErrors
    TGraphErrors *grScale = new TGraphErrors(betaValues.size(), betaValues.data(), scaleMean.data(), 0, scaleError.data());
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