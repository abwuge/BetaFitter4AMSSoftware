#include <TFile.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TLatex.h>
#include <iostream>

/**
 * Plot element abundance from AMS material scan file
 * This macro reads the g4mscan.root file and generates 2D histograms
 * of material amount for tracks parallel to z-axis
 *
 * @param sfn Path to the input ROOT file containing element abundance data
 * @param outputName Base name for output files (without extension)
 */
void plotElementAbundance(
    double zMin = 86.10,
    double zMax = 64.425,
    const char *outputName = "test_MaterialAbundanceMap.pdf",
    const char *sfn = "/cvmfs/ams.cern.ch/Offline/AMSDataDir/v6.00/LAPP/dEdxPDF/g4mscan.root")
{
    if (zMin > zMax)
        std::swap(zMin, zMax);

    // Set batch mode and style
    // ------------------------------------------------------------------------

    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);

    // Defines
    // ------------------------------------------------------------------------

    // Atomic mass of elements
    double A[9] = {1.0079, 12.0107, 14.0067, 15.9994, 18.9984, 22.9897, 26.9815, 28.0855, 207.2};

    // element names
    const char *elemNames[9] = {"H", "C", "N", "O", "F", "Na", "Al", "Si", "Pb"};

    // range
    double xMin = -60;
    double yMin = -60;
    double xMax = 60;
    double yMax = 60;

    // Open file
    // ------------------------------------------------------------------------

    TFile *file = TFile::Open(sfn);
    if (!file)
    {
        std::cerr << "Error: Cannot open file " << sfn << std::endl;
        return;
    }

    TH3F *hist[9];
    for (int i = 0; i < 9; i++)
        hist[i] = (TH3F *)file->Get(Form("hist%d", i + 1));

    // Create histograms
    // ------------------------------------------------------------------------

    TH2F *histTotal = new TH2F("totalMaterialAbundance",
                               "Total Material Abundance Map;X [cm];Y [cm];Radiation Length [g/cm^{2}]",
                               120, xMin, xMax, 120, yMin, yMax);

    TH2F *hist2D[9];
    for (int elem = 0; elem < 9; elem++)
        hist2D[elem] = new TH2F(Form("materialAbundance_%s", elemNames[elem]),
                                Form("Material Abundance Map for %s;X [cm];Y [cm];Radiation Length [g/cm^{2}]", elemNames[elem]),
                                120, xMin, xMax, 120, yMin, yMax);

    // Get bin range
    // ------------------------------------------------------------------------

    int xMinBin = hist[0]->GetXaxis()->FindBin(xMin);
    int xMaxBin = hist[0]->GetXaxis()->FindBin(xMax);

    int yMinBin = hist[0]->GetYaxis()->FindBin(yMin);
    int yMaxBin = hist[0]->GetYaxis()->FindBin(yMax);

    int zMinBin = hist[0]->GetZaxis()->FindBin(zMin);
    int zMaxBin = hist[0]->GetZaxis()->FindBin(zMax);

    // Calculate material amount for each (x,y) position
    // ------------------------------------------------------------------------

    for (int ix = xMinBin; ix <= xMaxBin; ix++)
    {
        double x = hist[0]->GetXaxis()->GetBinCenter(ix);
        for (int iy = yMinBin; iy <= yMaxBin; iy++)
        {
            double y = hist[0]->GetYaxis()->GetBinCenter(iy);

            double totalMs = 0;
            double ms[9]{};

            for (int i = 0; i < 9; ++i)
            {
                ms[i] = hist[i]->Integral(ix, ix, iy, iy, zMinBin, zMaxBin) * A[i];
                hist2D[i]->Fill(x, y, ms[i]);
                totalMs += ms[i];
            }

            histTotal->Fill(x, y, totalMs);
        }
    }

    // Draw histograms
    // ------------------------------------------------------------------------

    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "Material Abundance", 900, 800);
    c1->SetRightMargin(0.15);

    // Open PDF file
    c1->Print(Form("%s[", outputName));

    // Draw total abundance
    histTotal->SetMinimum(0);
    histTotal->Draw("COLZ");

    // Add title
    TLatex *title = new TLatex(0.5, 0.95, Form("Total Material Abundance (Z: %.2f to %.2f cm)", zMin, zMax));
    title->SetNDC();
    title->SetTextAlign(23);
    title->Draw();

    // Save total abundance image
    c1->Print(outputName);

    // Create individual element plots
    for (int elem = 0; elem < 9; elem++)
    {
        hist2D[elem]->SetMinimum(0);
        hist2D[elem]->Draw("COLZ");

        // Add title
        title->SetText(0.5, 0.95, Form("%s Material Abundance (Z: %.2f to %.2f cm)", elemNames[elem], zMin, zMax));
        title->Draw();

        // Save element abundance image
        c1->Print(outputName);
    }

    // Close PDF file
    c1->Print(Form("%s]", outputName));

    std::cout << "Combined material abundance plots generated and saved" << std::endl;
}