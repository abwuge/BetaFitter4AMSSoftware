#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <map>
#include <vector>
#include <utility>

/**
 * Plot magnetic field projections on XZ and YZ planes passing through specified points
 * The script creates histograms showing field components parallel and perpendicular to these planes
 * Uses distance-weighted interpolation from nearest points to improve accuracy
 * @param xzPlaneY Y coordinate of XZ plane
 * @param yzPlaneX X coordinate of YZ plane
 * @param inputFile Path to ROOT file containing magnetic field data
 * @param outFileName Path to save output PDF file
 */
void plotMagneticField(Double_t xzPlaneY = 0.0,
                       Double_t yzPlaneX = 0.0,
                       const char *inputFile = "test.root",
                       TString outFileName = "test_magnetic_field_projections.pdf")
{
    // Use batch mode to avoid GUI issues
    gROOT->SetBatch(true);

    // Open the file containing magnetic field data
    TFile *file = TFile::Open(inputFile);
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << inputFile << std::endl;
        return;
    }

    // Get the magnetic field tree
    TTree *magTree = (TTree *)file->Get("magfield");
    if (!magTree)
    {
        std::cerr << "Error: Could not find magfield tree in " << inputFile << std::endl;
        file->Close();
        return;
    }

    // Get coordinate ranges directly from the tree
    Double_t xMin = magTree->GetMinimum("x");
    Double_t xMax = magTree->GetMaximum("x");
    Double_t yMin = magTree->GetMinimum("y");
    Double_t yMax = magTree->GetMaximum("y");
    Double_t zMin = magTree->GetMinimum("z");
    Double_t zMax = magTree->GetMaximum("z");

    // Clamp the plane positions to valid ranges
    xzPlaneY = TMath::Range(yMin, yMax, xzPlaneY);
    yzPlaneX = TMath::Range(xMin, xMax, yzPlaneX);

    Double_t x, y, z, bx, by, bz;
    magTree->SetBranchAddress("x", &x);
    magTree->SetBranchAddress("y", &y);
    magTree->SetBranchAddress("z", &z);
    magTree->SetBranchAddress("bx", &bx);
    magTree->SetBranchAddress("by", &by);
    magTree->SetBranchAddress("bz", &bz);

    // Calculate step sizes (assuming uniform grid)
    Int_t nEntries = magTree->GetEntries();
    Double_t dx = (xMax - xMin) / TMath::Nint(TMath::Power(nEntries, 1. / 3));
    Double_t dy = (yMax - yMin) / TMath::Nint(TMath::Power(nEntries, 1. / 3));
    Double_t dz = (zMax - zMin) / TMath::Nint(TMath::Power(nEntries, 1. / 3));

    Int_t nBins = TMath::Nint((xMax - xMin) / dx); // Number of bins should match the data grid

    // Create histograms with ranges from data
    TString xzTitleParallel = TString::Format("B_{parallel} on XZ Plane (y=%.1f);X (cm);Z (cm)", xzPlaneY);
    TH2D *hXZ_parallel = new TH2D("hXZ_parallel", xzTitleParallel.Data(),
                                  nBins, xMin, xMax, nBins, zMin, zMax);

    TString xzTitlePerp = TString::Format("B_{perpendicular} on XZ Plane (y=%.1f);X (cm);Z (cm)", xzPlaneY);
    TH2D *hXZ_perp = new TH2D("hXZ_perp", xzTitlePerp.Data(),
                              nBins, xMin, xMax, nBins, zMin, zMax);

    TString yzTitleParallel = TString::Format("B_{parallel} on YZ Plane (x=%.1f);Y (cm);Z (cm)", yzPlaneX);
    TH2D *hYZ_parallel = new TH2D("hYZ_parallel", yzTitleParallel.Data(),
                                  nBins, yMin, yMax, nBins, zMin, zMax);

    TString yzTitlePerp = TString::Format("B_{perpendicular} on YZ Plane (x=%.1f);Y (cm);Z (cm)", yzPlaneX);
    TH2D *hYZ_perp = new TH2D("hYZ_perp", yzTitlePerp.Data(),
                              nBins, yMin, yMax, nBins, zMin, zMax);

    // Define data structures to store closest points for interpolation
    // Map key is (x,z) or (y,z) bin, value is vector of (y or x distance, field components)
    std::map<std::pair<Int_t, Int_t>, std::vector<std::tuple<Double_t, Double_t, Double_t>>> xzPlaneData;
    std::map<std::pair<Int_t, Int_t>, std::vector<std::tuple<Double_t, Double_t, Double_t>>> yzPlaneData;

    // Instead of processing all entries, process only points near the target planes
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        magTree->GetEntry(i);

        // Check if this point is within one step size of either plane
        bool nearXZPlane = std::abs(y - xzPlaneY) <= dy;
        bool nearYZPlane = std::abs(x - yzPlaneX) <= dx;

        if (!nearXZPlane && !nearYZPlane)
            continue;

        // Get field magnitudes
        if (nearXZPlane)
        {
            Double_t b_parallel_xz = TMath::Sqrt(bx * bx + bz * bz);
            Double_t b_perp_xz = std::abs(by);
            Int_t binX = hXZ_parallel->GetXaxis()->FindBin(x);
            Int_t binZ = hXZ_parallel->GetYaxis()->FindBin(z);
            xzPlaneData[std::make_pair(binX, binZ)].push_back(std::make_tuple(std::abs(y - xzPlaneY), b_parallel_xz, b_perp_xz));
        }

        if (nearYZPlane)
        {
            Double_t b_parallel_yz = TMath::Sqrt(by * by + bz * bz);
            Double_t b_perp_yz = std::abs(bx);
            Int_t binY = hYZ_parallel->GetXaxis()->FindBin(y);
            Int_t binZ = hYZ_parallel->GetYaxis()->FindBin(z);
            yzPlaneData[std::make_pair(binY, binZ)].push_back(std::make_tuple(std::abs(x - yzPlaneX), b_parallel_yz, b_perp_yz));
        }
    }

    // Function to perform distance-weighted interpolation
    auto interpolateField = [](const std::vector<std::tuple<Double_t, Double_t, Double_t>> &points,
                               Double_t &parallel, Double_t &perpendicular)
    {
        // If no data points, return zero
        if (points.empty())
        {
            parallel = 0.0;
            perpendicular = 0.0;
            return false;
        }

        // If only one point, use it directly
        if (points.size() == 1)
        {
            parallel = std::get<1>(points[0]);
            perpendicular = std::get<2>(points[0]);
            return true;
        }

        // Find the two closest points to the plane
        std::sort(const_cast<std::vector<std::tuple<Double_t, Double_t, Double_t>> &>(points).begin(),
                  const_cast<std::vector<std::tuple<Double_t, Double_t, Double_t>> &>(points).end(),
                  [](const std::tuple<Double_t, Double_t, Double_t> &a,
                     const std::tuple<Double_t, Double_t, Double_t> &b)
                  {
                      return std::get<0>(a) < std::get<0>(b);
                  });

        // Get the two closest points
        auto closest1 = points[0];
        auto closest2 = points[1];

        // Calculate distance-based weights
        Double_t d1 = std::get<0>(closest1);
        Double_t d2 = std::get<0>(closest2);

        // If both points are at the same distance, take average
        if (TMath::AreEqualAbs(d1, d2, 1e-10))
        {
            parallel = (std::get<1>(closest1) + std::get<1>(closest2)) / 2.0;
            perpendicular = (std::get<2>(closest1) + std::get<2>(closest2)) / 2.0;
            return true;
        }

        if (d1 <= 1e-10)
        {
            parallel = std::get<1>(closest1);
            perpendicular = std::get<2>(closest1);
            return true;
        }

        if (d2 <= 1e-10)
        {
            parallel = std::get<1>(closest2);
            perpendicular = std::get<2>(closest2);
            return true;
        }

        // Weight is inversely proportional to distance
        Double_t w1 = 1.0 / d1;
        Double_t w2 = 1.0 / d2;
        Double_t totalWeight = w1 + w2;

        // Weighted average
        parallel = (w1 * std::get<1>(closest1) + w2 * std::get<1>(closest2)) / totalWeight;
        perpendicular = (w1 * std::get<2>(closest1) + w2 * std::get<2>(closest2)) / totalWeight;

        return true;
    };

    // Fill histograms with interpolated field values
    for (int binx = 1; binx <= hXZ_parallel->GetNbinsX(); binx++)
    {
        for (int binz = 1; binz <= hXZ_parallel->GetNbinsY(); binz++)
        {
            // Interpolate for XZ plane
            Double_t xz_parallel, xz_perp;
            if (interpolateField(xzPlaneData[std::make_pair(binx, binz)], xz_parallel, xz_perp))
            {
                hXZ_parallel->SetBinContent(binx, binz, xz_parallel);
                hXZ_perp->SetBinContent(binx, binz, xz_perp);
            }

            // Interpolate for YZ plane
            Double_t yz_parallel, yz_perp;
            if (interpolateField(yzPlaneData[std::make_pair(binx, binz)], yz_parallel, yz_perp))
            {
                hYZ_parallel->SetBinContent(binx, binz, yz_parallel);
                hYZ_perp->SetBinContent(binx, binz, yz_perp);
            }
        }
    }

    // Style settings
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);

    // Create individual canvases for each plot to prevent memory issues
    TCanvas *c1 = new TCanvas("c1", "XZ Parallel", 800, 600);
    hXZ_parallel->SetTitle(xzTitleParallel);
    hXZ_parallel->Draw("COLZ");
    c1->Print(outFileName + "("); // Open PDF file with first plot
    delete c1;                    // Delete canvas to free memory

    TCanvas *c2 = new TCanvas("c2", "XZ Perpendicular", 800, 600);
    hXZ_perp->SetTitle(xzTitlePerp);
    hXZ_perp->Draw("COLZ");
    c2->Print(outFileName);
    delete c2;

    TCanvas *c3 = new TCanvas("c3", "YZ Parallel", 800, 600);
    hYZ_parallel->SetTitle(yzTitleParallel);
    hYZ_parallel->Draw("COLZ");
    c3->Print(outFileName);
    delete c3;

    TCanvas *c4 = new TCanvas("c4", "YZ Perpendicular", 800, 600);
    hYZ_perp->SetTitle(yzTitlePerp);
    hYZ_perp->Draw("COLZ");
    c4->Print(outFileName + ")"); // Close PDF file with last plot
    delete c4;

    std::cout << "Plots saved to " << outFileName << std::endl;

    // Clean up
    delete hXZ_parallel;
    delete hXZ_perp;
    delete hYZ_parallel;
    delete hYZ_perp;
    file->Close();
    delete file;
}