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

/**
 * Draw beta comparison plots from a ROOT file containing betaTree
 *
 * @param fileName Path to the input ROOT file
 * @param outputName Output file name (default: test_beta_comparison.pdf)
 */
void plotBetaComparison(const char *fileName = "test.root",
                        const char *outputName = "test_beta_comparison.pdf")
{
    // Set batch mode to avoid GUI related issues
    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetEndErrorSize(20);

    // Open the ROOT file
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return;
    }

    // Get the beta reconstruction tree
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

    // Use a more efficient approach to get beta range - default values
    double mcBetaMin = tree->GetMinimum("mcBeta");
    double mcBetaMax = tree->GetMaximum("mcBeta");

    // Define variables for reading from tree
    double mcBeta, linearBeta, nonlinearBeta;

    // Set branch addresses
    tree->SetBranchAddress("mcBeta", &mcBeta);
    tree->SetBranchAddress("linearBeta", &linearBeta);
    tree->SetBranchAddress("nonlinearBeta", &nonlinearBeta);

    // Create histograms for the 2D plots
    // Using reasonable beta range from 0.1 to 1.0
    int nBinsX = 40;                           // Number of bins in x direction
    int nBinsY = 100;                          // Number of bins in y direction for the residual plot
    double xMin = mcBetaMin, xMax = mcBetaMax; // Range for the beta values

    double yMinRes = -0.2, yMaxRes = 1.2;

    // Create histograms for the residuals plot (1/beta_rec - 1/beta_mc)
    TH2F *hNonlinearResVsMC = new TH2F("hNonlinearResVsMC",
                                       "Non-linear Reconstruction Residuals;#beta_{MC};1/#beta_{non-linear} - 1/#beta_{MC}",
                                       nBinsX, xMin, xMax, nBinsY, yMinRes, yMaxRes);
    TH2F *hLinearResVsMC = new TH2F("hLinearResVsMC",
                                    "Linear Reconstruction Residuals;#beta_{MC};1/#beta_{linear} - 1/#beta_{MC}",
                                    nBinsX, xMin, xMax, nBinsY, yMinRes, yMaxRes);

    hNonlinearResVsMC->SetMinimum(1);
    hLinearResVsMC->SetMinimum(1);

    // Fill histograms
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);

        // Calculate residuals
        double nonlinearRes = (1.0 / nonlinearBeta) - (1.0 / mcBeta);
        double linearRes = (1.0 / linearBeta) - (1.0 / mcBeta);

        // Fill histograms with residual values
        hNonlinearResVsMC->Fill(mcBeta, nonlinearRes);
        hLinearResVsMC->Fill(mcBeta, linearRes);
    }

    // Arrays to store fit results
    const int nProfiles = nBinsX;
    double mcBetaValues[nProfiles];
    double nonlinearMean[nProfiles], nonlinearError[nProfiles];
    double linearMean[nProfiles], linearError[nProfiles];

    // Print table header before the loop
    std::cout << std::endl;
    std::cout << "Bin\tBeta\t\tProj_NL\tProj_L\tNL_mean\t\tL_mean" << std::endl;
    std::cout << std::string(65, '-') << std::endl;

    double binWidth = (xMax - xMin) / nBinsX;
    for (int bin = 0; bin < nBinsX; ++bin)
    {
        double binCenter = xMin + (bin + 0.5) * binWidth;
        int binIdx = bin + 1; // ROOT histograms are 1-indexed

        // Calculate bin range for projection with wider window
        int binWindow = 1; // Use ±1 bins for better statistics
        int binLow = TMath::Max(1, binIdx - binWindow);
        int binHigh = TMath::Min(nBinsX, binIdx + binWindow);

        mcBetaValues[bin] = binCenter;

        // Process nonlinear beta residual distribution with explicit options
        TH1D *projNonlinear = hNonlinearResVsMC->ProjectionY(Form("projNL_%d", bin), binLow, binHigh, "e"); // 'e' option to keep errors
        Int_t projNLEntries = projNonlinear->GetEntries();

        if (projNLEntries) // Only fit if we have enough statistics
        {
            // Add reasonable initial parameters for the fit
            projNonlinear->Fit("gaus", "QN", "", yMinRes, yMaxRes); // Add fit range
            TF1 *fitNL = projNonlinear->GetFunction("gaus");
            if (fitNL && fitNL->GetProb() > 0.01) // Check if fit is reasonable
            {
                nonlinearMean[bin] = fitNL->GetParameter(1);
                nonlinearError[bin] = fitNL->GetParameter(2);
            }
            else
            {
                nonlinearMean[bin] = projNonlinear->GetMean(); // Use histogram mean if fit fails
                nonlinearError[bin] = projNonlinear->GetRMS();
            }
        }
        else
        {
            nonlinearMean[bin] = 0;
            nonlinearError[bin] = 0;
        }
        delete projNonlinear;

        // Process linear beta residual distribution with same improvements
        TH1D *projLinear = hLinearResVsMC->ProjectionY(Form("projL_%d", bin), binLow, binHigh, "e");
        Int_t projLEntries = projLinear->GetEntries();

        if (projLEntries)
        {
            projLinear->Fit("gaus", "QN", "", yMinRes, yMaxRes);
            TF1 *fitL = projLinear->GetFunction("gaus");
            if (fitL && fitL->GetProb() > 0.01)
            {
                linearMean[bin] = fitL->GetParameter(1);
                linearError[bin] = fitL->GetParameter(2);
            }
            else
            {
                linearMean[bin] = projLinear->GetMean();
                linearError[bin] = projLinear->GetRMS();
            }
        }
        else
        {
            linearMean[bin] = 0;
            linearError[bin] = 0;
        }
        delete projLinear;

        // Print data in aligned columns using tabs
        printf("%d\t%.6f\t%d\t%d\t%.6f\t%.6f\n",
               bin, binCenter, projNLEntries, projLEntries,
               nonlinearMean[bin], linearMean[bin]);
    }

    // Create canvas with multiple pages
    TCanvas *canvas = new TCanvas("canvas", "Beta Residuals Comparison", 3508, 2480);
    canvas->Print(Form("%s[", outputName)); // Open PDF file

    // First page - Non-linear beta residuals
    TCanvas *c1 = new TCanvas("c1", "Non-linear Beta Residuals", 3508, 2480);
    c1->SetLeftMargin(0.16);
    c1->SetRightMargin(0.11);
    c1->SetGridx();
    c1->SetGridy();
    c1->SetLogz();

    hNonlinearResVsMC->Draw("COLZ");

    // Perfect residual reference line (zero residual)
    TF1 *perfectLine1 = new TF1("perfectLine1", "0", xMin, xMax);
    perfectLine1->SetLineColor(kRed);
    perfectLine1->SetLineStyle(2);
    perfectLine1->Draw("SAME");

    // Create and draw TGraphErrors for nonlinear residuals
    TGraphErrors *grNonlinear = new TGraphErrors(nProfiles, mcBetaValues, nonlinearMean, 0, nonlinearError);
    grNonlinear->SetMarkerStyle(20);
    grNonlinear->SetMarkerColor(kBlack);
    grNonlinear->SetMarkerSize(3.0);
    grNonlinear->Draw("P");

    c1->Print(outputName);

    // Second page - Linear beta residuals
    TCanvas *c2 = new TCanvas("c2", "Linear Beta Residuals", 3508, 2480);
    c2->SetLeftMargin(0.16);
    c2->SetRightMargin(0.11);
    c2->SetGridx();
    c2->SetGridy();
    c2->SetLogz();

    hLinearResVsMC->Draw("COLZ");

    // Perfect residual reference line (zero residual)
    TF1 *perfectLine2 = new TF1("perfectLine2", "0", xMin, xMax);
    perfectLine2->SetLineColor(kRed);
    perfectLine2->SetLineStyle(2);
    perfectLine2->Draw("SAME");

    // Create and draw TGraphErrors for linear residuals
    TGraphErrors *grLinear = new TGraphErrors(nProfiles, mcBetaValues, linearMean, 0, linearError);
    grLinear->SetMarkerStyle(20);
    grLinear->SetMarkerColor(kBlack);
    grLinear->SetMarkerSize(3.0);
    grLinear->Draw("P");

    c2->Print(outputName);

    // Third page - Comparison of TGraphErrors
    TCanvas *c3 = new TCanvas("c3", "Beta Reconstruction Methods Comparison", 3508, 2480);
    c3->SetLeftMargin(0.16);
    c3->SetGridx();
    c3->SetGridy();

    // Create a frame for the comparison plot
    TH2F *hFrame = new TH2F("hFrame", "Beta Reconstruction Methods Comparison;#beta_{MC};1/#beta_{rec} - 1/#beta_{MC}",
                            100, xMin, xMax, 100, yMinRes, yMaxRes);
    hFrame->Draw();

    // Draw both graphs with different styles
    grNonlinear->SetMarkerStyle(20);
    grNonlinear->SetMarkerColor(kBlue);
    grNonlinear->SetLineColor(kBlue);
    grNonlinear->SetMarkerSize(3.0);
    grNonlinear->Draw("LP");

    grLinear->SetMarkerStyle(21);
    grLinear->SetMarkerColor(kRed);
    grLinear->SetLineColor(kRed);
    grLinear->SetMarkerSize(3.0);
    grLinear->Draw("LP");

    // Add legend
    TLegend *legend = new TLegend(0.55, 0.7, 0.85, 0.85);
    legend->AddEntry(grNonlinear, "Non-linear Method", "lp");
    legend->AddEntry(grLinear, "Linear Method", "lp");
    legend->SetBorderSize(0);
    legend->Draw();

    // Draw zero line
    TF1 *zeroLine = new TF1("zeroLine", "0", xMin, xMax);
    zeroLine->SetLineStyle(2);
    zeroLine->SetLineColor(kGray + 2);
    zeroLine->Draw("SAME");

    c3->Print(outputName);
    c3->Print(Form("%s]", outputName)); // Close PDF file

    // Clean up
    delete zeroLine;
    delete legend;
    delete hFrame;
    delete c3;
    delete grNonlinear;
    delete grLinear;
    delete perfectLine1;
    delete perfectLine2;
    delete c2;
    delete c1;
    delete hNonlinearResVsMC;
    delete hLinearResVsMC;
    delete canvas;

    // Reset tree branches to avoid dangling pointers
    tree->ResetBranchAddresses();

    // Close and delete the file after all references are gone
    file->Close();
    delete file;

    std::cout << "Beta residuals comparison plot saved to: " << outputName << std::endl;
}
