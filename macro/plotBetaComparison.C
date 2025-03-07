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

    // Use a more efficient approach to get beta range - default values
    double mcBetaMin = 0.1;
    double mcBetaMax = 1.0;
    
    // Get actual beta range from the tree using ROOT's statistical functions
    if (tree->GetEntries() > 0) {
        // Create a temporary histogram to get min/max efficiently
        // This won't be displayed and will be deleted right after use
        TH1D* tempHist = new TH1D("tempHist", "Temporary Histogram", 1000, 0, 0);
        
        // Fill the histogram with mcBeta values
        tree->Draw("mcBeta>>tempHist", "", "goff");
        
        // Get min and max from histogram statistics
        if (tempHist->GetEntries() > 0) {
            // For safety, in case of outliers, we take a range that covers most data points
            mcBetaMin = tempHist->GetXaxis()->GetXmin();
            mcBetaMax = tempHist->GetXaxis()->GetXmax();
            
            // Add a small margin
            double margin = (mcBetaMax - mcBetaMin) * 0.05;
            mcBetaMin = TMath::Max(0.1, mcBetaMin - margin);  
            mcBetaMax = TMath::Min(1.0, mcBetaMax + margin);
        }
        
        // Clean up temporary histogram
        delete tempHist;
    }

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

    // // Set color palette for better visualization
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);

    // Create canvas
    TCanvas *canvas = new TCanvas("canvas", "Beta Residuals Comparison", 800, 600);
    canvas->Divide(1, 2);

    // Arrays to store fit results
    const int nProfiles = nBinsX;
    double mcBetaValues[nProfiles];
    double nonlinearMean[nProfiles], nonlinearError[nProfiles];
    double linearMean[nProfiles], linearError[nProfiles];

    // Process upper plot (nonlinear beta residuals)
    canvas->cd(1);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogz();

    hNonlinearResVsMC->Draw("COLZ");

    // Perfect residual reference line (zero residual)
    TF1 *perfectLine1 = new TF1("perfectLine1", "0", xMin, xMax);
    perfectLine1->SetLineColor(kRed);
    perfectLine1->SetLineStyle(2);
    perfectLine1->SetLineWidth(2);
    perfectLine1->Draw("SAME");

    // Process lower plot (linear beta residuals)
    canvas->cd(2);
    gPad->SetGridx();
    gPad->SetGridy();
    gPad->SetLogz();

    hLinearResVsMC->Draw("COLZ");

    // Perfect residual reference line (zero residual)
    TF1 *perfectLine2 = new TF1("perfectLine2", "0", xMin, xMax);
    perfectLine2->SetLineColor(kRed);
    perfectLine2->SetLineStyle(2);
    perfectLine2->SetLineWidth(2);
    perfectLine2->Draw("SAME");

    // Print table header before the loop
    std::cout << std::endl;
    std::cout << "Bin\tBeta\t\tProj_NL\t\tProj_L\t\tNL_mean\t\tL_mean" << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    double binWidth = (xMax - xMin) / nBinsX;
    for (int bin = 0; bin < nBinsX; ++bin)
    {
        double binCenter = xMin + (bin + 0.5) * binWidth;
        int binIdx = bin + 1; // ROOT histograms are 1-indexed

        // Calculate bin range for projection with wider window
        int binWindow = 1; // Use ±1 bins for better statistics
        int binLow = std::max(1, binIdx - binWindow);
        int binHigh = std::min(nBinsX, binIdx + binWindow);

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
        printf("%d\t%.6f\t%d\t\t%d\t\t%.6f\t%.6f\n",
               bin, binCenter, projNLEntries, projLEntries,
               nonlinearMean[bin], linearMean[bin]);
    }

    // Create and draw TGraphErrors for the Gaussian means of the residuals
    canvas->cd(1);
    TGraphErrors *grNonlinear = new TGraphErrors(nProfiles, mcBetaValues, nonlinearMean, 0, nonlinearError);
    grNonlinear->SetMarkerStyle(20);
    grNonlinear->SetMarkerSize(0.8);
    grNonlinear->SetMarkerColor(kBlack);
    grNonlinear->Draw("P");

    canvas->cd(2);
    TGraphErrors *grLinear = new TGraphErrors(nProfiles, mcBetaValues, linearMean, 0, linearError);
    grLinear->SetMarkerStyle(20);
    grLinear->SetMarkerSize(0.8);
    grLinear->SetMarkerColor(kBlack);
    grLinear->Draw("P");

    // Save the canvas
    canvas->SaveAs(outputName);

    // Clean up - make sure to delete objects in the correct order
    delete grNonlinear;
    delete grLinear;
    delete perfectLine1;
    delete perfectLine2;
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
