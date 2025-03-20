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
#include <fstream>
#include <vector>
#include <TFitResult.h>
#include <string>
#include <array>
#include <algorithm>

const int nBinsX = 40;
const int nBinsY = 100;
const double yMin = -0.1;
const double yMax = 0.4;
const int colors[] = {
    kRed,
    kBlue,
    kGreen + 2,
    kMagenta,
    kCyan + 2,
    kOrange + 7,
    kViolet - 3,
    kPink + 2};

double xMin = 0;
double xMax = 0;
double binWidth = 0;

/**
 * Get Z value and energyLossScale for a file from README.md
 * @param fileName Name of the ROOT file
 * @return pair of Z value and energyLossScale. Returns {-1, -1} if not found
 */
double getParamFromReadme(const std::string &fileName)
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
    {
        targetFile = targetFile.substr(8);
    }
    else
    {
        std::cerr << "Warning: File name does not start with 'results/'"
                  << "\nCannot find the file in README.md" << std::endl;
        return -1;
    }

    double energyLossScale = -1.0;

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
            // Get energyLossScale value
            size_t scalePos = it->find("energyLossScale = ", filePos);
            if (scalePos != std::string::npos)
            {
                scalePos += 17; // Skip "energyLossScale = "
                size_t commaPos = it->find(",", scalePos);
                if (commaPos != std::string::npos)
                {
                    std::string scaleStr = it->substr(scalePos, commaPos - scalePos);
                    try
                    {
                        energyLossScale = std::stod(scaleStr);
                    }
                    catch (...)
                    {
                        continue;
                    }
                }
            }

            if (energyLossScale > 0)
                break;
        }
    }

    return energyLossScale;
}

/**
 * Structure to hold file information and data
 */
struct Data
{
    char *legend{nullptr};
    TH2D *hist{nullptr};
    TGraph *graph{nullptr};
    TGraph *graphErrors{nullptr};
};

Data getData(std::string fileName, const std::string branch = "nonlinearBeta", const char *const legend = nullptr)
{
    Data data;

    // Open file
    // ------------------------------------------------------------------------

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
            return {};
        }
    }

    // Get tree
    // ------------------------------------------------------------------------

    TTree *tree = (TTree *)file->Get("betaTree");
    if (!tree)
    {
        std::cerr << "Error: Could not find betaTree in " << fileName << std::endl;
        return {};
    }

    if (tree->GetEntries() == 0)
    {
        std::cerr << "Error: Empty tree in " << fileName << ", no data to plot" << std::endl;
        return {};
    }

    // Set legend
    // ------------------------------------------------------------------------

    if (legend)
        data.legend = strdup(legend);
    else
        data.legend = strdup(Form("#zeta = %.3f", getParamFromReadme(fileName)));

    // Set histogram
    // ------------------------------------------------------------------------

    static int histogramCounter = 0;

    if (histogramCounter == 0)
    {
        xMin = tree->GetMinimum("mcBeta");
        xMax = tree->GetMaximum("mcBeta");
        binWidth = (xMax - xMin) / nBinsX;
    }

    data.hist = new TH2D(Form("hResVsMC_%d", histogramCounter),
                         "Reconstruction Residuals",
                         nBinsX, xMin, xMax, nBinsY, yMin, yMax);

    tree->Draw(
        Form("1/%s - 1/mcBeta:mcBeta>>hResVsMC_%d",
             branch.c_str(), histogramCounter),
        "",
        "goff");

    // Fit
    // ------------------------------------------------------------------------
    
    // Arrays to store fit results
    std::vector<double> mcBetaValues;
    std::vector<double> means, errors;
    mcBetaValues.reserve(nBinsX);
    means.reserve(nBinsX);
    errors.reserve(nBinsX);

    // Fit each column
    for (int i = 0; i < nBinsX; ++i)
    {
        double binCenter = xMin + (i + 0.5) * binWidth;
        TH1D *proj = data.hist->ProjectionY(
            Form("proj_%d_%d", histogramCounter, i),
            i + 1,
            i + 1);
        if (proj->GetEntries() > 50)
        {
            TFitResultPtr fitResult = proj->Fit("gaus", "SQNR");
            if (fitResult->Status() == 0)
            {
                mcBetaValues.push_back(binCenter);
                means.push_back(fitResult->Parameter(1));
                errors.push_back(fitResult->Parameter(2));
            }
        }
    }

    // Set graphs
    //  ------------------------------------------------------------------------

    // Create TGraph for residuals
    int nPoints = mcBetaValues.size();
    data.graph = new TGraph(nPoints, mcBetaValues.data(), means.data());
    data.graph->SetMarkerStyle(20);
    data.graph->SetMarkerSize(2.0);
    data.graph->SetMarkerColor(colors[histogramCounter % 8]);
    data.graph->SetLineColor(colors[histogramCounter % 8]);

    // Create TGraph for errors
    data.graphErrors = new TGraph(nPoints, mcBetaValues.data(), errors.data());
    data.graphErrors->SetMarkerStyle(20);
    data.graphErrors->SetMarkerSize(2.0);
    data.graphErrors->SetMarkerColor(colors[histogramCounter % 8]);
    data.graphErrors->SetLineColor(colors[histogramCounter % 8]);

    // Increment histogram counter
    histogramCounter++;

    return data;
}

/**
 * Draw beta comparison plots from multiple ROOT files containing betaTree
 *
 * @param fileNames Array of paths to the input ROOT files
 * @param outputName Output file name (default: "test_multiple_beta_comparison.pdf")
 */
void plotMultipleBetaComparison(const std::vector<std::string> &fileNames,
                                const char *outputName = "test_multiple_beta_comparison.pdf")
{
    if (fileNames.empty())
    {
        std::cerr << "Error: No input files provided." << std::endl;
        return;
    }

    gROOT->SetBatch(true);
    gROOT->SetStyle("Pub");
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetEndErrorSize(20);

    // Open files and prepare data
    // ------------------------------------------------------------------------

    std::vector<Data> dataVector;
    dataVector.reserve(fileNames.size() + 1);

    dataVector.push_back(getData(fileNames[0], "linearBeta", "Linear"));
    for (auto &fileName : fileNames)
        dataVector.push_back(getData(fileName));

    if (dataVector.empty())
    {
        std::cerr << "Error: No valid files to process." << std::endl;
        return;
    }

    // Create canvas and draw comparison
    // ------------------------------------------------------------------------

    TCanvas *canvas = new TCanvas("canvas", "", 3508, 2480); // A4 landscape
    canvas->SetLeftMargin(0.16);
    canvas->SetRightMargin(0.12);
    canvas->SetGridx();
    canvas->SetGridy();
    canvas->Print(Form("%s[", outputName));

    // Create comparison plot
    TH2F *hComparison = new TH2F("hComparison",
                                 ";#beta_{MC};Residual Means",
                                 nBinsX, xMin, xMax, nBinsY, yMin, yMax);

    // Create perfect residual reference line (zero residual)
    TF1 *perfectResidual = new TF1("perfectResidual", "0", xMin, xMax);
    perfectResidual->SetLineColor(kBlack);
    perfectResidual->SetLineStyle(2);

    // Create legend
    TLegend *legend = new TLegend(0.65, 0.65, 0.85, 0.85);
    legend->SetBorderSize(0);

    // Draw all residuals
    hComparison->Draw();
    perfectResidual->Draw("SAME");
    for (auto &data : dataVector)
    {
        data.graph->Draw("LP SAME");

        legend->AddEntry(data.graph, data.legend, "lp");
    }
    legend->Draw();
    canvas->Print(outputName);

    // Draw all errors
    hComparison->GetYaxis()->SetTitle("Residual Errors");
    hComparison->SetAxisRange(0, 0.1, "Y");
    hComparison->Draw();
    for (auto &data : dataVector)
        data.graphErrors->Draw("LP SAME");
    legend->Draw();
    canvas->Print(outputName);

    canvas->Print(Form("%s]", outputName));

    std::cout << "Multiple beta comparison plot saved to: " << outputName << std::endl;
}

/**
 * Convenience wrapper that accepts exactly two filenames
 */
void plotMultipleBetaComparison(const char *fileName1, const char *fileName2,
                                const char *outputName = "test_multiple_beta_comparison.pdf")
{
    std::vector<std::string> fileNames = {fileName1, fileName2};
    plotMultipleBetaComparison(fileNames, outputName);
}

/**
 * Convenience wrapper that accepts exactly three filenames
 */
void plotMultipleBetaComparison(const char *fileName1, const char *fileName2, const char *fileName3,
                                const char *outputName = "test_multiple_beta_comparison.pdf")
{
    std::vector<std::string> fileNames = {fileName1, fileName2, fileName3};
    plotMultipleBetaComparison(fileNames, outputName);
}

void plotMultipleBetaComparison(std::initializer_list<const char *> fileNames, const char *outputName = "test_multiple_beta_comparison.pdf")
{
    std::vector<std::string> fileNamesString;
    for (const char *fileName : fileNames)
        fileNamesString.push_back(fileName);

    plotMultipleBetaComparison(fileNamesString, outputName);
}