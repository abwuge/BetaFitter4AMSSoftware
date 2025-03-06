#include <iostream>
#include <string>
#include <vector>

#include "amschain.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"

#include "Util.hh"
#include "BetaFitter.hh"
#include "ParticlePropagator.hh"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <inputFile> <outputFile>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];

    // Load particle data from input file
    // /home/ams/hxwu/AMSSoft/amsd69n/amsd69n/test_list_-1.root
    std::vector<ParticleData> particles = Util::loadParticleData(inputFile);
    if (particles.empty())
    {
        std::cerr << "Error: No particles loaded from input file" << std::endl;
        return 1;
    }

    // Create vectors for graph points
    std::vector<double> mcBeta;
    std::vector<double> linearBeta;
    std::vector<double> nonlinearBeta;

    TH1F *hBetaDiff = new TH1F("hBetaDiff", "1/#beta_{rec} - 1/#beta_{MC};1/#beta_{rec} - 1/#beta_{MC};Events", 100, -0.5, 0.5);

    // Draw and save plot
    TCanvas *c1 = new TCanvas("c1", "Beta Reconstruction Comparison", 800, 600);
    c1->cd();

    // Process each particle
    for (const auto &particle : particles)
    {
        // Skip invalid particles
        if (!particle.isMC)
            continue;
        // if (particle.mass < 0)
        //     continue;
        if (particle.mcCoo[0] == -1000)
            continue;

        // Setup particle propagator with initial state
        ParticlePropagator propagator(particle);

        // Prepare arrays for beta reconstruction
        double measuredTimes[4], timeErrors[4];
        for (int i = 0; i < 4; ++i)
        {
            measuredTimes[i] = particle.TOF_hitTime[i];
            timeErrors[i] = particle.TOF_hitTimeError[i];
        }

        // Reconstruct beta with nonlinear algorithm
        double nonlinear_beta_rec = 1 / BetaFitter::reconstructBeta(&particle, propagator, measuredTimes, timeErrors);
        if (nonlinear_beta_rec <= 0)
            continue;

        // Get linear beta from particle data
        double linear_beta_rec = particle.betaLinear;
        if (linear_beta_rec <= 0)
            continue;

        // Store values for plotting
        mcBeta.push_back(particle.mcBeta);
        linearBeta.push_back(linear_beta_rec);
        nonlinearBeta.push_back(nonlinear_beta_rec); // Store nonlinear beta direct

        hBetaDiff->Fill(1.0 / nonlinear_beta_rec - 1.0 / particle.mcBeta);
    }

    // Create graphs
    TGraph *grLinear = new TGraph(mcBeta.size(), &mcBeta[0], &linearBeta[0]);
    TGraph *grNonlinear = new TGraph(mcBeta.size(), &mcBeta[0], &nonlinearBeta[0]);

    // Set graph styles
    grLinear->SetMarkerStyle(20);
    grLinear->SetMarkerColor(kBlue);
    grLinear->SetLineColor(kBlue);
    grLinear->SetMarkerSize(0.5);

    grNonlinear->SetMarkerStyle(21);
    grNonlinear->SetMarkerColor(kRed);
    grNonlinear->SetLineColor(kRed);
    grNonlinear->SetMarkerSize(0.5);

    // Set up drawing
    grLinear->SetTitle("Beta Reconstruction Comparison");
    grLinear->GetXaxis()->SetTitle("MC Beta");
    grLinear->GetXaxis()->SetRangeUser(0, 1);
    grLinear->GetYaxis()->SetTitle("Reconstructed Beta");
    grLinear->GetYaxis()->SetRangeUser(0, 1);

    // Draw graphs
    grLinear->Draw("AP");
    grNonlinear->Draw("P SAME");

    // Add legend
    TLegend *legend = new TLegend(0.15, 0.7, 0.45, 0.85);
    legend->AddEntry(grLinear, "Linear Reconstruction", "p");
    legend->AddEntry(grNonlinear, "Nonlinear Reconstruction", "p");
    legend->Draw();

    // Draw diagonal reference line (perfect reconstruction)
    TGraph *grReference = new TGraph(2);
    grReference->SetPoint(0, 0, 0);
    grReference->SetPoint(1, 1, 1);
    grReference->SetLineStyle(2);
    grReference->SetLineColor(kGreen);
    grReference->Draw("L SAME");

    // Save plot
    c1->SaveAs((outputFile + "_betaComparison.pdf").c_str());

    hBetaDiff->Draw();
    c1->SaveAs((outputFile + "_betaDiff.pdf").c_str());

    // Clean up
    delete hBetaDiff;
    delete grLinear;
    delete grNonlinear;
    delete grReference;
    delete legend;
    delete c1;

    return 0;
}