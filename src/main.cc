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
#include "DataProcessor.hh"
#include "BetaFitter.hh"
#include "ParticlePropagator.hh"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <inputFile> <outputFile.root> <maxEventsNumber>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];
    int maxEventsNumber = argc > 3 ? std::stoi(argv[3]) : -1;

    AMSChain chain;
    if (!Util::addInputFile(inputFile, chain))
    {
        std::cerr << "Error: could not add input file" << std::endl;
        return 1;
    }

    DataProcessor* processor = new DataProcessor(outputFile);
    processor->processEvents(chain, maxEventsNumber);

    delete processor;

    // Load particle data from input file
    std::vector<ParticleData> particles = Util::loadParticleData(outputFile);
    if (particles.empty()) {
        std::cerr << "Error: No particles loaded from input file" << std::endl;
        return 1;
    }

    // Create vectors for graph points
    std::vector<double> mcBeta;
    std::vector<double> linearBeta;
    std::vector<double> nonlinearBeta;

    // Process each particle
    for (const auto& particle : particles) {
        // Skip invalid particles
        if (particle.momentum == 0 || !particle.isMC) continue;
        if (particle.mass < 0) continue;

        // Create initial state for propagation
        AMSPoint pos(particle.hitX[0], particle.hitY[0], particle.hitZ[0]);
        AMSDir dir;
        dir.SetTheta(particle.Theta);
        dir.SetPhi(particle.Phi);

        // Setup particle propagator with initial state
        ParticlePropagator propagator(pos, dir, particle.momentum, particle.mass, particle.charge);

        // Prepare arrays for beta reconstruction
        double measuredTimes[4], timeErrors[4];
        for (int i = 0; i < 4; ++i) {
            measuredTimes[i] = particle.hitTime[i];
            timeErrors[i] = particle.hitTimeError[i];
        }

        // Create temporary ParticleR object for reconstruction
        ParticleR tempParticle;
        tempParticle.Mass = particle.mass;
        tempParticle.Charge = particle.charge;
        tempParticle.Momentum = particle.momentum;
        tempParticle.Theta = particle.Theta;
        tempParticle.Phi = particle.Phi;
        for (int i = 0; i < 4; ++i) {
            tempParticle.TOFCoo[i][0] = particle.hitX[i];
            tempParticle.TOFCoo[i][1] = particle.hitY[i];
            tempParticle.TOFCoo[i][2] = particle.hitZ[i];
            tempParticle.TOFTLength[i] = std::sqrt(pow(particle.hitX[i] - particle.hitX[0], 2) + 
                                                  pow(particle.hitY[i] - particle.hitY[0], 2) +
                                                  pow(particle.hitZ[i] - particle.hitZ[0], 2));
        }
        
        // Reconstruct beta with nonlinear algorithm
        double nonlinear_beta_rec = BetaFitter::reconstructBeta(&tempParticle, propagator, measuredTimes, timeErrors);
        if (nonlinear_beta_rec <= 0) continue;
        
        // Get linear beta from particle data
        double linear_beta_rec = particle.beta;
        if (linear_beta_rec <= 0) continue;
        
        // Store values for plotting
        mcBeta.push_back(particle.mcBeta);
        linearBeta.push_back(linear_beta_rec);
        nonlinearBeta.push_back(1.0/nonlinear_beta_rec); // Convert from 1/beta to beta
    }

    // Draw and save plot
    TCanvas *c1 = new TCanvas("c1", "Beta Reconstruction Comparison", 800, 600);
    
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
    grLinear->GetYaxis()->SetTitle("Reconstructed Beta");
    
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
    
    // Clean up
    delete grLinear;
    delete grNonlinear;
    delete grReference;
    delete legend;
    delete c1;

    return 0;
}