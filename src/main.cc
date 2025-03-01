#include <iostream>
#include <string>
#include <vector>

#include "amschain.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"

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

    // AMSChain chain;
    // if (!Util::addInputFile(inputFile, chain))
    // {
    //     std::cerr << "Error: could not add input file" << std::endl;
    //     return 1;
    // }

    // DataProcessor* processor = new DataProcessor(outputFile);
    // processor->processEvents(chain, maxEventsNumber);

    // delete processor;

    // Load particle data from input file
    std::vector<ParticleData> particles = Util::loadParticleData(outputFile);
    if (particles.empty()) {
        std::cerr << "Error: No particles loaded from input file" << std::endl;
        return 1;
    }

    // Create output file and histogram
    TH1F *hBetaDiff = new TH1F("hBetaDiff", "1/#beta_{rec} - 1/#beta_{MC};1/#beta_{rec} - 1/#beta_{MC};Events", 100, -0.5, 0.5);

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
        
        // Reconstruct beta
        double beta_rec = BetaFitter::reconstructBeta(&tempParticle, propagator, measuredTimes, timeErrors);
        if (beta_rec <= 0) continue;

        // Fill histogram with difference
        double inv_beta_diff = beta_rec - 1.0/particle.mcBeta; // beta_rec已经是1/beta了
        hBetaDiff->Fill(inv_beta_diff);
    }

    // Draw and save histogram
    TCanvas *c1 = new TCanvas("c1", "Beta Difference", 800, 600);
    gStyle->SetOptStat(111111);
    hBetaDiff->Draw();
    c1->SaveAs((outputFile + "_betaDiff.pdf").c_str());

    return 0;
}