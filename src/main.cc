#include <iostream>
#include <string>
#include <vector>

#include "amschain.h"
#include "TFile.h"
#include "TTree.h"

#include "Util.hh"
#include "BetaFitter.hh"
#include "ParticlePropagator.hh"

int main(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " <inputFile.root> <outputFile.root> [<fitOption>]" << std::endl;
        std::cout << "fitOption: " << std::endl;
        std::cout << "  -1: Save magnetic field information to ROOT file" << std::endl;
        std::cout << "   0: Only TOF hits" << std::endl;
        std::cout << "   1: TOF hits + Tracker hits" << std::endl;
        std::cout << "   2: TOF hits + Tracker hits + Energy loss scale" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    std::string outputFile = argv[2];

    /**
     * Fit option:
     * -1: Save magnetic field information
     *  0: Only TOF hits
     *  1: TOF hits + Tracker hits
     *  2: TOF hits + Tracker hits + Energy loss scale
     */
    BetaFitter::fitOption = argc > 3 ? atoi(argv[3]) : 0;

    // Handle magnetic field information saving
    if (BetaFitter::fitOption == -1) {
        return Util::saveMagneticField(outputFile) ? 0 : 1;
    }

    // Normal beta reconstruction workflow
    if (BetaFitter::fitOption < 0 || BetaFitter::fitOption > 2)
        BetaFitter::fitOption = 0;

    // Load particle data from input file
    std::vector<ParticleData> particles = Util::loadParticleData(inputFile);
    if (particles.empty())
    {
        std::cerr << "Error: No particles loaded from input file" << std::endl;
        return 1;
    }

    // Create output ROOT file
    TFile *outFile = new TFile(outputFile.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie())
    {
        std::cerr << "Error: Could not create output file" << std::endl;
        return 1;
    }

    // Create TTree
    TTree *tree = new TTree("betaTree", "Beta Reconstruction Results");
    double mcBeta = 0;
    double linearBeta = 0;
    double nonlinearBeta = 0;

    // Set up branches
    tree->Branch("mcBeta", &mcBeta, "mcBeta/D");
    tree->Branch("linearBeta", &linearBeta, "linearBeta/D");
    tree->Branch("nonlinearBeta", &nonlinearBeta, "nonlinearBeta/D");

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

        // Get beta values
        mcBeta = particle.mcBeta;
        linearBeta = particle.betaLinear;
        nonlinearBeta = 1 / BetaFitter::reconstructBeta(&particle, propagator);

        // Skip invalid reconstructions
        if (nonlinearBeta <= 0 || linearBeta <= 0)
            continue;

        // Fill tree
        tree->Fill();
    }

    // Write and close
    tree->Write();
    outFile->Close();
    delete outFile;

    return 0;
}