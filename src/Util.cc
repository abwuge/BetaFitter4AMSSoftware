#include "Util.hh"
#include "amschain.h"
#include <fstream>
#include <TFile.h>
#include <TTree.h>

bool Util::addInputFile(std::string &inputFile, AMSChain &chain)
{
    try
    {
        if (inputFile.substr(inputFile.find_last_of(".") + 1) == "root")
        {
            chain.Add(inputFile.c_str());
            return true;
        }

        std::ifstream file(inputFile);
        if (!file.is_open())
            return false;

        std::string line;
        while (std::getline(file, line))
        {
            if (!line.empty())
            {
                chain.Add(line.c_str());
            }
        }

        return true;
    }
    catch (...)
    {
        return false;
    }
}

std::vector<ParticleData> Util::loadParticleData(const std::string &inputFile)
{
    std::vector<ParticleData> particles;

    // Open ROOT file
    TFile *file = TFile::Open(inputFile.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        return particles;
    }

    // Get tree
    TTree *tree = (TTree *)file->Get("particle");
    if (!tree)
    {
        file->Close();
        return particles;
    }

    // Pre-allocate space for particles based on the number of entries
    Long64_t nEntries = tree->GetEntries();
    particles.reserve(nEntries);

    // Create a particle data object and set branch addresses
    ParticleData data;

    // Set branch addresses for particle properties
    tree->SetBranchAddress("mass", &data.mass);
    tree->SetBranchAddress("charge", &data.charge);
    tree->SetBranchAddress("momentum", &data.momentum);
    tree->SetBranchAddress("beta", &data.beta);
    tree->SetBranchAddress("Theta", &data.Theta);
    tree->SetBranchAddress("Phi", &data.Phi);

    // Set branch addresses for hit information
    tree->SetBranchAddress("hitX", data.hitX);
    tree->SetBranchAddress("hitY", data.hitY);
    tree->SetBranchAddress("hitZ", data.hitZ);
    tree->SetBranchAddress("hitTime", data.hitTime);
    tree->SetBranchAddress("hitTimeError", data.hitTimeError);

    // Set branch addresses for MC truth information
    tree->SetBranchAddress("mcBeta", &data.mcBeta);
    tree->SetBranchAddress("mcMomentum", &data.mcMomentum);
    tree->SetBranchAddress("mcMass", &data.mcMass);
    tree->SetBranchAddress("mcCharge", &data.mcCharge);
    tree->SetBranchAddress("mcPdgId", &data.mcPdgId);
    tree->SetBranchAddress("isMC", &data.isMC);

    // Read all entries
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);
        particles.push_back(data);
    }

    file->Close();
    return particles;
}