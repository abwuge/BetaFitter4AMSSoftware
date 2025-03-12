#include "Util.hh"
#include "amschain.h"
#include "ParticlePropagator.hh"
#include <fstream>
#include <TFile.h>
#include <TTree.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TCanvas.h>
#include <TView.h>
#include <TText.h>
#include "TDatabasePDG.h"

// Implementation of getMass function
float Util::getMass(int pdgId, double charge)
{
    switch (pdgId)
    {
    case 69: // O16
        return 14.899169;
    default:
        return 2 * charge * 0.9314941; // Suppose its' number of neutron equals to number of proton
    }
}

std::vector<ParticleData> Util::loadParticleData(const std::string &inputFile)
{
    std::vector<ParticleData> particles;

    // Open ROOT file
    TFile *file = TFile::Open(inputFile.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << inputFile << std::endl;
        return particles;
    }

    // Get amstreea tree
    TTree *tree = (TTree *)file->Get("amstreea");
    if (!tree)
    {
        std::cerr << "Error: Could not find amstreea tree in " << inputFile << std::endl;
        file->Close();
        return particles;
    }

    // Variables to read from tree
    bool isMC = false;
    int mpar = 0;
    float mch = 0.0f;
    float mevcoo1[21][3] = {0};
    float mevdir1[21][3] = {0};
    float mevmom1[21] = {0};
    float tof_betah = 0.0f;
    float tof_tl[4] = {0};
    float tof_pos[4][3] = {0};
    float tof_dir[4][3] = {0};
    float tof_edep[4] = {0};
    float tk_pos[9][3] = {0};
    float tk_q[2] = {0};
    float tk_qin[2][3] = {0};

    // Set branch addresses
    if (tree->GetBranch("mpar"))
    {
        isMC = true;
        tree->SetBranchAddress("mpar", &mpar);
        tree->SetBranchAddress("mch", &mch);
        tree->SetBranchAddress("mevcoo1", mevcoo1);
        tree->SetBranchAddress("mevdir1", mevdir1);
        tree->SetBranchAddress("mevmom1", mevmom1);
    }
    tree->SetBranchAddress("tof_betah", &tof_betah);
    tree->SetBranchAddress("tof_tl", tof_tl);
    tree->SetBranchAddress("tof_pos", tof_pos);
    tree->SetBranchAddress("tof_dir", tof_dir);
    tree->SetBranchAddress("tof_edep", tof_edep);
    tree->SetBranchAddress("tk_pos", tk_pos);
    tree->SetBranchAddress("tk_q", tk_q);
    tree->SetBranchAddress("tk_qin", tk_qin);

    // Read all entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);

        ParticleData data;
        if (isMC)
        {
            data.isMC = true;
            data.mcGeantId = mpar;
            data.mcCharge = mch;
            float mmass = getMass(mpar, mch);
            data.mcMass = mmass;

            // Index 4 --> z = 65.95 cm --> utof(65.2 cm)
            data.mcInitCoo[0] = mevcoo1[4][0];
            data.mcInitCoo[1] = mevcoo1[4][1];
            data.mcInitCoo[2] = mevcoo1[4][2];

            data.mcInitDir[0] = mevdir1[4][0];
            data.mcInitDir[1] = mevdir1[4][1];
            data.mcInitDir[2] = mevdir1[4][2];

            float mmom = mevmom1[4];
            data.mcMomentum = mmom;
            data.mcBeta = mmom / sqrt(mmom * mmom + mmass * mmass);
        }

        data.initCoo[0] = tof_pos[0][0];
        data.initCoo[1] = tof_pos[0][1];
        data.initCoo[2] = tof_pos[0][2];

        data.initDir[0] = tof_dir[0][0];
        data.initDir[1] = tof_dir[0][1];
        data.initDir[2] = tof_dir[0][2];

        float minTime = *std::min_element(tof_tl, tof_tl + 4);
        for (int j = 0; j < ParticleData::TOF_MAX_HITS; ++j)
        {
            // TODO: this position might not be correct
            data.TOF_hitZ[j] = tof_pos[j][2];
            data.TOF_hitTime[j] = tof_tl[j] == -1 ? -1 : tof_tl[j] - minTime;
            // TODO: SO FAR it's constant, but it might not be correct
            data.TOF_hitTimeError[j] = 0.1544809;
            data.TOF_hitEdep[j] = tof_edep[j] * 1e-3;
        }

        for (int j = 0; j < ParticleData::TRACKER_MAX_HITS; ++j)
        {
            data.TRACKER_hitX[j] = tk_pos[j][0];
            data.TRACKER_hitY[j] = tk_pos[j][1];
            data.TRACKER_hitZ[j] = tk_pos[j][2];
            // TODO: Though here it's constant, but it's NOT totally correct
            data.TRACKER_hitError[j] = 6.3e-4;
        }

        data.charge = (int)((tk_qin[0][2] < 2.5 ? tk_q[1] : tk_qin[0][2]) + 0.5);
        data.mass = 2 * data.charge * 0.9314941; // Suppose its' number of neutron equals to number of proton
        data.betaLinear = tof_betah;
        data.momentum = data.mass * data.betaLinear / sqrt(1 - data.betaLinear * data.betaLinear);

        // TODO: This is a temporary method for filtering out fragmental particles
        if (data.charge < 6)
            continue;

        particles.push_back(data);
    }

    file->Close();
    return particles;
}

bool Util::saveMagneticField(const std::string &outputFile)
{
    // Create output ROOT file
    TFile *outFile = new TFile(outputFile.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie())
    {
        std::cerr << "Error: Could not create output file " << outputFile << std::endl;
        return false;
    }

    // Define the dimensions and ranges for magnetic field sampling
    const int nx = 101;
    const int ny = 101;
    const int nz = 101;
    const double xmin = -100, xmax = 100; // cm
    const double ymin = -100, ymax = 100; // cm
    const double zmin = -100, zmax = 100; // cm

    // Create TTree to store magnetic field data
    TTree *magTree = new TTree("magfield", "Magnetic Field Data");
    double x, y, z;
    double bx, by, bz;
    double b_magnitude;

    // Set up branches
    magTree->Branch("x", &x, "x/D");
    magTree->Branch("y", &y, "y/D");
    magTree->Branch("z", &z, "z/D");
    magTree->Branch("bx", &bx, "bx/D");
    magTree->Branch("by", &by, "by/D");
    magTree->Branch("bz", &bz, "bz/D");
    magTree->Branch("b_magnitude", &b_magnitude, "b_magnitude/D");

    // Sample magnetic field at grid points
    double dx = (xmax - xmin) / (nx - 1);
    double dy = (ymax - ymin) / (ny - 1);
    double dz = (zmax - zmin) / (nz - 1);
    double bf[3];

    for (int ix = 0; ix < nx; ix++)
    {
        x = xmin + ix * dx;
        for (int iy = 0; iy < ny; iy++)
        {
            y = ymin + iy * dy;
            for (int iz = 0; iz < nz; iz++)
            {
                z = zmin + iz * dz;

                // Get magnetic field at this point using AMS software
                TrFit::GuFld(x, y, z, bf);

                bx = bf[0];
                by = bf[1];
                bz = bf[2];
                b_magnitude = sqrt(bx * bx + by * by + bz * bz);

                magTree->Fill();
            }
        }
    }

    // Write and close
    magTree->Write();
    outFile->Close();
    delete outFile;

    return true;
}

bool Util::calculateTOFEnergyLoss(const ParticleData &particle, double beta, double energyLoss[4])
{
    // not implemented

    // // Calculate momentum and rigidity from beta
    // double gamma = 1.0 / std::sqrt(1.0 - beta * beta);
    // double momentum = gamma * beta * particle.mass;
    // double rigidity = momentum / particle.charge;

    // // Create particle propagator
    // ParticlePropagator propagator(particle);
    // propagator.resetPropagator(beta);

    // // Calculate energy loss for each TOF layer
    // for (int i = 0; i < ParticleData::TOF_MAX_HITS; ++i)
    // {
    //     // Calculate energy loss
    //     double layer_length = propagator.Propagate(particle.TOF_hitZ[i]);
    //     if (layer_length < 0)
    //         return false;

    //     double layer_time = layer_length / (beta * SPEED_OF_LIGHT);
    //     double layer_energy_loss = particle.mass * (1.0 / beta - 1.0 / gamma) * SPEED_OF_LIGHT * layer_time;
    //     energyLoss[i] = layer_energy_loss;
    // }

    return false;
}