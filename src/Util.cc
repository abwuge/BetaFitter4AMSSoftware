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

bool Util::saveEnergyLoss(const std::string &inputFile, const std::string &outputFile)
{
    TFile *fileIn = TFile::Open(inputFile.c_str(), "READ");
    if (!fileIn || fileIn->IsZombie())
    {
        std::cerr << "Error: Could not open file " << inputFile << std::endl;
        return false;
    }

    TTree *treeIn = (TTree *)fileIn->Get("amstreea");
    if (!treeIn)
    {
        std::cerr << "Error: Could not find amstreea tree in " << inputFile << std::endl;
        fileIn->Close();
        return false;
    }

    if (!treeIn->GetBranch("mpar"))
    {
        std::cerr << "Error: Could not find Monte Carlo information in " << inputFile << std::endl;
        fileIn->Close();
        return false;
    }

    int mpar = 0;
    float mch = 0.0f;
    float mevmom1[21] = {0};
    float tof_edep[4] = {0};

    treeIn->SetBranchAddress("mpar", &mpar);
    treeIn->SetBranchAddress("mch", &mch);
    treeIn->SetBranchAddress("mevmom1", mevmom1);
    treeIn->SetBranchAddress("tof_edep", tof_edep);

    TFile *fileOut = new TFile(outputFile.c_str(), "RECREATE");
    if (!fileOut || fileOut->IsZombie())
    {
        std::cerr << "Error: Could not create output file " << outputFile << std::endl;
        return false;
    }

    TTree *treeOut = new TTree("energyLoss", "Energy Loss Information");

    float energyDepositedS1S2 = 0.0;  // energy deposited from before S1 to after S2
    float energyDepositedTotal = 0.0; // total energy deposited
    float energyLoss_S1S2_ = 0.0;     // energy loss from before S1 to after S2
    float energyLoss_S1S4_ = 0.0;     // energy loss from before S1 to after S4
    float energyLossScaleS1S2 = 0.0;  // energy loss scale factor from before S1 to after S2
    float energyLossScaleTotal = 0.0; // energy loss scale factor from before S1 to after S4
    float energyLossS2__S3 = 0.0;     // energy loss from after S2 to before S3
    float energyLossS2S3_Total = 0.0; // energy loss from after S2 to before S3, normalized to total energy loss
    float mcBeta = 0.0;               // Monte Carlo beta

    treeOut->Branch("energyDepositedS1S2", &energyDepositedS1S2, "energyDepositedS1S2/F");
    treeOut->Branch("energyDepositedTotal", &energyDepositedTotal, "energyDepositedTotal/F");
    treeOut->Branch("energyLoss_S1S2_", &energyLoss_S1S2_, "energyLoss_S1S2_/F");
    treeOut->Branch("energyLoss_S1S4_", &energyLoss_S1S4_, "energyLoss_S1S4_/F");
    treeOut->Branch("energyLossScaleS1S2", &energyLossScaleS1S2, "energyLossScaleS1S2/F");
    treeOut->Branch("energyLossScaleTotal", &energyLossScaleTotal, "energyLossScaleTotal/F");
    treeOut->Branch("energyLossS2__S3", &energyLossS2__S3, "energyLossS2__S3/F");
    treeOut->Branch("energyLossS2S3_Total", &energyLossS2S3_Total, "energyLossS2S3_Total/F");
    treeOut->Branch("mcBeta", &mcBeta, "mcBeta/F");

    // Read all entries
    Long64_t nEntries = treeIn->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        treeIn->GetEntry(i);

        if (mch != 8)
            continue;

        if (mevmom1[4] == -1000 || mevmom1[6] == -1000 || mevmom1[7] == -1000 || mevmom1[15] == -1000 || mevmom1[14] == -1000 || mevmom1[17] == -1000)
            continue;

        double mass = getMass(mpar, mch);

        /**
         * Kinetic energy index - z-position:
         * Index     Z (cm)     MC Index
         *   0       65.97          4
         *   -       65.20          -        (TOF S1)
         *   1       62.87          6
         *   -       62.10          -        (TOF S2)
         *   2       53.06          7
         *   3      -61.33         15
         *   -      -62.10          -        (TOF S3)
         *   4      -63.27         14
         *   -      -65.20          -        (TOF S4)
         *   5      -69.98         17
         */
        double kineticEnergy[6]{};
        kineticEnergy[0] = sqrt(mevmom1[4] * mevmom1[4] + mass * mass) - mass;
        kineticEnergy[1] = sqrt(mevmom1[6] * mevmom1[6] + mass * mass) - mass;
        kineticEnergy[2] = sqrt(mevmom1[7] * mevmom1[7] + mass * mass) - mass;
        kineticEnergy[3] = sqrt(mevmom1[15] * mevmom1[15] + mass * mass) - mass;
        kineticEnergy[4] = sqrt(mevmom1[14] * mevmom1[14] + mass * mass) - mass;
        kineticEnergy[5] = sqrt(mevmom1[17] * mevmom1[17] + mass * mass) - mass;

        energyDepositedS1S2 = (tof_edep[0] + tof_edep[1]) * 1e-3;
        energyDepositedTotal = (tof_edep[0] + tof_edep[1] + tof_edep[2] + tof_edep[3]) * 1e-3;
        energyLoss_S1S2_ = kineticEnergy[0] - kineticEnergy[2];
        energyLoss_S1S4_ = kineticEnergy[0] - kineticEnergy[5];
        energyLossScaleS1S2 = energyLoss_S1S2_ / energyDepositedS1S2;
        energyLossScaleTotal = energyLoss_S1S4_ / energyDepositedTotal;
        energyLossS2__S3 = kineticEnergy[2] - kineticEnergy[3];
        energyLossS2S3_Total = energyLossS2__S3 / energyLoss_S1S4_;
        mcBeta = mevmom1[4] / sqrt(mevmom1[4] * mevmom1[4] + mass * mass);

        treeOut->Fill();
    }

    treeOut->Write();
    fileOut->Close();

    fileIn->Close();
    return true;
}