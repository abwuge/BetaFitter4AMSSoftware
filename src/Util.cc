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

// Implementation of getMassFromPDG function
float Util::getMassFromPDG(int pdgId, double charge)
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
    float mmom = 0.0f, mch = 0.0f, mmass = 0.0f;
    float mevcoo1[21][3] = {0};
    float mevdir1[21][3] = {0};
    float mevmom1[21] = {0};
    float tof_betah = 0.0f;
    float tof_tl[4] = {0};
    float tof_pos[4][3] = {0};
    float tof_edep[4] = {0};
    float tk_pos[9][3] = {0};
    float tk_dir[9][3] = {0};
    float tk_q[2] = {0};
    float tk_qin[2][3] = {0};

    // Set branch addresses
    if (tree->GetBranch("mpar"))
    {
        isMC = true;
        tree->SetBranchAddress("mpar", &mpar);
        tree->SetBranchAddress("mmom", &mmom);
        tree->SetBranchAddress("mch", &mch);
        tree->SetBranchAddress("mevcoo1", mevcoo1);
        tree->SetBranchAddress("mevdir1", mevdir1);
        tree->SetBranchAddress("mevmom1", mevmom1);
    }
    tree->SetBranchAddress("tof_betah", &tof_betah);
    tree->SetBranchAddress("tof_tl", tof_tl);
    tree->SetBranchAddress("tof_pos", tof_pos);
    tree->SetBranchAddress("tof_edep", tof_edep);
    tree->SetBranchAddress("tk_pos", tk_pos);
    tree->SetBranchAddress("tk_dir", tk_dir);
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
            data.mcPdgId = mpar;
            data.mcCharge = mch;
            data.mcMomentum = mmom;
            mmass = getMassFromPDG(mpar, mch);
            data.mcMass = mmass;
            data.mcBeta = mmom / sqrt(mmom * mmom + mmass * mmass);

            data.mcCoo[0] = mevcoo1[0][0];
            data.mcCoo[1] = mevcoo1[0][1];
            data.mcCoo[2] = mevcoo1[0][2];

            data.mcDir[0] = mevdir1[0][0];
            data.mcDir[1] = mevdir1[0][1];
            data.mcDir[2] = mevdir1[0][2];
        }

        for (int j = 0; j < ParticleData::TOF_MAX_HITS; ++j)
        {
            // TODO: this position might not be correct
            data.TOF_hitZ[j] = tof_pos[j][2];
            data.TOF_hitTime[j] = tof_tl[j];
            // TODO: SO FAR it's constant, but it might not be correct
            data.TOF_hitTimeError[j] = 0.1544809;
            data.TOF_hitEdep[j] = tof_edep[j] * 1e-3;
        }

        data.TRACKER_dir[0] = tk_dir[0][0];
        data.TRACKER_dir[1] = tk_dir[0][1];
        data.TRACKER_dir[2] = tk_dir[0][2];

        for (int j = 0; j < ParticleData::TRACKER_MAX_HITS; ++j)
        {
            data.TRACKER_hitX[j] = tk_pos[j][0];
            data.TRACKER_hitY[j] = tk_pos[j][1];
            data.TRACKER_hitZ[j] = tk_pos[j][2];
            // TODO: Though here it's constant, but it's NOT totally correct
            data.TRACKER_hitError[j] = 6.3e-4;
        }

        data.charge = (tk_qin[0][2] < 2.5 ? tk_q[1] : tk_qin[0][2]);
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