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

        if (data.charge < 6)
            continue;

        particles.push_back(data);
    }

    file->Close();
    return particles;
}

bool Util::drawTrajectory(const ParticleData &particle, const std::string &outputPath, int nPoints)
{
#if false
    try
    {
        // Create canvas and 3D view
        TCanvas *c = new TCanvas("c", "Particle Trajectory", 800, 600);
        c->SetFillColor(kWhite);

        // Create 3D view
        TView *view = TView::CreateView(1);
        view->SetRange(-50, -50, -66, 50, 50, 66);
        view->ShowAxis();

        // Create TOF hit points
        TPolyMarker3D *tofHits = new TPolyMarker3D(4);
        for (int i = 0; i < 4; ++i)
        {
            tofHits->SetPoint(i, particle.hitX[i], particle.hitY[i], particle.hitZ[i]);
        }
        tofHits->SetMarkerStyle(20);
        tofHits->SetMarkerColor(kRed);
        tofHits->SetMarkerSize(1.2);
        tofHits->Draw();

        // Use MC truth position and direction if available
        AMSPoint initPos;
        AMSDir initDir;

        if (particle.isMC)
        {
            // Use MC truth position and direction
            initPos = AMSPoint(particle.mcCoo[0], particle.mcCoo[1], particle.mcCoo[2]);
            // initPos = AMSPoint(particle.hitX[0], particle.hitY[0], particle.hitZ[0]);
            initDir = AMSDir(particle.mcDir[0], particle.mcDir[1], particle.mcDir[2]);
            // initDir = AMSDir(particle.mcDir[0], particle.mcDir[1], particle.mcDir[2]);
            initDir.SetTheta(particle.Theta);
            initDir.SetPhi(particle.Phi); 
        }
        else
        {
            // Fallback to reconstructed values
            initPos = AMSPoint(0, 0, 100); // Start from above top TOF
            double px = particle.momentum * sin(particle.Theta) * cos(particle.Phi);
            double py = particle.momentum * sin(particle.Theta) * sin(particle.Phi);
            double pz = particle.momentum * cos(particle.Theta);
            initDir = AMSDir(px, py, pz);
        }

        ParticlePropagator prop(initPos, initDir, particle.momentum,
                                particle.mass, particle.charge);

        // Create trajectory points
        TPolyLine3D *track = new TPolyLine3D(nPoints);
        double zStep = (ParticlePropagator::TOF_Z[3] - ParticlePropagator::TOF_Z[0]) / (nPoints - 1);

        for (int i = 0; i < nPoints; ++i)
        {
            double z = ParticlePropagator::TOF_Z[0] + i * zStep;
            double beta = prop.PropagateToZ(z);
            if (beta < 0)
                continue;

            AMSPoint pos = prop.GetP0();
            track->SetPoint(i, pos.x(), pos.y(), pos.z());
        }

        track->SetLineColor(kBlue);
        track->SetLineWidth(2);
        track->Draw();

        // Add labels for TOF layers
        for (int i = 0; i < 4; ++i)
        {
            TText *label = new TText(particle.hitX[i] + 10, particle.hitY[i] + 10,
                                     Form("TOF %d", i));
            label->SetTextSize(0.02);
            label->Draw();
        }

        c->Update();
        c->Print(outputPath.c_str(), "pdf");

        delete c; // This will delete all drawn objects
        return true;
    }
    catch (...)
    {
        return false;
    }
#endif
    return false;
}