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

std::vector<ParticleData> Util::loadParticleData(const std::string &inputFile)
{
    std::vector<ParticleData> particles;

    // Open ROOT file
    TFile *file = TFile::Open(inputFile.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        return particles;
    }

    // Get amstreea tree
    TTree *tree = (TTree *)file->Get("amstreea");
    if (!tree)
    {
        file->Close();
        return particles;
    }

    // Variables to read from tree
    int mpar;
    float mmom, mch, mmass;
    float mevcoo1[21][3];
    float mevdir1[21][3];
    float mevmom1[21];
    float tof_betah, tof_betahmc;
    float tof_tl[4];
    float tof_etl[4];

    // Set branch addresses
    tree->SetBranchAddress("mpar", &mpar);
    tree->SetBranchAddress("mmom", &mmom);
    tree->SetBranchAddress("mch", &mch);
    tree->SetBranchAddress("mmass", &mmass);
    tree->SetBranchAddress("mevcoo1", mevcoo1);
    tree->SetBranchAddress("mevdir1", mevdir1);
    tree->SetBranchAddress("mevmom1", mevmom1);
    tree->SetBranchAddress("tof_betah", &tof_betah);
    tree->SetBranchAddress("tof_betahmc", &tof_betahmc);
    tree->SetBranchAddress("tof_tl", tof_tl);
    tree->SetBranchAddress("tof_etl", tof_etl);

    // Read all entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);
        
        ParticleData data;
        data.isMC = true;
        data.mcMomentum = mmom;
        data.mcCharge = mch;
        data.mcMass = mmass;
        data.mcBeta = tof_betahmc;
        data.beta = tof_betah;
        
        data.mcCoo[0] = mevcoo1[0][0];
        data.mcCoo[1] = mevcoo1[0][1];
        data.mcCoo[2] = mevcoo1[0][2];
        
        data.mcDir[0] = mevdir1[0][0];
        data.mcDir[1] = mevdir1[0][1];
        data.mcDir[2] = mevdir1[0][2];

        for (int j = 0; j < 4; ++j) {
            data.hitTime[j] = tof_tl[j];
            data.hitTimeError[j] = tof_etl[j];
        }
        
        particles.push_back(data);
    }

    file->Close();
    return particles;
}

bool Util::drawTrajectory(const ParticleData &particle, const std::string &outputPath, int nPoints)
{
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
}