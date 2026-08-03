#include "Util.hh"
#include "GlobalZetaFitter.hh"
#include "amschain.h"
#include "ParticlePropagator.hh"
#include <fstream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TCanvas.h>
#include <TView.h>
#include <TText.h>
#include <chrono>

namespace
{
bool addGlobalZetaInput(TChain &chain, const std::string &inputPath, int &filesAdded)
{
    filesAdded = 0;
    const bool isList = inputPath.size() >= 5 &&
                        inputPath.substr(inputPath.size() - 5) == ".list";
    if (!isList)
    {
        filesAdded = chain.Add(inputPath.c_str());
        return filesAdded > 0;
    }

    std::ifstream input(inputPath.c_str());
    if (!input)
        return false;
    std::string path;
    while (std::getline(input, path))
    {
        const std::string::size_type first = path.find_first_not_of(" \t\r");
        if (first == std::string::npos || path[first] == '#')
            continue;
        const std::string::size_type last = path.find_last_not_of(" \t\r");
        filesAdded += chain.Add(path.substr(first, last - first + 1).c_str());
    }
    return filesAdded > 0;
}

struct MomentumPoint
{
    double z;
    double momentum;
};

double momentumAt(const std::vector<MomentumPoint> &points, double z)
{
    if (points.empty() || z > points.front().z || z < points.back().z)
        return std::numeric_limits<double>::quiet_NaN();
    for (size_t index = 0; index + 1 < points.size(); ++index)
    {
        if (points[index].z + 1e-9 < z || z < points[index + 1].z - 1e-9)
            continue;
        const double deltaZ = points[index].z - points[index + 1].z;
        if (std::abs(deltaZ) < 1e-10)
            return 0.5 * (points[index].momentum + points[index + 1].momentum);
        const double fraction = (points[index].z - z) / deltaZ;
        return points[index].momentum * (1 - fraction) +
               points[index + 1].momentum * fraction;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double averageInverseBeta(const std::vector<MomentumPoint> &points,
                          double zFirst, double zSecond, double mass)
{
    if (zFirst < zSecond)
        std::swap(zFirst, zSecond);
    if (points.empty() || zFirst > points.front().z || zSecond < points.back().z ||
        !(zFirst > zSecond))
        return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> positions;
    positions.push_back(zFirst);
    positions.push_back(zSecond);
    for (const MomentumPoint &point : points)
        if (point.z < zFirst && point.z > zSecond)
            positions.push_back(point.z);
    std::sort(positions.begin(), positions.end(), std::greater<double>());
    positions.erase(std::unique(positions.begin(), positions.end(),
                                [](double left, double right)
                                { return std::abs(left - right) < 1e-6; }),
                    positions.end());

    const auto inverseBeta = [mass](double momentum)
    {
        return std::sqrt(1.0 + (mass / momentum) * (mass / momentum));
    };
    double integral = 0;
    for (size_t index = 0; index + 1 < positions.size(); ++index)
    {
        const double first = positions[index];
        const double second = positions[index + 1];
        const double middle = 0.5 * (first + second);
        const double momentumFirst = momentumAt(points, first);
        const double momentumSecond = momentumAt(points, second);
        const double momentumMiddle = momentumAt(points, middle);
        if (!(momentumFirst > 0) || !(momentumSecond > 0) || !(momentumMiddle > 0))
            return std::numeric_limits<double>::quiet_NaN();
        integral += (first - second) *
                    (inverseBeta(momentumFirst) + 4 * inverseBeta(momentumMiddle) +
                     inverseBeta(momentumSecond)) /
                    6.0;
    }
    return integral / (zFirst - zSecond);
}

bool buildCheckpointTruthTimes(const float momentum[21], const float coordinate[21][3],
                               const float tofPosition[4][3], const std::array<float, 4> &pathLength,
                               double mass, std::array<float, 4> &truthTime)
{
    std::vector<MomentumPoint> points;
    for (int index = 0; index < 21; ++index)
        if (momentum[index] > 0 && std::isfinite(momentum[index]) &&
            std::isfinite(coordinate[index][2]))
            points.push_back(MomentumPoint{coordinate[index][2], momentum[index]});
    std::sort(points.begin(), points.end(),
              [](const MomentumPoint &left, const MomentumPoint &right)
              { return left.z > right.z; });
    std::vector<MomentumPoint> uniquePoints;
    for (const MomentumPoint &point : points)
    {
        if (uniquePoints.empty() || std::abs(point.z - uniquePoints.back().z) > 1e-5)
            uniquePoints.push_back(point);
        else
            uniquePoints.back().momentum = 0.5 * (uniquePoints.back().momentum + point.momentum);
    }

    truthTime.fill(0);
    for (int station = 1; station < 4; ++station)
    {
        const double inverseBeta = averageInverseBeta(
            uniquePoints, tofPosition[station - 1][2], tofPosition[station][2], mass);
        const double segmentLength = pathLength[station] - pathLength[station - 1];
        if (!std::isfinite(inverseBeta) || !(segmentLength > 0 && segmentLength < 300))
            return false;
        truthTime[station] = static_cast<float>(truthTime[station - 1] +
                                                segmentLength / 29.9792458 * inverseBeta);
    }
    return true;
}
}

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

std::vector<ParticleData> Util::loadParticleData(const std::string &inputFile,
                                                 BetaReferencePoint referencePoint,
                                                 bool requireTrackerEnergyLoss)
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
    bool isMC{};
    int mpar{};
    float mch{};
    float mevcoo1[21][3]{};
    float mevdir1[21][3]{};
    float mevmom1[21]{};
    float tof_betah{};
    float tof_tl[4]{};
    float tof_pos[4][3]{};
    float tof_dir[4][3]{};
    float tof_edep[4]{};
    float tk_pos[9][3]{};
    float tk_edep[9]{};
    float tk_q[2]{};
    float tk_qin[2][3]{};
    float tk_rigidity1[3][3][7]{};
    float tof_leng[4]{};

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
    if (requireTrackerEnergyLoss && !tree->GetBranch("tk_edep"))
    {
        std::cerr << "Error: Missing required tk_edep branch in " << inputFile << std::endl;
        file->Close();
        return particles;
    }
    if (tree->GetBranch("tk_edep"))
        tree->SetBranchAddress("tk_edep", tk_edep);
    tree->SetBranchAddress("tk_q", tk_q);
    tree->SetBranchAddress("tk_qin", tk_qin);
    tree->SetBranchAddress("tk_rigidity1", tk_rigidity1);
    tree->SetBranchAddress("tof_leng", tof_leng);

    // Read all entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        tree->GetEntry(i);

        const int mcReferenceIndex = referencePoint == BetaReferencePoint::AMSCenter ? 10 : 4;
        if (mevmom1[mcReferenceIndex] == -1000)
            continue;

        ParticleData data;
        if (isMC)
        {
            data.isMC = true;
            data.mcGeantId = mpar;
            data.mcCharge = mch;
            float mmass = getMass(mpar, mch);
            data.mcMass = mmass;

            data.mcInitCoo[0] = mevcoo1[mcReferenceIndex][0];
            data.mcInitCoo[1] = mevcoo1[mcReferenceIndex][1];
            data.mcInitCoo[2] = mevcoo1[mcReferenceIndex][2];

            data.mcInitDir[0] = mevdir1[mcReferenceIndex][0];
            data.mcInitDir[1] = mevdir1[mcReferenceIndex][1];
            data.mcInitDir[2] = mevdir1[mcReferenceIndex][2];

            float mmom = mevmom1[mcReferenceIndex];
            data.mcMomentum = mmom;
            data.mcBeta = mmom / sqrt(mmom * mmom + mmass * mmass);
        }

        data.initCoo[0] = tof_pos[0][0];
        data.initCoo[1] = tof_pos[0][1];
        data.initCoo[2] = tof_pos[0][2];

        data.initDir[0] = tof_dir[0][0];
        data.initDir[1] = tof_dir[0][1];
        data.initDir[2] = tof_dir[0][2];

        for (int j = 0; j < ParticleData::TOF_MAX_HITS; ++j)
        {
            // TODO: verify that tof_pos is the hit position used for path calculation
            data.TOF_hitX[j] = tof_pos[j][0];
            data.TOF_hitY[j] = tof_pos[j][1];
            data.TOF_hitZ[j] = tof_pos[j][2];
            data.TOF_hitTime[j] = tof_tl[j];
            // TODO: SO FAR it's constant, but it might not be correct
            data.TOF_hitTimeError[j] = 0.1544809;
            data.TOF_hitEdep[j] = tof_edep[j];
            data.TOF_length[j] = tof_leng[j];
        }

        for (int j = 0; j < ParticleData::TRACKER_MAX_HITS; ++j)
        {
            data.TRACKER_hitX[j] = tk_pos[j][0];
            data.TRACKER_hitY[j] = tk_pos[j][1];
            data.TRACKER_hitZ[j] = tk_pos[j][2];
            data.TRACKER_hitEdep[j] = tk_edep[j];
            // TODO: Though here it's constant, but it's NOT totally correct
            data.TRACKER_hitError[j] = 6.3e-4;
        }

        data.charge = (int)((tk_qin[0][2] < 2.5 ? tk_q[1] : tk_qin[0][2]) + 0.5);
        data.innerRigidity = tk_rigidity1[1][2][1];
        data.mass = 2 * data.charge * 0.9314941; // Suppose its' number of neutron equals to number of proton
        data.betaLinear = tof_betah;
        float momentumRigidity = data.innerRigidity * data.charge;
        data.betaRigidity = momentumRigidity / TMath::Sqrt(momentumRigidity * momentumRigidity + data.mass * data.mass);
        data.momentum = data.mass * data.betaLinear / sqrt(1 - data.betaLinear * data.betaLinear);

        particles.push_back(data);
    }

    file->Close();
    return particles;
}

bool Util::saveBeta(const std::string &inputFile,
                    const std::string &outputFile,
                    double energyLossScale,
                    EnergyLossScaleMode energyLossScaleMode,
                    BetaReferencePoint referencePoint)
{
    // Load particle data from input file
    std::vector<ParticleData> particles =
        Util::loadParticleData(inputFile, referencePoint, true);
    if (particles.empty())
    {
        std::cerr << "Error: No particles loaded from input file" << std::endl;
        return false;
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
    double innerRigidity = 0;
    double betaRigidity = 0;
    double linearBeta = 0;
    double nonlinearBetaLegacy = 0;
    double nonlinearBeta = 0;
    int Z = 0;

    // Set up branches
    tree->Branch("mcBeta", &mcBeta, "mcBeta/D");
    tree->Branch("innerRigidity", &innerRigidity, "innerRigidity/D");
    tree->Branch("betaRigidity", &betaRigidity, "betaRigidity/D");
    tree->Branch("linearBeta", &linearBeta, "linearBeta/D");
    tree->Branch("nonlinearBetaLegacy", &nonlinearBetaLegacy, "nonlinearBetaLegacy/D");
    tree->Branch("nonlinearBeta", &nonlinearBeta, "nonlinearBeta/D");
    tree->Branch("Z", &Z, "Z/I");

    // Process each particle
    for (const auto &particle : particles)
    {
        bool trackerComplete = true;
        for (int layer = 1; layer <= 7; ++layer)
            trackerComplete = trackerComplete && particle.TRACKER_hitEdep[layer] > 0 &&
                              std::isfinite(particle.TRACKER_hitEdep[layer]) &&
                              std::isfinite(particle.TRACKER_hitX[layer]) &&
                              std::isfinite(particle.TRACKER_hitY[layer]) &&
                              std::isfinite(particle.TRACKER_hitZ[layer]);
        if (!trackerComplete || !(particle.TRACKER_hitZ[4] > 0) ||
            !(particle.TRACKER_hitZ[5] < 0))
            continue;

        float tofPosition[ParticleData::TOF_MAX_HITS][3];
        float trackerPosition[ParticleData::TRACKER_MAX_HITS][3];
        for (int station = 0; station < ParticleData::TOF_MAX_HITS; ++station)
        {
            tofPosition[station][0] = particle.TOF_hitX[station];
            tofPosition[station][1] = particle.TOF_hitY[station];
            tofPosition[station][2] = particle.TOF_hitZ[station];
        }
        for (int layer = 0; layer < ParticleData::TRACKER_MAX_HITS; ++layer)
        {
            trackerPosition[layer][0] = particle.TRACKER_hitX[layer];
            trackerPosition[layer][1] = particle.TRACKER_hitY[layer];
            trackerPosition[layer][2] = particle.TRACKER_hitZ[layer];
        }

        // Get beta values
        mcBeta = particle.mcBeta;
        innerRigidity = particle.innerRigidity;
        betaRigidity = particle.betaRigidity;
        linearBeta = particle.betaLinear;
        nonlinearBetaLegacy =
            BetaNL(
                BetaNLPars(
                    particle.betaLinear,
                    particle.mass,
                    particle.TOF_hitEdep,
                    particle.TOF_hitTime,
                    particle.TOF_hitTimeError,
                    particle.TOF_length),
                energyLossScale,
                energyLossScaleMode,
                referencePoint)
                .Beta();
        nonlinearBeta =
            BetaNL(
                BetaNLPars(
                    particle.betaLinear,
                    particle.mass,
                    particle.TOF_hitEdep,
                    particle.TOF_hitTime,
                    particle.TOF_hitTimeError,
                    particle.TOF_length,
                    tofPosition,
                    particle.TRACKER_hitEdep,
                    trackerPosition),
                energyLossScale,
                energyLossScaleMode,
                referencePoint)
                .Beta();
        Z = particle.charge;

        // Fill tree
        tree->Fill();
    }

    // Write and close
    tree->Write();
    outFile->Close();
    delete outFile;

    return true;
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

    std::cout << "Magnetic field data saved to " << outputFile << std::endl;

    return true;
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

    int tof_qs = 0; // Q Status (1111: all unoverlapped, 0000: all overlapped, left to right: S1, S2, S3, S4)
    int mpar = 0;
    float mch = 0.0f;
    float mevmom1[21]{};
    float tof_edep[4]{};
    float mevcoo1[21][3]{};
    float mevdir1[21][3]{};

    treeIn->SetBranchAddress("tof_qs", &tof_qs);
    treeIn->SetBranchAddress("mpar", &mpar);
    treeIn->SetBranchAddress("mch", &mch);
    treeIn->SetBranchAddress("mevmom1", mevmom1);
    treeIn->SetBranchAddress("tof_edep", tof_edep);
    treeIn->SetBranchAddress("mevcoo1", mevcoo1);
    treeIn->SetBranchAddress("mevdir1", mevdir1);

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
    float energyLoss_S1L3 = 0.0;      // energy loss from before S1 to L3
    float energyLoss_S1L4 = 0.0;      // energy loss from before S1 to L4
    float energyLoss_S1S4_ = 0.0;     // energy loss from before S1 to after S4
    float energyLossScaleS1S2 = 0.0;  // energy loss scale factor from before S1 to after S2
    float energyLossScaleTotal = 0.0; // energy loss scale factor from before S1 to after S4
    float energyLossS2__S3 = 0.0;     // energy loss from after S2 to before S3
    float energyLossS2S3_Total = 0.0; // energy loss from after S2 to before S3, normalized to total energy loss
    float mcBeta = 0.0;               // Monte Carlo beta
    float position[3]{};
    float direction[3]{};

    treeOut->Branch("energyDepositedS1S2", &energyDepositedS1S2, "energyDepositedS1S2/F");
    treeOut->Branch("energyDepositedTotal", &energyDepositedTotal, "energyDepositedTotal/F");
    treeOut->Branch("energyLoss_S1S2_", &energyLoss_S1S2_, "energyLoss_S1S2_/F");
    treeOut->Branch("energyLoss_S1L3", &energyLoss_S1L3, "energyLoss_S1L3/F");
    treeOut->Branch("energyLoss_S1L4", &energyLoss_S1L4, "energyLoss_S1L4/F");
    treeOut->Branch("energyLoss_S1S4_", &energyLoss_S1S4_, "energyLoss_S1S4_/F");
    treeOut->Branch("energyLossScaleS1S2", &energyLossScaleS1S2, "energyLossScaleS1S2/F");
    treeOut->Branch("energyLossScaleTotal", &energyLossScaleTotal, "energyLossScaleTotal/F");
    treeOut->Branch("energyLossS2__S3", &energyLossS2__S3, "energyLossS2__S3/F");
    treeOut->Branch("energyLossS2S3_Total", &energyLossS2S3_Total, "energyLossS2S3_Total/F");
    treeOut->Branch("mcBeta", &mcBeta, "mcBeta/F");
    treeOut->Branch("tof_qs", &tof_qs, "tof_qs/I");
    treeOut->Branch("position", position, "position[3]/F");
    treeOut->Branch("direction", direction, "direction[3]/F");

    // Read all entries
    Long64_t nEntries = treeIn->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        treeIn->GetEntry(i);

        if (mevmom1[4] == -1000 || mevmom1[6] == -1000 || mevmom1[7] == -1000 || mevmom1[15] == -1000 || mevmom1[14] == -1000 || mevmom1[17] == -1000)
            continue;

        double mass = getMass(mpar, mch);

        /**
         * Kinetic energy index - z-position:
         * Index     Z (cm)     MC Index
         *   -       64.43          3
         *   0       65.97          4
         *   -       65.20          -        (TOF S1)
         *   1       62.87          6
         *   -       62.10          -        (TOF S2)
         *   2       53.06          7        (TR  L2)
         *   3       29.23          8        (TR  L3)
         *   4       25.21          9        (TR  L4)
         *   5      -61.33         15
         *   -      -62.10          -        (TOF S3)
         *   6      -63.27         14
         *   -      -65.20          -        (TOF S4)
         *   7      -69.98         17
         */
        double kineticEnergy[8]{};
        kineticEnergy[0] = sqrt(mevmom1[4] * mevmom1[4] + mass * mass) - mass;
        kineticEnergy[1] = sqrt(mevmom1[6] * mevmom1[6] + mass * mass) - mass;
        kineticEnergy[2] = sqrt(mevmom1[7] * mevmom1[7] + mass * mass) - mass;
        kineticEnergy[3] = sqrt(mevmom1[8] * mevmom1[8] + mass * mass) - mass;
        kineticEnergy[4] = sqrt(mevmom1[9] * mevmom1[9] + mass * mass) - mass;
        kineticEnergy[5] = sqrt(mevmom1[15] * mevmom1[15] + mass * mass) - mass;
        kineticEnergy[6] = sqrt(mevmom1[14] * mevmom1[14] + mass * mass) - mass;
        kineticEnergy[7] = sqrt(mevmom1[17] * mevmom1[17] + mass * mass) - mass;

        energyDepositedS1S2 = (tof_edep[0] + tof_edep[1]) * 1e-3;
        energyDepositedTotal = (tof_edep[0] + tof_edep[1] + tof_edep[2] + tof_edep[3]) * 1e-3;
        energyLoss_S1S2_ = kineticEnergy[0] - kineticEnergy[2];
        energyLoss_S1L3 = kineticEnergy[0] - kineticEnergy[3];
        energyLoss_S1L4 = kineticEnergy[0] - kineticEnergy[4];
        energyLoss_S1S4_ = kineticEnergy[0] - kineticEnergy[7];
        energyLossScaleS1S2 = energyLoss_S1S2_ / energyDepositedS1S2;
        energyLossScaleTotal = energyLoss_S1S4_ / energyDepositedTotal;
        energyLossS2__S3 = kineticEnergy[2] - kineticEnergy[5];
        energyLossS2S3_Total = energyLossS2__S3 / energyLoss_S1S4_;
        mcBeta = mevmom1[4] / sqrt(mevmom1[4] * mevmom1[4] + mass * mass);

        for (int i = 0; i < 3; ++i)
        {
            position[i] = mevcoo1[4][i];
            direction[i] = mevdir1[4][i];
        }

        treeOut->Fill();
    }

    if (treeOut->GetEntries() == 0)
    {
        std::cerr << "Error: No valid data saved to " << outputFile << std::endl;
        fileOut->Close();
        return false;
    }
    std::cout << "Total entries: " << treeOut->GetEntries() << std::endl;
    treeOut->Write();
    std::cout << "Energy loss information saved to " << outputFile << std::endl;

    fileOut->Close();
    fileIn->Close();

    return true;
}

bool Util::saveEnergyLossScale(const std::string &inputFile,
                               const std::string &outputFile,
                               EnergyLossScaleMode energyLossScaleMode,
                               BetaReferencePoint referencePoint)
{
    // Load particle data from input file
    std::vector<ParticleData> particles = Util::loadParticleData(inputFile, referencePoint);
    if (particles.empty())
    {
        std::cerr << "Error: No particles loaded from input file" << std::endl;
        return false;
    }

    // Create output ROOT file
    TFile *outFile = new TFile(outputFile.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie())
    {
        std::cerr << "Error: Could not create output file" << std::endl;
        return 1;
    }

    // Create TTree
    TTree *tree = new TTree("scaleTree", "Energy Loss Scale Factor");

    float energyLossScale;
    float mcBeta;
    float position[3];
    float direction[3];

    tree->Branch("energyLossScale", &energyLossScale, "energyLossScale/F");
    tree->Branch("mcBeta", &mcBeta, "mcBeta/F");
    tree->Branch("position", position, "position[3]/F");
    tree->Branch("direction", direction, "direction[3]/F");

    for (const auto &particle : particles)
    {
        // Skip invalid particles
        if (!particle.isMC)
            continue;

        if (particle.mcInitCoo[0] == -1000)
            continue;

        // Get energy loss scale factor
        energyLossScale = BetaNL(
                              BetaNLPars(
                                  particle.betaLinear,
                                  particle.mass,
                                  particle.TOF_hitEdep,
                                  particle.TOF_hitTime,
                                  particle.TOF_hitTimeError,
                                  particle.TOF_length),
                              2,
                              energyLossScaleMode,
                              referencePoint)
                              .EnergyLossScale(particle.mcBeta);
        mcBeta = particle.mcBeta;
        position[0] = particle.initCoo[0];
        position[1] = particle.initCoo[1];
        position[2] = particle.initCoo[2];
        direction[0] = particle.initDir[0];
        direction[1] = particle.initDir[1];
        direction[2] = particle.initDir[2];

        // Fill tree
        tree->Fill();
    }

    // Write and close
    tree->Write();
    outFile->Close();
    delete outFile;

    return true;
}

bool Util::saveGlobalEnergyLossScale(const std::string &inputFile,
                                     const std::string &outputFile,
                                     double betaMax,
                                     double zetaMin,
                                     double zetaMax,
                                     EnergyLossScaleMode energyLossScaleMode,
                                     BetaReferencePoint referencePoint)
{
    if (!(betaMax > 0.1 && betaMax < 1.0))
    {
        std::cerr << "Error: Global fit beta maximum must be between 0.1 and 1" << std::endl;
        return false;
    }
    if (!(zetaMin < zetaMax))
    {
        std::cerr << "Error: Global fit zeta minimum must be below maximum" << std::endl;
        return false;
    }

    TChain chain("amstreea");
    int filesAdded = 0;
    if (!addGlobalZetaInput(chain, inputFile, filesAdded) || chain.GetEntries() <= 0)
    {
        std::cerr << "Error: Could not load global-fit input " << inputFile << std::endl;
        return false;
    }
    chain.LoadTree(0);
    const char *requiredBranches[] = {
        "mpar", "mch", "mevmom1", "mevcoo1", "tof_tl", "tof_edep",
        "tof_leng", "tof_pos", "tk_edep", "tk_pos"};
    for (const char *branch : requiredBranches)
        if (!chain.GetBranch(branch))
        {
            std::cerr << "Error: Missing MC branch " << branch << std::endl;
            return false;
        }

    int mpar = 0;
    float mch = 0;
    float mevmom1[21] = {};
    float mevcoo1[21][3] = {};
    float tofTime[4] = {};
    float tofEnergyDeposited[4] = {};
    float tofLength[4] = {};
    float tofPosition[4][3] = {};
    float trackerEnergyDeposited[9] = {};
    float trackerPosition[9][3] = {};
    chain.SetBranchStatus("*", 0);
    for (const char *branch : requiredBranches)
        chain.SetBranchStatus(branch, 1);
    chain.SetBranchAddress("mpar", &mpar);
    chain.SetBranchAddress("mch", &mch);
    chain.SetBranchAddress("mevmom1", mevmom1);
    chain.SetBranchAddress("mevcoo1", mevcoo1);
    chain.SetBranchAddress("tof_tl", tofTime);
    chain.SetBranchAddress("tof_edep", tofEnergyDeposited);
    chain.SetBranchAddress("tof_leng", tofLength);
    chain.SetBranchAddress("tof_pos", tofPosition);
    chain.SetBranchAddress("tk_edep", trackerEnergyDeposited);
    chain.SetBranchAddress("tk_pos", trackerPosition);

    const Long64_t inputEntries = chain.GetEntries();
    Long64_t validFourHitEntries = 0;
    Long64_t trackerCompleteEntries = 0;
    Long64_t betaSelectedEntries = 0;
    Long64_t truthIntegrableEntries = 0;
    std::vector<GlobalZetaEvent> stableEvents;
    stableEvents.reserve(std::max<Long64_t>(10000, inputEntries / 4));
    const std::vector<GlobalZetaEvent> emptyEvents;
    GlobalZetaFitter validator(emptyEvents, energyLossScaleMode, referencePoint);
    int fittedCharge = 0;
    bool mixedCharge = false;
    const int mcReferenceIndex = referencePoint == BetaReferencePoint::AMSCenter ? 10 : 4;
    for (Long64_t entry = 0; entry < inputEntries; ++entry)
    {
        chain.GetEntry(entry);
        const double mass = getMass(mpar, mch);
        const double momentum = mevmom1[mcReferenceIndex];
        if (!(mass > 0) || !(momentum > 0) || !std::isfinite(momentum))
            continue;

        GlobalZetaEvent event;
        bool valid = true;
        for (int station = 0; station < 4; ++station)
        {
            if (tofTime[station] == -1 || !std::isfinite(tofTime[station]) ||
                !(tofEnergyDeposited[station] > 0) || !std::isfinite(tofEnergyDeposited[station]) ||
                !std::isfinite(tofLength[station]))
            {
                valid = false;
                break;
            }
            event.hitTime[station] = tofTime[station];
            event.energyDeposited[station] = tofEnergyDeposited[station] * 1e-3;
            event.pathLength[station] = tofLength[station];
            for (int coordinate = 0; coordinate < 3; ++coordinate)
                event.tofPosition[station][coordinate] = tofPosition[station][coordinate];
        }
        if (!valid)
            continue;
        ++validFourHitEntries;

        for (int layer = 1; layer <= 7; ++layer)
        {
            if (!(trackerEnergyDeposited[layer] > 0) ||
                !std::isfinite(trackerEnergyDeposited[layer]))
            {
                valid = false;
                break;
            }
            event.trackerEnergyDeposited[layer] = trackerEnergyDeposited[layer] * 1e-3;
            for (int coordinate = 0; coordinate < 3; ++coordinate)
            {
                if (!std::isfinite(trackerPosition[layer][coordinate]))
                    valid = false;
                event.trackerPosition[layer][coordinate] = trackerPosition[layer][coordinate];
            }
        }
        if (!valid)
            continue;
        ++trackerCompleteEntries;

        event.mass = mass;
        event.mcBeta = momentum / std::sqrt(momentum * momentum + mass * mass);
        if (!(event.mcBeta > 0.2 && event.mcBeta < betaMax))
            continue;
        ++betaSelectedEntries;

        if (!buildCheckpointTruthTimes(mevmom1, mevcoo1, tofPosition,
                                       event.pathLength, mass, event.checkpointTruthTime))
            continue;
        ++truthIntegrableEntries;
        if (validator.IsValidAt(event, zetaMin) && validator.IsValidAt(event, zetaMax))
            stableEvents.push_back(event);

        const int charge = static_cast<int>(mch + 0.5);
        if (fittedCharge == 0)
            fittedCharge = charge;
        else if (charge != fittedCharge)
            mixedCharge = true;
    }

    chain.Reset();

    GlobalZetaFitter measuredFitter(stableEvents, energyLossScaleMode, referencePoint,
                                    GlobalZetaTarget::MeasuredTime);
    GlobalZetaFitter truthFitter(stableEvents, energyLossScaleMode, referencePoint,
                                 GlobalZetaTarget::CheckpointTruth);
    const GlobalZetaResult measuredResult = measuredFitter.Fit(zetaMin, zetaMax);
    const GlobalZetaResult truthResult = truthFitter.Fit(zetaMin, zetaMax);
    if (!measuredResult.valid || !truthResult.valid)
    {
        std::cerr << "Error: Global energy-loss scale fit failed" << std::endl;
        return false;
    }

    TFile output(outputFile.c_str(), "RECREATE");
    if (output.IsZombie())
    {
        std::cerr << "Error: Could not create output file " << outputFile << std::endl;
        return false;
    }

    double fittedZeta = measuredResult.zeta;
    double fittedZetaMeasured = measuredResult.zeta;
    double fittedZetaTruth = truthResult.zeta;
    double chi2 = measuredResult.chi2;
    double chi2PerEvent = measuredResult.chi2PerEvent;
    double chi2Truth = truthResult.chi2;
    double chi2TruthPerEvent = truthResult.chi2PerEvent;
    double storedBetaMax = betaMax;
    double storedZetaMin = zetaMin;
    double storedZetaMax = zetaMax;
    Long64_t storedInputEntries = inputEntries;
    Long64_t storedFourHitEntries = validFourHitEntries;
    Long64_t storedTrackerCompleteEntries = trackerCompleteEntries;
    Long64_t storedBetaSelectedEntries = betaSelectedEntries;
    Long64_t storedTruthIntegrableEntries = truthIntegrableEntries;
    Long64_t storedStableEntries = measuredResult.entries;
    int storedFilesAdded = filesAdded;
    int storedCharge = mixedCharge ? -1 : fittedCharge;
    int storedScaleMode = static_cast<int>(energyLossScaleMode);
    int storedReferencePoint = static_cast<int>(referencePoint);
    const double boundaryTolerance = (zetaMax - zetaMin) / 40.0;
    int atBoundary = fittedZeta <= zetaMin + boundaryTolerance ||
                     fittedZeta >= zetaMax - boundaryTolerance;
    int truthAtBoundary = fittedZetaTruth <= zetaMin + boundaryTolerance ||
                          fittedZetaTruth >= zetaMax - boundaryTolerance;

    TTree summary("globalScaleTree", "Global Energy Loss Scale Fit");
    summary.Branch("energyLossScale", &fittedZeta, "energyLossScale/D");
    summary.Branch("energyLossScaleMeasured", &fittedZetaMeasured, "energyLossScaleMeasured/D");
    summary.Branch("energyLossScaleTruth", &fittedZetaTruth, "energyLossScaleTruth/D");
    summary.Branch("chi2", &chi2, "chi2/D");
    summary.Branch("chi2PerEvent", &chi2PerEvent, "chi2PerEvent/D");
    summary.Branch("chi2Truth", &chi2Truth, "chi2Truth/D");
    summary.Branch("chi2TruthPerEvent", &chi2TruthPerEvent, "chi2TruthPerEvent/D");
    summary.Branch("betaMax", &storedBetaMax, "betaMax/D");
    summary.Branch("zetaMin", &storedZetaMin, "zetaMin/D");
    summary.Branch("zetaMax", &storedZetaMax, "zetaMax/D");
    summary.Branch("inputEntries", &storedInputEntries, "inputEntries/L");
    summary.Branch("fourHitEntries", &storedFourHitEntries, "fourHitEntries/L");
    summary.Branch("trackerCompleteEntries", &storedTrackerCompleteEntries, "trackerCompleteEntries/L");
    summary.Branch("betaSelectedEntries", &storedBetaSelectedEntries, "betaSelectedEntries/L");
    summary.Branch("truthIntegrableEntries", &storedTruthIntegrableEntries, "truthIntegrableEntries/L");
    summary.Branch("stableEntries", &storedStableEntries, "stableEntries/L");
    summary.Branch("filesAdded", &storedFilesAdded, "filesAdded/I");
    summary.Branch("charge", &storedCharge, "charge/I");
    summary.Branch("energyLossScaleMode", &storedScaleMode, "energyLossScaleMode/I");
    summary.Branch("referencePoint", &storedReferencePoint, "referencePoint/I");
    summary.Branch("atBoundary", &atBoundary, "atBoundary/I");
    summary.Branch("truthAtBoundary", &truthAtBoundary, "truthAtBoundary/I");
    summary.Fill();
    summary.Write();

    TTree profile("globalScaleProfile", "Global Energy Loss Scale Profile");
    double profileZeta = 0;
    double profileChi2 = 0;
    double deltaChi2 = 0;
    double profileChi2Truth = 0;
    double deltaChi2Truth = 0;
    profile.Branch("energyLossScale", &profileZeta, "energyLossScale/D");
    profile.Branch("chi2", &profileChi2, "chi2/D");
    profile.Branch("deltaChi2", &deltaChi2, "deltaChi2/D");
    profile.Branch("chi2Truth", &profileChi2Truth, "chi2Truth/D");
    profile.Branch("deltaChi2Truth", &deltaChi2Truth, "deltaChi2Truth/D");
    const int profilePoints = 41;
    for (int point = 0; point < profilePoints; ++point)
    {
        profileZeta = zetaMin + (zetaMax - zetaMin) * point / (profilePoints - 1);
        profileChi2 = measuredFitter.Chi2(profileZeta);
        deltaChi2 = profileChi2 - measuredResult.chi2;
        profileChi2Truth = truthFitter.Chi2(profileZeta);
        deltaChi2Truth = profileChi2Truth - truthResult.chi2;
        profile.Fill();
    }
    profile.Write();
    output.Close();

    std::cout << "Global measured-time zeta fit: " << fittedZetaMeasured
              << ", checkpoint-truth zeta fit: " << fittedZetaTruth
              << ", measured chi2/event: " << chi2PerEvent
              << ", stable entries: " << storedStableEntries
              << ", beta-selected entries: " << storedBetaSelectedEntries
              << std::endl;
    return true;
}

bool Util::benchmarkBetaNL(const std::string &inputFile,
                           const std::string &outputFile,
                           double energyLossScale,
                           EnergyLossScaleMode energyLossScaleMode,
                           BetaReferencePoint referencePoint)
{
    // Load particle data from input file
    std::vector<ParticleData> particles = Util::loadParticleData(inputFile, referencePoint);
    if (particles.empty())
    {
        std::cerr << "Error: No particles loaded from input file" << std::endl;
        return false;
    }

    std::cout << "Starting BetaNL::Beta() benchmark with " << particles.size() << " particles..." << std::endl;
    std::cout << "Energy loss scale: " << energyLossScale << std::endl;

    double totalTimeS = 0;
    double averageTimeMs = 0;

    for (const auto &particle : particles)
    {
        auto start = std::chrono::high_resolution_clock::now();

        double betaValue =
            BetaNL(
                BetaNLPars(
                    particle.betaLinear,
                    particle.mass,
                    particle.TOF_hitEdep,
                    particle.TOF_hitTime,
                    particle.TOF_hitTimeError,
                    particle.TOF_length),
                energyLossScale,
                energyLossScaleMode,
                referencePoint)
                .Beta();

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        std::cout << "Final beta value: " << betaValue << std::endl;

        totalTimeS += elapsed.count();
    }

    averageTimeMs = totalTimeS * 1000 / particles.size();

    // Display results
    std::cout << "Benchmark completed with " << particles.size() << " particles" << std::endl;
    std::cout << "Total time: " << totalTimeS << " seconds" << std::endl;
    std::cout << "Average time per Beta() call: " << averageTimeMs << " milliseconds" << std::endl;

    return true;
}

bool Util::saveBetaDiff(const std::string &inputFile,
                        const std::string &outputFile,
                        double energyLossScale,
                        EnergyLossScaleMode energyLossScaleMode,
                        BetaReferencePoint referencePoint)
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

    float tof_betah{};
    float mevmom1[21]{};
    float mevcoo1[21][3]{};
    float tk_q[2]{};
    float tk_qin[2][3]{};
    float tof_edep[4]{};
    float tof_tl[4]{};
    float tof_tres[4]{};
    float tof_leng[4]{};

    treeIn->SetBranchAddress("tof_betah", &tof_betah);
    treeIn->SetBranchAddress("mevmom1", mevmom1);
    treeIn->SetBranchAddress("mevcoo1", mevcoo1);
    treeIn->SetBranchAddress("tk_q", tk_q);
    treeIn->SetBranchAddress("tk_qin", tk_qin);
    treeIn->SetBranchAddress("tof_edep", tof_edep);
    treeIn->SetBranchAddress("tof_tl", tof_tl);
    treeIn->SetBranchAddress("tof_tres", tof_tres);
    treeIn->SetBranchAddress("tof_leng", tof_leng);

    TFile *fileOut = new TFile(outputFile.c_str(), "RECREATE");
    if (!fileOut || fileOut->IsZombie())
    {
        std::cerr << "Error: Could not create output file " << outputFile << std::endl;
        return false;
    }

    TTree *treeOut = new TTree("betaDiff", "beta difference");

    float linearBeta{};
    float nonlinearBeta{};
    float mcBeta[21]{};

    treeOut->Branch("linearBeta", &linearBeta, "linearBeta/F");
    treeOut->Branch("nonlinearBeta", &nonlinearBeta, "nonlinearBeta/F");
    treeOut->Branch("mcBeta", mcBeta, "mcBeta[21]/F");
    treeOut->Branch("mevmom1", mevmom1, "mevmom1[21]/F");
    treeOut->Branch("mevcoo1", mevcoo1, "mevcoo1[21][3]/F");
    treeOut->Branch("tof_tres", tof_tres, "tof_tres[4]/F");
    float tof_etl[4] = {0.1544809, 0.1544809, 0.1544809, 0.1544809};

    for (Long64_t i = 0; i < treeIn->GetEntries(); ++i)
    {
        treeIn->GetEntry(i);

        float charge = (int)((tk_qin[0][2] < 2.5 ? tk_q[1] : tk_qin[0][2]) + 0.5);
        float mass = 2 * charge * 0.9314941;

        linearBeta = tof_betah;
        nonlinearBeta =
            BetaNL(
                BetaNLPars(
                    tof_betah,
                    mass,
                    tof_edep,
                    tof_tl,
                    tof_etl,
                    tof_leng),
                energyLossScale,
                energyLossScaleMode,
                referencePoint)
                .Beta();

        for (int j = 0; j < 21; ++j)
            if (mevmom1[j] != -1000)
                mcBeta[j] = mevmom1[j] / sqrt(mevmom1[j] * mevmom1[j] + mass * mass);
            else
                mcBeta[j] = -1000;

        treeOut->Fill();
    }

    treeOut->Write();
    fileOut->Close();
    fileIn->Close();

    return true;
}

void Util::test()
{
    AMSChain *ch0 = new AMSChain();
    // ch0->Add("root/be10.82318769.00000001.root");
    ch0->Add("/lustre02/data/ams/MC/AMS02/2022.v6.00/Li.B1308/li6.pl1.l1.36000.6_05/1686725052.00000001.root");
    std::cout << ch0->GetEntries() << std::endl;
    for (int ievt = 0; ievt < ch0->GetEntries(); ievt++)
    {
        AMSEventR *evt = ch0->GetEvent(ievt);

        int ntrack = evt->NTrTrack();
        if (ntrack != 2)
            continue;

        int nparticle = evt->NParticle();
        if (nparticle != 1)
            continue;

        auto particles = evt->Particle();
        if (particles[0].iCharge() != 3)
            continue;

        std::cout << evt->Event() << std::endl;
    }
}
