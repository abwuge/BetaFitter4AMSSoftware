#include "Util.hh"
#include "amschain.h"
#include "HybridEnergyLoss.h"
#include "ParticlePropagator.hh"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TPolyMarker3D.h>
#include <TPolyLine3D.h>
#include <TCanvas.h>
#include <TView.h>
#include <TText.h>
#include <chrono>

namespace
{
bool getHybridElementAbundance(
    const ParticleData &particle,
    double abundance[HybridEnergyLossCalibration::nSegments]
                    [HybridEnergyLossCalibration::nElements])
{
    if (!(std::abs(particle.innerRigidity) > 0) ||
        !std::isfinite(particle.innerRigidity))
        return false;

    double boundaryZ[5]{};
    double tofCenterZ[4]{};
    if (!HybridEnergyLossCalibration::RegionGeometry(
            particle.TOF_bar, boundaryZ, tofCenterZ))
        return false;

    const AMSPoint point(particle.initCoo[0], particle.initCoo[1],
                         particle.initCoo[2]);
    const AMSDir direction(particle.initDir[0], particle.initDir[1],
                           particle.initDir[2]);
    for (int segment = 0; segment < HybridEnergyLossCalibration::nSegments;
         ++segment)
    {
        const int status = Betalhd::betalhd_GetElementAbundance(
            point, direction, particle.innerRigidity,
            boundaryZ[segment], boundaryZ[segment + 1], abundance[segment]);
        if (status != 0)
            return false;
    }
    return true;
}
}

float Util::getAverageMass(int charge)
{
    // Fixed isotope mixtures used by the AMS analysis. Values are average
    // atomic masses in u; multiplying by u converts them to GeV/c^2.
    static const double averageAtomicMassU[] = {
        0.0,
        1.00782500,  // Z=1:  100% 1H
        3.90394781,  // Z=2:   10% 3He + 90% 4He
        6.51556315,  // Z=3:   50% 6Li + 50% 7Li
        7.91516552,  // Z=4:   60% 7Be + 30% 9Be + 10% 10Be
        10.71039495, // Z=5:   30% 10B + 70% 11B
        12.10033548, // Z=6:   90% 12C + 10% 13C
        14.50159145, // Z=7:   50% 14N + 50% 15N
        15.99491460, // Z=8:  100% 16O
        18.99840320, // Z=9:  100% 19F
        20.79201832  // Z=10: 60% 20Ne + 40% 22Ne
    };
    static const double atomicMassUnitGeV = 0.9314941;
    const int maximumCharge = sizeof(averageAtomicMassU) / sizeof(averageAtomicMassU[0]) - 1;
    if (charge < 1 || charge > maximumCharge)
        return 0;
    return averageAtomicMassU[charge] * atomicMassUnitGeV;
}

std::vector<ParticleData> Util::loadParticleData(const std::string &inputFile,
                                                  BetaReferencePoint referencePoint)
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
    float mmass{};
    float mevcoo1[21][3]{};
    float mevdir1[21][3]{};
    float mevmom1[21]{};
    float tof_betah{};
    float tof_tl[4]{};
    float tof_pos[4][3]{};
    float tof_dir[4][3]{};
    float tof_edep[4]{};
    int tof_barid[4] = {-1, -1, -1, -1};
    float tk_pos[9][3]{};
    float tk_edep[9]{};
    float tk_q[2]{};
    float tk_qin[2][3]{};
    float tk_rigidity{};
    float tof_leng[4]{};

    // Set branch addresses
    if (tree->GetBranch("mpar"))
    {
        if (!tree->GetBranch("mmass"))
        {
            std::cerr << "Error: Missing required MC truth mmass branch in "
                      << inputFile << std::endl;
            file->Close();
            return particles;
        }
        isMC = true;
        tree->SetBranchAddress("mpar", &mpar);
        tree->SetBranchAddress("mch", &mch);
        tree->SetBranchAddress("mmass", &mmass);
        tree->SetBranchAddress("mevcoo1", mevcoo1);
        tree->SetBranchAddress("mevdir1", mevdir1);
        tree->SetBranchAddress("mevmom1", mevmom1);
    }
    tree->SetBranchAddress("tof_betah", &tof_betah);
    tree->SetBranchAddress("tof_tl", tof_tl);
    tree->SetBranchAddress("tof_pos", tof_pos);
    tree->SetBranchAddress("tof_dir", tof_dir);
    tree->SetBranchAddress("tof_edep", tof_edep);
    if (tree->GetBranch("tof_barid"))
        tree->SetBranchAddress("tof_barid", tof_barid);
    tree->SetBranchAddress("tk_pos", tk_pos);
    if (!tree->GetBranch("tk_edep"))
    {
        std::cerr << "Error: Missing required tk_edep branch in " << inputFile << std::endl;
        file->Close();
        return particles;
    }
    tree->SetBranchAddress("tk_edep", tk_edep);
    tree->SetBranchAddress("tk_q", tk_q);
    tree->SetBranchAddress("tk_qin", tk_qin);
    if (!tree->GetBranch("tk_rigidity"))
    {
        std::cerr << "Error: Missing required amsd69n tk_rigidity branch in "
                  << inputFile << std::endl;
        file->Close();
        return particles;
    }
    tree->SetBranchAddress("tk_rigidity", &tk_rigidity);
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
            data.TOF_bar[j] = tof_barid[j];
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
        data.innerRigidity = tk_rigidity;
        data.mass = getAverageMass(data.charge);
        if (!(data.mass > 0))
            continue;
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
                    double trackerEnergyLossScale,
                    EnergyLossScaleMode energyLossScaleMode,
                    BetaReferencePoint referencePoint,
                    const std::string &scaleConfig,
                    const std::string &hybridConfig,
                    double materialStep)
{
    std::map<int, std::pair<double, double>> chargeScales;
    if (scaleConfig != "none")
    {
        std::ifstream config(scaleConfig.c_str());
        if (!config)
        {
            std::cerr << "Error: Could not open scale configuration " << scaleConfig << std::endl;
            return false;
        }

        std::string line;
        int lineNumber = 0;
        while (std::getline(config, line))
        {
            ++lineNumber;
            const std::string::size_type first = line.find_first_not_of(" \t\r");
            if (first == std::string::npos || line[first] == '#')
                continue;

            int charge;
            double tofScale;
            double trackerScale;
            std::string extra;
            std::istringstream values(line);
            if (!(values >> charge >> tofScale >> trackerScale) ||
                (values >> extra) || charge <= 0 ||
                !std::isfinite(tofScale) || !std::isfinite(trackerScale) ||
                tofScale <= 0 || trackerScale <= 0 || chargeScales.count(charge))
            {
                std::cerr << "Error: Invalid scale configuration at " << scaleConfig
                          << ":" << lineNumber << std::endl;
                return false;
            }
            chargeScales[charge] = std::make_pair(tofScale, trackerScale);
        }
        if (chargeScales.empty())
        {
            std::cerr << "Error: No charge scales found in " << scaleConfig << std::endl;
            return false;
        }
    }

    std::shared_ptr<HybridEnergyLossCalibration> hybridCalibration;
    if (hybridConfig != "none")
    {
        hybridCalibration = std::make_shared<HybridEnergyLossCalibration>();
        if (!hybridCalibration->Load(hybridConfig))
            return false;
    }

    // Load particle data from input file
    std::vector<ParticleData> particles =
        Util::loadParticleData(inputFile, referencePoint);
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
    double hybridBeta = std::numeric_limits<double>::quiet_NaN();
    int hybridFitStatus = -1;
    double hybridChi2 = 0;
    int hybridSteps = 0;
    int Z = 0;

    // Set up branches
    tree->Branch("mcBeta", &mcBeta, "mcBeta/D");
    tree->Branch("innerRigidity", &innerRigidity, "innerRigidity/D");
    tree->Branch("betaRigidity", &betaRigidity, "betaRigidity/D");
    tree->Branch("linearBeta", &linearBeta, "linearBeta/D");
    tree->Branch("nonlinearBetaLegacy", &nonlinearBetaLegacy, "nonlinearBetaLegacy/D");
    tree->Branch("nonlinearBeta", &nonlinearBeta, "nonlinearBeta/D");
    tree->Branch("hybridBeta", &hybridBeta, "hybridBeta/D");
    tree->Branch("hybridFitStatus", &hybridFitStatus, "hybridFitStatus/I");
    tree->Branch("hybridChi2", &hybridChi2, "hybridChi2/D");
    tree->Branch("hybridSteps", &hybridSteps, "hybridSteps/I");
    tree->Branch("Z", &Z, "Z/I");
    int mcCharge = 0;
    tree->Branch("mcCharge", &mcCharge, "mcCharge/I");

    // Process each particle
    for (const auto &particle : particles)
    {
        double particleTOFScale = energyLossScale;
        double particleTrackerScale = trackerEnergyLossScale;
        if (scaleConfig != "none")
        {
            const auto scale = chargeScales.find(particle.charge);
            if (scale == chargeScales.end())
                continue;
            particleTOFScale = scale->second.first;
            particleTrackerScale = scale->second.second;
        }

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
                particleTOFScale,
                1.0,
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
                particleTOFScale,
                particleTrackerScale,
                energyLossScaleMode,
                referencePoint)
                .Beta();
        hybridBeta = std::numeric_limits<double>::quiet_NaN();
        hybridFitStatus = -1; hybridChi2 = std::numeric_limits<double>::quiet_NaN(); hybridSteps=0;
        if (hybridCalibration)
        {
            double elementAbundance[HybridEnergyLossCalibration::nSegments]
                                   [HybridEnergyLossCalibration::nElements]{};
            HybridMaterialProfile profile;
            bool materialValid = false;
            if (materialStep > 0)
            {
                profile = HybridEnergyLossCalibration::BuildProfile(
                    particle.initCoo, particle.initDir, particle.innerRigidity,
                    particle.TOF_bar, particle.TOF_length,
                    particle.TRACKER_hitZ, materialStep);
                materialValid = profile.valid;
                for(int s=0;s<4;++s) for(int e=0;e<9;++e)
                    elementAbundance[s][e]=profile.abundance[s][e];
            }
            else materialValid = getHybridElementAbundance(particle, elementAbundance);
            if (materialValid && hybridCalibration->SupportsCharge(particle.charge))
            {
                BetaNLPars hybridPars(
                                     particle.betaLinear,
                                     particle.mass,
                                     particle.TOF_hitEdep,
                                     particle.TOF_hitTime,
                                     particle.TOF_hitTimeError,
                                     particle.TOF_length,
                                     tofPosition,
                                     particle.TRACKER_hitEdep,
                                     trackerPosition,
                                     particle.charge,
                                     elementAbundance);
                if (materialStep>0) hybridPars.SetMaterialProfile(profile);
                BetaNL hybridFit(hybridPars,
                                 particleTOFScale,
                                 particleTrackerScale,
                                 energyLossScaleMode,
                                 referencePoint,
                                 hybridCalibration);
                hybridBeta=hybridFit.Beta();
                hybridFitStatus=hybridFit.FitStatus();
                hybridChi2=hybridFit.Chi2();
                hybridSteps=profile.steps.size();
            }
        }
        Z = particle.charge;
        mcCharge = particle.mcCharge;

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
    float mmass = 0.0f;
    float mevmom1[21]{};
    float tof_edep[4]{};
    float mevcoo1[21][3]{};
    float mevdir1[21][3]{};

    treeIn->SetBranchAddress("tof_qs", &tof_qs);
    treeIn->SetBranchAddress("mpar", &mpar);
    treeIn->SetBranchAddress("mch", &mch);
    if (!treeIn->GetBranch("mmass"))
    {
        std::cerr << "Error: Missing required MC truth mmass branch in " << inputFile << std::endl;
        fileIn->Close();
        return false;
    }
    treeIn->SetBranchAddress("mmass", &mmass);
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

        double mass = mmass;

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
    TTree *tree = new TTree("scaleTree", "Per-event Common Energy Loss Scale Fit");

    float energyLossScale;
    float timeOffset;
    float chi2;
    int fitStatus;
    int fitValid;
    float mcBeta;
    float position[3];
    float direction[3];

    tree->Branch("energyLossScale", &energyLossScale, "energyLossScale/F");
    tree->Branch("timeOffset", &timeOffset, "timeOffset/F");
    tree->Branch("chi2", &chi2, "chi2/F");
    tree->Branch("fitStatus", &fitStatus, "fitStatus/I");
    tree->Branch("fitValid", &fitValid, "fitValid/I");
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

        const EnergyLossFitResult fit = BetaNL(
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
            2,
            1.0,
            energyLossScaleMode,
            referencePoint)
            .EnergyLossScales(particle.mcBeta);
        energyLossScale = fit.energyLossScale;
        timeOffset = fit.timeOffset;
        chi2 = fit.chi2;
        fitStatus = fit.status;
        fitValid = fit.valid;
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

bool Util::saveHybridEnergyLossTraining(const std::string &inputFile,
                                        const std::string &outputFile,
                                        double materialStep)
{
    TFile *fileIn = TFile::Open(inputFile.c_str(), "READ");
    if (!fileIn || fileIn->IsZombie())
    {
        std::cerr << "Error: Could not open file " << inputFile << std::endl;
        return false;
    }

    TTree *treeIn = static_cast<TTree *>(fileIn->Get("amstreea"));
    if (!treeIn)
    {
        std::cerr << "Error: Could not find amstreea tree in " << inputFile << std::endl;
        fileIn->Close();
        return false;
    }

    const char *requiredBranches[] = {
        "mpar", "mch", "mmass", "mevmom1", "mevcoo1", "mevdir1",
        "tof_pos", "tof_dir", "tof_edep", "tof_barid", "tk_edep",
        "tk_q", "tk_qin", "tk_rigidity", "tof_leng", "tk_pos",
        "mtofdep", "mtofpri", "mtofbar", "mtrmom", "mtrpri", "mtrcoo1"};
    for (const char *branch : requiredBranches)
    {
        if (!treeIn->GetBranch(branch))
        {
            std::cerr << "Error: Missing required hybrid-training branch "
                      << branch << " in " << inputFile << std::endl;
            fileIn->Close();
            return false;
        }
    }

    int mpar = 0;
    float mch = 0;
    float mmass = 0;
    float mevmom1[21]{};
    float mevcoo1[21][3]{};
    float mevdir1[21][3]{};
    float tofPosition[ParticleData::TOF_MAX_HITS][3]{};
    float tofDirection[ParticleData::TOF_MAX_HITS][3]{};
    float tofEdepInput[ParticleData::TOF_MAX_HITS]{};
    int tofBar[ParticleData::TOF_MAX_HITS]{};
    float mcTofDepositInput[ParticleData::TOF_MAX_HITS]{};
    int mcTofPrimary[ParticleData::TOF_MAX_HITS]{};
    int mcTofBar[ParticleData::TOF_MAX_HITS]{};
    float trackerEdepInput[ParticleData::TRACKER_MAX_HITS]{};
    float trackerCharge[2]{};
    float innerTrackerCharge[2][3]{};
    float trackerRigidity = 0;
    float tofLength[4]{};
    float trackerPosition[9][3]{};
    float mcTrackerMomentum[ParticleData::TRACKER_MAX_HITS]{};
    int mcTrackerPrimary[ParticleData::TRACKER_MAX_HITS]{};
    float mcTrackerCenter[ParticleData::TRACKER_MAX_HITS][3]{};

    treeIn->SetBranchAddress("mpar", &mpar);
    treeIn->SetBranchAddress("mch", &mch);
    treeIn->SetBranchAddress("mmass", &mmass);
    treeIn->SetBranchAddress("mevmom1", mevmom1);
    treeIn->SetBranchAddress("mevcoo1", mevcoo1);
    treeIn->SetBranchAddress("mevdir1", mevdir1);
    treeIn->SetBranchAddress("tof_pos", tofPosition);
    treeIn->SetBranchAddress("tof_dir", tofDirection);
    treeIn->SetBranchAddress("tof_edep", tofEdepInput);
    treeIn->SetBranchAddress("tof_barid", tofBar);
    treeIn->SetBranchAddress("mtofdep", mcTofDepositInput);
    treeIn->SetBranchAddress("mtofpri", mcTofPrimary);
    treeIn->SetBranchAddress("mtofbar", mcTofBar);
    treeIn->SetBranchAddress("tk_edep", trackerEdepInput);
    treeIn->SetBranchAddress("tk_q", trackerCharge);
    treeIn->SetBranchAddress("tk_qin", innerTrackerCharge);
    treeIn->SetBranchAddress("tk_rigidity", &trackerRigidity);
    treeIn->SetBranchAddress("tof_leng", tofLength);
    treeIn->SetBranchAddress("tk_pos", trackerPosition);
    treeIn->SetBranchAddress("mtrmom", mcTrackerMomentum);
    treeIn->SetBranchAddress("mtrpri", mcTrackerPrimary);
    treeIn->SetBranchAddress("mtrcoo1", mcTrackerCenter);

    TFile *fileOut = TFile::Open(outputFile.c_str(), "RECREATE");
    if (!fileOut || fileOut->IsZombie())
    {
        std::cerr << "Error: Could not create output file " << outputFile << std::endl;
        fileIn->Close();
        return false;
    }

    TTree *treeOut = new TTree(
        "hybridTraining",
        "Material-mean plus measured-fluctuation training samples");

    Long64_t sourceEntry = -1;
    int geantId = 0;
    int charge = 0;
    float mass = 0;
    float mcBeta = 0;
    float position[3]{};
    float direction[3]{};
    float segmentZ[5]{};
    float segmentBeta[4]{};
    float materialMean[4]{};
    float truthLoss[4]{};
    float measuredDeposit[4]{};
    float tofDeposit[ParticleData::TOF_MAX_HITS]{};
    float mcTofDeposit[ParticleData::TOF_MAX_HITS]{};
    float truthBoundaryKinetic[5]{};
    float truthCheckpointZ[5]{};
    int truthCheckpoint[5]{};
    float trackerDeposit[ParticleData::TRACKER_MAX_HITS]{};
    double elementAbundance[HybridEnergyLossCalibration::nSegments]
                           [HybridEnergyLossCalibration::nElements]{};

    treeOut->Branch("sourceEntry", &sourceEntry, "sourceEntry/L");
    treeOut->Branch("geantId", &geantId, "geantId/I");
    treeOut->Branch("Z", &charge, "Z/I");
    treeOut->Branch("mass", &mass, "mass/F");
    treeOut->Branch("mcBeta", &mcBeta, "mcBeta/F");
    treeOut->Branch("position", position, "position[3]/F");
    treeOut->Branch("direction", direction, "direction[3]/F");
    treeOut->Branch("segmentZ", segmentZ, "segmentZ[5]/F");
    treeOut->Branch("segmentBeta", segmentBeta, "segmentBeta[4]/F");
    treeOut->Branch("materialMean", materialMean, "materialMean[4]/F");
    treeOut->Branch("truthLoss", truthLoss, "truthLoss[4]/F");
    treeOut->Branch("measuredDeposit", measuredDeposit, "measuredDeposit[4]/F");
    treeOut->Branch("tofDeposit", tofDeposit, "tofDeposit[4]/F");
    treeOut->Branch("tofBar", tofBar, "tofBar[4]/I");
    treeOut->Branch("mcTofDeposit", mcTofDeposit, "mcTofDeposit[4]/F");
    treeOut->Branch("mcTofBar", mcTofBar, "mcTofBar[4]/I");
    treeOut->Branch("mcTofPrimary", mcTofPrimary, "mcTofPrimary[4]/I");
    treeOut->Branch("truthBoundaryKinetic", truthBoundaryKinetic,
                    "truthBoundaryKinetic[5]/F");
    treeOut->Branch("truthCheckpointZ", truthCheckpointZ,
                    "truthCheckpointZ[5]/F");
    treeOut->Branch("truthCheckpoint", truthCheckpoint,
                    "truthCheckpoint[5]/I");
    treeOut->Branch("trackerDeposit", trackerDeposit, "trackerDeposit[9]/F");
    treeOut->Branch("elementAbundance", elementAbundance,
                    "elementAbundance[4][9]/D");


    for (Long64_t entry = 0; entry < treeIn->GetEntries(); ++entry)
    {
        treeIn->GetEntry(entry);

        if (!(std::abs(mch) > 0) || !(mmass > 0) ||
            !std::isfinite(mch) || !std::isfinite(mmass))
            continue;
        const int truthCharge = static_cast<int>(std::lround(std::abs(mch)));
        const int reconstructedCharge = static_cast<int>(
            (innerTrackerCharge[0][2] < 2.5 ? trackerCharge[1]
                                            : innerTrackerCharge[0][2]) +
            0.5);
        if (reconstructedCharge != truthCharge)
            continue;

        bool valid = true;
        double boundaryZ[5]{};
        double tofCenterZ[4]{};
        if (!HybridEnergyLossCalibration::RegionGeometry(
                tofBar, boundaryZ, tofCenterZ))
            continue;
        for (int boundary = 0; boundary < 5; ++boundary)
            segmentZ[boundary] = boundaryZ[boundary];

        // Training uses a single, auditable bar identity for geometry, MC
        // surface labels, and active-detector deposits.  Mismatched reco/MC
        // bars are excluded rather than mixing two different material paths.
        for (int station = 0; station < ParticleData::TOF_MAX_HITS; ++station)
        {
            valid = valid && tofBar[station] == mcTofBar[station] &&
                    mcTofPrimary[station] == 1 &&
                    mcTofDepositInput[station] >= 0 &&
                    std::isfinite(mcTofDepositInput[station]) &&
                    tofEdepInput[station] > 0 &&
                    std::isfinite(tofEdepInput[station]);
            tofDeposit[station] = tofEdepInput[station] * 1e-3;
            mcTofDeposit[station] = mcTofDepositInput[station] * 1e-3;
        }
        for (int coordinate = 0; coordinate < 3; ++coordinate)
            valid = valid && std::isfinite(tofPosition[0][coordinate]) &&
                    std::isfinite(tofDirection[0][coordinate]);
        valid = valid && std::isfinite(trackerRigidity) &&
                std::abs(trackerRigidity) > 0;
        if (!valid)
            continue;

        // AMS G4 acceptance planes at the actual staggered TOF bar centers.
        // The recorded checkpoint is the post-step state after crossing that
        // plane, so truthCheckpointZ is retained for explicit audit.
        truthCheckpoint[0] = (tofBar[0] % 2 == 0) ? 3 : 4;   // S1 center
        truthCheckpoint[1] = (tofBar[1] % 2 == 0) ? 5 : 6;   // S2 center
        truthCheckpoint[2] = -1;                              // primary L5 TrMCCluster
        truthCheckpoint[3] = (tofBar[2] % 2 == 0) ? 14 : 15; // S3 center
        truthCheckpoint[4] = 17;                              // post-S4 G4 acceptance plane

        double checkpointKinetic[5]{};
        for (int boundary = 0; boundary < 5; ++boundary)
        {
            if (boundary == 2)
            {
                valid = valid && mcTrackerPrimary[4] == 1 &&
                        mcTrackerMomentum[4] > 0 &&
                        std::isfinite(mcTrackerMomentum[4]) &&
                        std::isfinite(mcTrackerCenter[4][2]) &&
                        std::abs(mcTrackerCenter[4][2] - boundaryZ[2]) < 0.03;
                if (!valid)
                    break;
                checkpointKinetic[boundary] =
                    std::sqrt(mcTrackerMomentum[4] * mcTrackerMomentum[4] +
                              mmass * mmass) - mmass;
                truthCheckpointZ[boundary] = mcTrackerCenter[4][2];
                continue;
            }

            const int index = truthCheckpoint[boundary];
            valid = valid && mevmom1[index] > 0 &&
                    std::isfinite(mevmom1[index]);
            for (int coordinate = 0; coordinate < 3; ++coordinate)
                valid = valid && std::isfinite(mevcoo1[index][coordinate]) &&
                        std::isfinite(mevdir1[index][coordinate]);
            if (!valid)
                break;
            checkpointKinetic[boundary] =
                std::sqrt(mevmom1[index] * mevmom1[index] + mmass * mmass) -
                mmass;
            truthCheckpointZ[boundary] = mevcoo1[index][2];
        }
        if (!valid)
            continue;

        // Move the S1/S2/S3 center-checkpoint energies to the actual bar
        // surfaces using the primary MC deposit and the checkpoint's true z.
        // This is a local within-scintillator interpolation; unlike the old
        // fixed checkpoints, it respects each event's staggered hit bar.
        const auto traversedFraction = [](double zHigh, double zLow,
                                          double zPoint, double &fraction)
        {
            constexpr double tolerance = 0.02; // cm, only guards float/step edge noise
            if (!(zHigh > zLow) || zPoint > zHigh + tolerance ||
                zPoint < zLow - tolerance)
                return false;
            fraction = (zHigh - zPoint) / (zHigh - zLow);
            fraction = std::max(0.0, std::min(1.0, fraction));
            return std::isfinite(fraction);
        };

        const double s1High = boundaryZ[0];
        const double s1Low = 2.0 * tofCenterZ[0] - s1High;
        const double s2Low = boundaryZ[1];
        const double s2High = 2.0 * tofCenterZ[1] - s2Low;
        const double s3High = boundaryZ[3];
        const double s3Low = 2.0 * tofCenterZ[2] - s3High;
        double s1Traversed = 0;
        double s2Traversed = 0;
        double s3Traversed = 0;
        valid = traversedFraction(s1High, s1Low, truthCheckpointZ[0],
                                  s1Traversed) &&
                traversedFraction(s2High, s2Low, truthCheckpointZ[1],
                                  s2Traversed) &&
                traversedFraction(s3High, s3Low, truthCheckpointZ[3],
                                  s3Traversed);
        if (!valid)
            continue;

        truthBoundaryKinetic[0] = checkpointKinetic[0] +
                                  s1Traversed * mcTofDeposit[0];
        truthBoundaryKinetic[1] = checkpointKinetic[1] -
                                  (1.0 - s2Traversed) * mcTofDeposit[1];
        truthBoundaryKinetic[2] = checkpointKinetic[2];
        truthBoundaryKinetic[3] = checkpointKinetic[3] +
                                  s3Traversed * mcTofDeposit[2];
        truthBoundaryKinetic[4] = checkpointKinetic[4];
        for (int boundary = 0; boundary < 5; ++boundary)
            valid = valid && truthBoundaryKinetic[boundary] >= 0 &&
                    std::isfinite(truthBoundaryKinetic[boundary]);
        if (!valid)
            continue;

        for (int layer = 0; layer < ParticleData::TRACKER_MAX_HITS; ++layer)
        {
            valid = valid && std::isfinite(trackerEdepInput[layer]);
            trackerDeposit[layer] = trackerEdepInput[layer] * 1e-3;
        }
        for (int layer = 1; layer <= 7; ++layer)
            valid = valid && trackerEdepInput[layer] > 0;
        if (!valid)
            continue;

        measuredDeposit[0] = tofDeposit[0] + tofDeposit[1];
        // Region 1/2 meet at the L5 sensor reference plane, so split the L5
        // active deposit symmetrically instead of assigning the full layer to
        // either side of a boundary that runs through its center.
        measuredDeposit[1] = 0.5 * trackerDeposit[4];
        for (int layer = 1; layer <= 3; ++layer)
            measuredDeposit[1] += trackerDeposit[layer];
        measuredDeposit[2] = 0.5 * trackerDeposit[4];
        for (int layer = 5; layer <= 7; ++layer)
            measuredDeposit[2] += trackerDeposit[layer];
        measuredDeposit[3] = tofDeposit[2] + tofDeposit[3];

        HybridMaterialProfile profile;
        if (materialStep > 0)
        {
            float trackerZ[9];
            for (int i = 0; i < 9; ++i)
                trackerZ[i] = trackerPosition[i][2];
            profile = HybridEnergyLossCalibration::BuildProfile(
                tofPosition[0], tofDirection[0], trackerRigidity, tofBar,
                tofLength, trackerZ, materialStep);
            if (!profile.valid)
                continue;
        }

        const AMSPoint point(tofPosition[0][0], tofPosition[0][1],
                             tofPosition[0][2]);
        const AMSDir dir(tofDirection[0][0], tofDirection[0][1],
                         tofDirection[0][2]);
        for (int segment = 0; segment < 4; ++segment)
        {
            truthLoss[segment] = truthBoundaryKinetic[segment] -
                                 truthBoundaryKinetic[segment + 1];
            if (!(truthLoss[segment] >= 0) || !std::isfinite(truthLoss[segment]))
            {
                valid = false;
                break;
            }

            int abundanceStatus = 0;
            if (materialStep > 0)
                for (int element = 0; element < 9; ++element)
                    elementAbundance[segment][element] =
                        profile.abundance[segment][element];
            else
                abundanceStatus = Betalhd::betalhd_GetElementAbundance(
                    point, dir, trackerRigidity, segmentZ[segment],
                    segmentZ[segment + 1], elementAbundance[segment]);

            const double totalEnergy = truthBoundaryKinetic[segment] + mmass;
            const double betaSquared =
                1.0 - mmass * mmass / (totalEnergy * totalEnergy);
            if (!(betaSquared > 0) || !(betaSquared < 1) ||
                !std::isfinite(betaSquared))
            {
                valid = false;
                break;
            }
            segmentBeta[segment] = std::sqrt(betaSquared);
            materialMean[segment] = HybridEnergyLossCalibration::MaterialMean(
                segmentBeta[segment], mmass, reconstructedCharge,
                elementAbundance[segment]);
            valid = valid && abundanceStatus == 0 && materialMean[segment] >= 0 &&
                    std::isfinite(materialMean[segment]);
        }
        if (!valid)
            continue;

        sourceEntry = entry;
        geantId = mpar;
        charge = reconstructedCharge;
        mass = mmass;
        mcBeta = segmentBeta[0];
        for (int coordinate = 0; coordinate < 3; ++coordinate)
        {
            position[coordinate] = tofPosition[0][coordinate];
            direction[coordinate] = tofDirection[0][coordinate];
        }
        treeOut->Fill();
    }

    if (treeOut->GetEntries() == 0)
    {
        std::cerr << "Error: No valid hybrid training samples saved to "
                  << outputFile << std::endl;
        fileOut->Close();
        fileIn->Close();
        delete fileOut;
        return false;
    }

    treeOut->Write();
    std::cout << "Hybrid training samples: " << treeOut->GetEntries() << std::endl;
    std::cout << "Hybrid training data saved to " << outputFile << std::endl;
    fileOut->Close();
    fileIn->Close();
    delete fileOut;
    return true;
}

bool Util::trainHybridEnergyLoss(const std::string &inputFile,
                                 const std::string &outputFile)
{
    std::unique_ptr<TFile> file;
    std::unique_ptr<TChain> chain;
    TTree *tree = nullptr;
    if (inputFile.size()>5 && inputFile.substr(inputFile.size()-5)==".list")
    {
        chain.reset(new TChain("hybridTraining"));
        std::ifstream list(inputFile);
        std::string path;
        while(std::getline(list,path))
            if(!path.empty() && path[0]!='#' && !chain->Add(path.c_str(),0))
            { std::cerr<<"Invalid training part: "<<path<<std::endl;return false; }
        chain->LoadTree(0);
        tree=chain.get();
    }
    else
    {
        file.reset(TFile::Open(inputFile.c_str(), "READ"));
        if(file && !file->IsZombie()) tree=static_cast<TTree *>(file->Get("hybridTraining"));
    }
    if (!tree)
    {
        std::cerr << "Error: Could not open hybrid training tree "
                  << inputFile << std::endl;
        return false;
    }

    int charge = 0;
    float beta = 0;
    float mass = 0;
    float segmentBeta[HybridEnergyLossCalibration::nSegments]{};
    float materialMean[HybridEnergyLossCalibration::nSegments]{};
    float truthLoss[HybridEnergyLossCalibration::nSegments]{};
    float measuredDeposit[HybridEnergyLossCalibration::nSegments]{};
    tree->SetBranchAddress("Z", &charge);
    tree->SetBranchAddress("mcBeta", &beta);
    tree->SetBranchAddress("mass", &mass);
    const bool hasSegmentBeta = tree->GetBranch("segmentBeta");
    if (hasSegmentBeta)
        tree->SetBranchAddress("segmentBeta", segmentBeta);
    tree->SetBranchAddress("materialMean", materialMean);
    tree->SetBranchAddress("truthLoss", truthLoss);
    tree->SetBranchAddress("measuredDeposit", measuredDeposit);

    struct TrainingRecord
    {
        int charge;
        double beta[HybridEnergyLossCalibration::nSegments];
        double material[HybridEnergyLossCalibration::nSegments];
        double truth[HybridEnergyLossCalibration::nSegments];
        double deposit[HybridEnergyLossCalibration::nSegments];
    };
    std::vector<TrainingRecord> records;
    typedef std::pair<int, int> ChargeSegment;
    std::map<ChargeSegment, std::vector<double>> betaValues;
    records.reserve(tree->GetEntries());

    for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry)
    {
        tree->GetEntry(entry);
        TrainingRecord record{};
        record.charge = charge;
        bool valid = charge > 0 && beta > 0 && beta < 1.0 &&
                     mass > 0 && std::isfinite(beta) && std::isfinite(mass);
        double boundaryEnergy = valid
                                    ? mass / std::sqrt(1.0 - beta * beta)
                                    : 0;
        for (int segment = 0;
             segment < HybridEnergyLossCalibration::nSegments; ++segment)
        {
            record.beta[segment] = hasSegmentBeta
                                       ? segmentBeta[segment]
                                       : (boundaryEnergy > mass
                                              ? std::sqrt(1.0 - mass * mass /
                                                                    (boundaryEnergy *
                                                                     boundaryEnergy))
                                              : 0);
            record.material[segment] = materialMean[segment];
            record.truth[segment] = truthLoss[segment];
            record.deposit[segment] = measuredDeposit[segment];
            valid = valid && record.beta[segment] > 0 &&
                    record.beta[segment] < 1.000001 &&
                    std::isfinite(record.beta[segment]) &&
                    materialMean[segment] > 0 &&
                    truthLoss[segment] >= 0 && measuredDeposit[segment] > 0 &&
                    std::isfinite(materialMean[segment]) &&
                    std::isfinite(truthLoss[segment]) &&
                    std::isfinite(measuredDeposit[segment]);
            boundaryEnergy -= truthLoss[segment];
        }
        if (!valid)
            continue;
        records.push_back(record);
        for (int segment = 0;
             segment < HybridEnergyLossCalibration::nSegments; ++segment)
            betaValues[ChargeSegment(charge, segment)].push_back(
                record.beta[segment]);
    }
    if (file) file->Close();

    if (records.empty())
    {
        std::cerr << "Error: No valid hybrid training records in "
                  << inputFile << std::endl;
        return false;
    }

    std::map<ChargeSegment, std::vector<double>> betaEdges;
    for (auto &segmentValues : betaValues)
    {
        std::vector<double> &values = segmentValues.second;
        std::sort(values.begin(), values.end());
        const int requestedBins = std::min(8, std::max(1,
            static_cast<int>(values.size() / 500)));
        std::vector<double> edges;
        edges.push_back(0.0);
        for (int bin = 1; bin < requestedBins; ++bin)
        {
            const size_t index = static_cast<size_t>(
                static_cast<double>(bin) * values.size() / requestedBins);
            const double edge = values[std::min(index, values.size() - 1)];
            if (edge > edges.back() + 1e-6)
                edges.push_back(edge);
        }
        edges.push_back(1.000001);
        betaEdges[segmentValues.first] = edges;
    }

    struct FitAccumulator
    {
        long long entries = 0;
        double sumMaterial = 0;
        double sumTruth = 0;
        double sumDeposit = 0;
        double materialScale = 0;
        double meanDeposit = 0;
        double sumX2 = 0;
        double sumResidual2 = 0;
        double sumXResidual = 0;
    };
    std::map<ChargeSegment, std::vector<FitAccumulator>> fit;
    for (const auto &edges : betaEdges)
        fit[edges.first].resize(edges.second.size() - 1);

    const auto findBin = [&](int recordCharge, int recordSegment,
                             double recordBeta)
    {
        const ChargeSegment key(recordCharge, recordSegment);
        const std::vector<double> &edges = betaEdges[key];
        const auto upper = std::upper_bound(edges.begin(), edges.end(),
                                            recordBeta);
        int bin = static_cast<int>(upper - edges.begin()) - 1;
        return std::max(0, std::min(bin,
                                   static_cast<int>(edges.size()) - 2));
    };

    std::map<ChargeSegment, std::vector<std::vector<double>>> ratioValues;
    for (const auto &edges : betaEdges)
        ratioValues[edges.first].resize(edges.second.size() - 1);
    for (const auto &record : records)
    {
        for (int segment = 0;
             segment < HybridEnergyLossCalibration::nSegments; ++segment)
        {
            const int bin = findBin(record.charge, segment,
                                    record.beta[segment]);
            ratioValues[ChargeSegment(record.charge, segment)][bin].push_back(
                record.truth[segment] / record.material[segment]);
        }
    }

    struct RobustRange
    {
        double low = 0;
        double high = std::numeric_limits<double>::max();
    };
    std::map<ChargeSegment, std::vector<RobustRange>> robustRanges;
    const auto median = [](std::vector<double> values)
    {
        std::sort(values.begin(), values.end());
        const size_t middle = values.size() / 2;
        return values.size() % 2
                   ? values[middle]
                   : 0.5 * (values[middle - 1] + values[middle]);
    };
    for (const auto &segmentRatios : ratioValues)
    {
        robustRanges[segmentRatios.first].resize(segmentRatios.second.size());
        for (size_t bin = 0; bin < segmentRatios.second.size(); ++bin)
        {
            const std::vector<double> &ratios = segmentRatios.second[bin];
            if (ratios.empty())
                continue;
            const double center = median(ratios);
            std::vector<double> deviations;
            deviations.reserve(ratios.size());
            for (double ratio : ratios)
                deviations.push_back(std::abs(ratio - center));
            const double sigma = 1.4826 * median(deviations);
            const double halfWidth = sigma > 0 ? 10 * sigma
                                               : std::max(1.0, center);
            robustRanges[segmentRatios.first][bin].low =
                std::max(0.0, center - halfWidth);
            robustRanges[segmentRatios.first][bin].high =
                center + halfWidth;
        }
    }

    const auto accepted = [&](const TrainingRecord &record, int bin,
                              int segment)
    {
        const double ratio =
            record.truth[segment] / record.material[segment];
        const RobustRange &range =
            robustRanges[ChargeSegment(record.charge, segment)][bin];
        return ratio >= range.low && ratio <= range.high;
    };

    for (const auto &record : records)
    {
        for (int segment = 0;
             segment < HybridEnergyLossCalibration::nSegments; ++segment)
        {
            const ChargeSegment key(record.charge, segment);
            const int bin = findBin(record.charge, segment,
                                    record.beta[segment]);
            if (!accepted(record, bin, segment))
                continue;
            FitAccumulator &accumulator = fit[key][bin];
            ++accumulator.entries;
            accumulator.sumMaterial += record.material[segment];
            accumulator.sumTruth += record.truth[segment];
            accumulator.sumDeposit += record.deposit[segment];
        }
    }

    for (auto &segmentFit : fit)
        for (auto &accumulator : segmentFit.second)
        {
            if (accumulator.entries == 0 ||
                !(accumulator.sumMaterial > 0))
                continue;
            accumulator.materialScale =
                accumulator.sumTruth / accumulator.sumMaterial;
            accumulator.meanDeposit =
                accumulator.sumDeposit / accumulator.entries;
        }

    for (const auto &record : records)
    {
        for (int segment = 0;
             segment < HybridEnergyLossCalibration::nSegments; ++segment)
        {
            const ChargeSegment key(record.charge, segment);
            const int bin = findBin(record.charge, segment,
                                    record.beta[segment]);
            if (!accepted(record, bin, segment))
                continue;
            FitAccumulator &accumulator = fit[key][bin];
            const double centeredDeposit =
                record.deposit[segment] - accumulator.meanDeposit;
            const double residual = record.truth[segment] -
                                    accumulator.materialScale *
                                        record.material[segment];
            accumulator.sumX2 += centeredDeposit * centeredDeposit;
            accumulator.sumResidual2 += residual * residual;
            accumulator.sumXResidual += centeredDeposit * residual;
        }
    }

    std::ofstream output(outputFile.c_str());
    if (!output)
    {
        std::cerr << "Error: Could not create hybrid calibration "
                  << outputFile << std::endl;
        return false;
    }
    output << "# betaFitter hybrid energy-loss calibration\n"
           << "# Z segment beta_low beta_high entries material_scale "
              "mean_deposit fluctuation_scale rms_before rms_after correlation\n"
           << std::setprecision(12);

    int writtenBins = 0;
    for (const auto &segmentFit : fit)
    {
        const int recordCharge = segmentFit.first.first;
        const int segment = segmentFit.first.second;
        const std::vector<double> &edges = betaEdges[segmentFit.first];
        for (size_t bin = 0; bin < segmentFit.second.size(); ++bin)
        {
                const FitAccumulator &accumulator = segmentFit.second[bin];
                if (accumulator.entries < 20 || !(accumulator.sumX2 > 0))
                    continue;
                const double fittedFluctuationScale =
                    accumulator.sumXResidual / accumulator.sumX2;
                const double fittedCorrelation =
                    accumulator.sumResidual2 > 0
                        ? accumulator.sumXResidual /
                              std::sqrt(accumulator.sumX2 *
                                        accumulator.sumResidual2)
                        : 0;
                // Deposited energy is only useful as a fluctuation estimator
                // where it has a reproducible positive relation to truth.
                // Suppress statistically weak and sign-flipped regressions;
                // otherwise they merely transfer detector noise into beta.
                const double fluctuationScale = fittedCorrelation > 0.1
                                                    ? fittedFluctuationScale
                                                    : 0.0;
                const double residualAfter2 = std::max(
                    0.0,
                    accumulator.sumResidual2 -
                        2 * fluctuationScale * accumulator.sumXResidual +
                        fluctuationScale * fluctuationScale *
                            accumulator.sumX2);
                const double rmsBefore = std::sqrt(
                    accumulator.sumResidual2 / accumulator.entries);
                const double rmsAfter = std::sqrt(
                    residualAfter2 / accumulator.entries);
                output << recordCharge << ' ' << segment << ' '
                       << edges[bin] << ' ' << edges[bin + 1] << ' '
                       << accumulator.entries << ' '
                       << accumulator.materialScale << ' '
                       << accumulator.meanDeposit << ' '
                       << fluctuationScale << ' ' << rmsBefore << ' '
                       << rmsAfter << ' ' << fittedCorrelation << '\n';
                ++writtenBins;
        }
    }

    if (writtenBins == 0)
    {
        std::cerr << "Error: Hybrid training sample is too small" << std::endl;
        return false;
    }
    std::cout << "Hybrid calibration records: " << records.size() << std::endl;
    std::cout << "Hybrid calibration bins: " << writtenBins << std::endl;
    std::cout << "Hybrid calibration saved to " << outputFile << std::endl;
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
                1.0,
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
    float mmass{};
    float mevmom1[21]{};
    float mevcoo1[21][3]{};
    float tk_q[2]{};
    float tk_qin[2][3]{};
    float tof_edep[4]{};
    float tof_tl[4]{};
    float tof_tres[4]{};
    float tof_leng[4]{};

    treeIn->SetBranchAddress("tof_betah", &tof_betah);
    if (!treeIn->GetBranch("mmass"))
    {
        std::cerr << "Error: Missing required MC truth mmass branch in " << inputFile << std::endl;
        fileIn->Close();
        return false;
    }
    treeIn->SetBranchAddress("mmass", &mmass);
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
        float mass = getAverageMass(charge);
        if (!(mass > 0))
            continue;

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
                1.0,
                energyLossScaleMode,
                referencePoint)
                .Beta();

        for (int j = 0; j < 21; ++j)
            if (mevmom1[j] != -1000)
                mcBeta[j] = mevmom1[j] / sqrt(mevmom1[j] * mevmom1[j] + mmass * mmass);
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
