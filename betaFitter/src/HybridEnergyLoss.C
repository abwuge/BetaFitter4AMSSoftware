#include "HybridEnergyLoss.h"

#include "Betalhd.h"
#include "Tofdbc.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <TFile.h>
#include <TH3.h>
#include <memory>

bool HybridEnergyLossCalibration::RegionGeometry(
    const int tofBar[4], double boundaryZ[5], double tofCenterZ[4])
{
    if (!tofBar || !boundaryZ || !tofCenterZ)
        return false;

    double surfaceLow[4]{};
    double surfaceHigh[4]{};
    for (int station = 0; station < 4; ++station)
    {
        const int bar = tofBar[station];
        if (bar < 0 || bar >= TOFGeom::Nbar[station])
            return false;
        float edge[3][2]{};
        TOFGeom::GetBarEdge(station, bar, edge);
        surfaceLow[station] = edge[2][0];
        surfaceHigh[station] = edge[2][1];
        tofCenterZ[station] = TOFGeom::Sci_pz[station][bar];
    }

    // Down-going definition.  Regions 1 and 2 are deliberately separated
    // from active TOF scintillator by the actual hit-bar surfaces.
    //   0: S1 scintillator entrance -> S2 scintillator exit
    //   1: S2 scintillator exit     -> Tracker L5 reference plane
    //   2: Tracker L5 reference     -> S3 scintillator entrance
    //   3: S3 scintillator entrance -> post-S4 G4 acceptance plane
    boundaryZ[0] = surfaceHigh[0];
    boundaryZ[1] = surfaceLow[1];
    boundaryZ[2] = 1.698;   // AMS G4 acceptance plane 10, Tracker L5
    boundaryZ[3] = surfaceHigh[2];
    boundaryZ[4] = -69.975; // AMS G4 acceptance plane 17, downstream of S4

    // The lower package boundary must remain below the full S4 scintillator.
    if (!(surfaceLow[3] > boundaryZ[4]))
        return false;
    for (int i = 0; i < 4; ++i)
        if (!(boundaryZ[i] > boundaryZ[i + 1]))
            return false;
    return true;
}

HybridMaterialProfile HybridEnergyLossCalibration::BuildProfile(
    const float position[3], const float direction[3], double rigidity,
    const int tofBar[4], const float tofLength[4],
    const float trackerZ[9], double maxStep)
{
    HybridMaterialProfile result;
    if (!std::isfinite(rigidity) || rigidity == 0 || !(maxStep > 0))
        return result;
    for (int i = 0; i < 3; ++i)
        if (!std::isfinite(position[i]) || !std::isfinite(direction[i]))
            return result;
    for (int i = 0; i < 4; ++i)
        if (!std::isfinite(tofLength[i]))
            return result;

    double boundary[5]{};
    double tofZ[4]{};
    if (!RegionGeometry(tofBar, boundary, tofZ))
        return result;
    for (int i = 0; i < 5; ++i)
        result.boundaryZ[i] = boundary[i];
    for (int i = 0; i < 4; ++i)
        result.tofCenterZ[i] = tofZ[i];
    result.centerZ = boundary[2];

    static std::unique_ptr<TFile> map;
    static TH3 *hist[9] = {};
    if (!map)
    {
        TDirectory *saved = gDirectory;
        std::string path = AMSDATADIR.amsdatadir;
        path += "LAPP/dEdxPDF/g4mscan.root";
        map.reset(TFile::Open(path.c_str(), "READ"));
        if (map && !map->IsZombie())
            for (int j = 0; j < 9; ++j)
                hist[j] = dynamic_cast<TH3 *>(map->Get(Form("hist%d", j+1)));
        if (saved) saved->cd();
    }
    for (int j = 0; j < 9; ++j)
        if (!hist[j]) return result;

    std::vector<double> nodes(boundary, boundary+5);
    for (int i = 0; i < 4; ++i) nodes.push_back(tofZ[i]);
    for (int i = 1; i <= 7; ++i)
        if (trackerZ[i] < boundary[0] && trackerZ[i] > boundary[4])
            nodes.push_back(trackerZ[i]);
    // Include every native z voxel face: neighboring pieces partition exactly.
    TAxis *axis = hist[0]->GetZaxis();
    for (int bin = 1; bin <= axis->GetNbins()+1; ++bin)
    {
        const double z = axis->GetBinLowEdge(bin);
        if (z < boundary[0] && z > boundary[4]) nodes.push_back(z);
    }
    std::sort(nodes.begin(), nodes.end(), std::greater<double>());
    nodes.erase(std::unique(nodes.begin(), nodes.end(),
                           [](double a, double b){return std::abs(a-b)<1e-7;}),
                nodes.end());
    std::vector<double> fine;
    for (size_t i = 0; i+1 < nodes.size(); ++i)
    {
        const int count = std::max(1, static_cast<int>(std::ceil((nodes[i]-nodes[i+1])/maxStep)));
        for (int k = 0; k < count; ++k)
            fine.push_back(nodes[i] + (nodes[i+1]-nodes[i])*k/count);
    }
    fine.push_back(nodes.back());
    TrProp track(AMSPoint(position[0],position[1],position[2]),
                 AMSDir(direction[0],direction[1],direction[2]),rigidity);
    track.Propagate(fine[0]);
    AMSPoint previous = track.GetP0();
    for (size_t i = 0; i+1 < fine.size(); ++i)
    {
        track.Propagate(fine[i+1]);
        const AMSPoint next = track.GetP0();
        HybridMaterialStep cell;
        cell.zHigh = fine[i]; cell.zLow = fine[i+1];
        cell.length = (next-previous).norm();
        if (!(cell.length > 0) || !std::isfinite(cell.length)) return result;
        const double middle = 0.5*(cell.zHigh+cell.zLow);
        while (cell.segment < 3 && middle < boundary[cell.segment+1]) ++cell.segment;
        for (int station = 0; station < 4; ++station)
        {
            if (std::abs(cell.zHigh-tofZ[station])<1e-6) cell.tofAtHigh=station;
            if (std::abs(cell.zLow-tofZ[station])<1e-6) cell.tofAtLow=station;
        }
        // Five transverse samples per cell, matching the old map quadrature.
        for (int k = 0; k < 5; ++k)
        {
            const AMSPoint p = previous + (next-previous)*((k+0.5)/5.0);
            const int ix=hist[0]->GetXaxis()->FindBin(p.x());
            const int iy=hist[0]->GetYaxis()->FindBin(p.y());
            const int iz=axis->FindBin(p.z());
            if (ix<1 || ix>hist[0]->GetNbinsX() || iy<1 || iy>hist[0]->GetNbinsY() || iz<1 || iz>axis->GetNbins())
                return result;
            for (int j=0;j<9;++j)
                cell.abundance[j] += std::max(0.0,hist[j]->GetBinContent(ix,iy,iz))*cell.length/5.0;
        }
        for (int j=0;j<9;++j) result.abundance[cell.segment][j]+=cell.abundance[j];
        result.steps.push_back(cell);
        previous=next;
    }
    // tof_leng is evaluated at the actual cluster coordinate, so normalize
    // each inter-layer geometry interval using the hit-bar center planes too.
    double geometry[3]={};
    for (const auto &cell:result.steps)
        for(int gap=0;gap<3;++gap)
            if(cell.zHigh<=tofZ[gap]+1e-6 && cell.zLow>=tofZ[gap+1]-1e-6)
                geometry[gap]+=cell.length;
    for(auto &cell:result.steps)
        for(int gap=0;gap<3;++gap)
            if(cell.zHigh<=tofZ[gap]+1e-6 && cell.zLow>=tofZ[gap+1]-1e-6)
            {
                if (!(geometry[gap]>0)) return result;
                cell.length *= std::abs(tofLength[gap+1]-tofLength[gap])/geometry[gap];
            }
    result.valid=true;
    return result;
}

bool HybridEnergyLossCalibration::SupportsCharge(int charge) const
{
    for (int segment=0;segment<nSegments;++segment)
        if (!FindBin(segment,charge,0.9)) return false;
    return true;
}

bool HybridEnergyLossCalibration::Load(const std::string &path)
{
    _bins.clear();
    std::ifstream input(path.c_str());
    if (!input)
    {
        std::cerr << "Error: Could not open hybrid calibration " << path << std::endl;
        return false;
    }

    std::string line;
    int lineNumber = 0;
    while (std::getline(input, line))
    {
        ++lineNumber;
        const std::string::size_type first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#')
            continue;

        HybridEnergyLossBin bin;
        double rmsBefore = 0;
        double rmsAfter = 0;
        double correlation = 0;
        std::string extra;
        std::istringstream values(line);
        if (!(values >> bin.charge >> bin.segment >> bin.betaLow >> bin.betaHigh >>
              bin.entries >> bin.materialScale >> bin.meanDeposit >>
              bin.fluctuationScale >> rmsBefore >> rmsAfter >> correlation) ||
            (values >> extra) || bin.charge <= 0 || bin.segment < 0 ||
            bin.segment >= nSegments || !(bin.betaHigh > bin.betaLow) ||
            bin.entries <= 0 || !std::isfinite(bin.materialScale) ||
            !std::isfinite(bin.meanDeposit) ||
            !std::isfinite(bin.fluctuationScale))
        {
            std::cerr << "Error: Invalid hybrid calibration at " << path << ":"
                      << lineNumber << std::endl;
            _bins.clear();
            return false;
        }
        _bins.push_back(bin);
    }

    if (_bins.empty())
    {
        std::cerr << "Error: No hybrid calibration bins found in " << path << std::endl;
        return false;
    }
    return true;
}

const HybridEnergyLossBin *HybridEnergyLossCalibration::FindBin(
    int segment, int charge, double beta) const
{
    const HybridEnergyLossBin *nearest = nullptr;
    double nearestDistance = std::numeric_limits<double>::max();
    for (const auto &bin : _bins)
    {
        if (bin.segment != segment || bin.charge != charge)
            continue;
        if (beta >= bin.betaLow && beta < bin.betaHigh)
            return &bin;

        const double center = 0.5 * (bin.betaLow + bin.betaHigh);
        const double distance = std::abs(beta - center);
        if (distance < nearestDistance)
        {
            nearest = &bin;
            nearestDistance = distance;
        }
    }
    return nearest;
}

double HybridEnergyLossCalibration::MaterialMean(
    double beta, double mass, int charge,
    const double elementAbundance[nElements])
{
    if (!(beta > 0) || !(mass > 0) || charge == 0)
        return std::numeric_limits<double>::quiet_NaN();
    beta = std::min(beta, 1.0 - 1e-12);
    // BetheEquation already contains the projectile-charge z^2 factor.  Pass
    // the event charge explicitly here; do not reproduce the legacy AMS
    // CalculateEnergyLoss call that hard-codes z=1 and leaves charge scaling
    // to an empirical calibration factor.
    const int projectileCharge = std::abs(charge);

    static const double chargeToMass[nElements] = {
        0.99212, 0.49954, 0.49976, 0.50002, 0.47372,
        0.47847, 0.48181, 0.49848, 0.39575};
    static const double excitationEnergy[nElements] = {
        19.2e-6, 78.0e-6, 82.0e-6, 95.0e-6, 115.0e-6,
        149.0e-6, 166.0e-6, 173.0e-6, 823.0e-6};
    static const double atomicMass[nElements] = {
        1.0079, 12.0107, 14.0067, 15.9994, 18.9984,
        22.9897, 26.9815, 28.0855, 207.2};

    double loss = 0;
    for (int element = 0; element < nElements; ++element)
    {
        if (!(elementAbundance[element] >= 0) ||
            !std::isfinite(elementAbundance[element]))
            return std::numeric_limits<double>::quiet_NaN();
        const double stoppingPower = 0.001 * Betalhd::BetheEquation(
            beta, chargeToMass[element], excitationEnergy[element],
            projectileCharge, 1000 * mass);
        loss += stoppingPower * elementAbundance[element] * atomicMass[element];
    }
    return loss;
}

double HybridEnergyLossCalibration::EstimateLoss(
    int segment, int charge, double beta, double mass,
    const double elementAbundance[nElements], double measuredDeposit) const
{
    const HybridEnergyLossBin *bin = FindBin(segment, charge, beta);
    if (!bin)
        return std::numeric_limits<double>::quiet_NaN();

    const double materialMean =
        MaterialMean(beta, mass, charge, elementAbundance);
    if (!std::isfinite(materialMean) || !std::isfinite(measuredDeposit))
        return std::numeric_limits<double>::quiet_NaN();

    const double baseline = bin->materialScale * materialMean;
    const double estimate = baseline +
                            bin->fluctuationScale *
                                (measuredDeposit - bin->meanDeposit);
    // Catastrophic nuclear interactions create deposit tails that are outside
    // the ionisation-loss population used for this calibration.  Keep their
    // detector signal from turning a timing correction into a multi-GeV trial
    // loss while retaining fluctuations up to twice the material baseline.
    return std::max(0.0, std::min(estimate, 3.0 * baseline));
}
