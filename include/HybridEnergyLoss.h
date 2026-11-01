#ifndef __HYBRIDENERGYLOSS_H__
#define __HYBRIDENERGYLOSS_H__

#include <string>
#include <vector>
#include <array>

struct HybridMaterialStep
{
    double zHigh = 0;
    double zLow = 0;
    double length = 0;
    int segment = 0;
    int tofAtHigh = -1;
    int tofAtLow = -1;
    std::array<double, 9> abundance = {};
};

struct HybridMaterialProfile
{
    std::vector<HybridMaterialStep> steps;
    std::array<std::array<double, 9>, 4> abundance = {};
    std::array<double, 5> boundaryZ = {};
    std::array<double, 4> tofCenterZ = {};
    double centerZ = 1.698;
    bool valid = false;
};

struct HybridEnergyLossBin
{
    int charge = 0;
    int segment = 0;
    double betaLow = 0;
    double betaHigh = 0;
    long long entries = 0;
    double materialScale = 1;
    double meanDeposit = 0;
    double fluctuationScale = 0;
};

class HybridEnergyLossCalibration
{
public:
    static constexpr int nSegments = 4;
    static constexpr int nElements = 9;

    bool Load(const std::string &path);
    bool Empty() const { return _bins.empty(); }
    bool SupportsCharge(int charge) const;

    // Integrate the native map without rounding endpoints to voxel centers.
    static bool RegionGeometry(const int tofBar[4], double boundaryZ[5],
                               double tofCenterZ[4]);

    static HybridMaterialProfile BuildProfile(
        const float position[3], const float direction[3], double rigidity,
        const int tofBar[4], const float tofLength[4],
        const float trackerZ[9], double maxStep = 1.0);

    double EstimateLoss(int segment, int charge, double beta, double mass,
                        const double elementAbundance[nElements],
                        double measuredDeposit) const;

    static double MaterialMean(double beta, double mass, int charge,
                               const double elementAbundance[nElements]);

private:
    const HybridEnergyLossBin *FindBin(int segment, int charge,
                                      double beta) const;

    std::vector<HybridEnergyLossBin> _bins;
};

#endif
