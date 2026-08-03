#ifndef __GLOBALZETAFITTER_HH__
#define __GLOBALZETAFITTER_HH__

#include <array>
#include <cstddef>
#include <vector>

#include "BetaNL.h"

struct GlobalZetaEvent
{
    float mcBeta = 0;
    float mass = 0;
    std::array<float, 4> energyDeposited = {};
    std::array<float, 9> trackerEnergyDeposited = {};
    std::array<float, 4> hitTime = {};
    std::array<float, 4> checkpointTruthTime = {};
    std::array<float, 4> pathLength = {};
    std::array<std::array<float, 3>, 4> tofPosition = {};
    std::array<std::array<float, 3>, 9> trackerPosition = {};
};

enum class GlobalZetaTarget
{
    MeasuredTime,
    CheckpointTruth
};

struct GlobalZetaResult
{
    bool valid = false;
    double zeta = 0;
    double trackerEnergyLossScale = 1;
    double chi2 = 0;
    double chi2PerEvent = 0;
    long long entries = 0;
};

class GlobalZetaFitter
{
public:
    GlobalZetaFitter(const std::vector<GlobalZetaEvent> &events,
                     EnergyLossScaleMode energyLossScaleMode,
                     BetaReferencePoint referencePoint,
                     GlobalZetaTarget target = GlobalZetaTarget::MeasuredTime,
                     double trackerEnergyLossScale = 1.0);

    GlobalZetaResult Fit(double zetaMin, double zetaMax) const;
    GlobalZetaResult FitJoint(double zetaMin, double zetaMax,
                              double trackerScaleMin, double trackerScaleMax) const;
    double Chi2(double zeta) const;
    double Chi2(double zeta, double trackerEnergyLossScale) const;
    bool IsValidAt(const GlobalZetaEvent &event, double zeta) const;
    bool IsValidAt(const GlobalZetaEvent &event, double zeta,
                   double trackerEnergyLossScale) const;
    const std::vector<GlobalZetaEvent> &Events() const { return _events; }

private:
    bool PredictHitTimes(const GlobalZetaEvent &event, double zeta,
                         double trackerEnergyLossScale,
                         std::array<double, 4> &hitTimes) const;
    double ProfiledEventChi2(const GlobalZetaEvent &event, double zeta,
                             double trackerEnergyLossScale) const;

    const std::vector<GlobalZetaEvent> &_events;
    EnergyLossScaleMode _energyLossScaleMode;
    BetaReferencePoint _referencePoint;
    GlobalZetaTarget _target;
    double _trackerEnergyLossScale;
};

#endif // __GLOBALZETAFITTER_HH__
