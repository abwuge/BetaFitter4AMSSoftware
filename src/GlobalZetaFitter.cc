#include "GlobalZetaFitter.hh"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
constexpr double kInverseSpeedOfLight = 1.0 / 29.9792458; // ns/cm
}

GlobalZetaFitter::GlobalZetaFitter(const std::vector<GlobalZetaEvent> &events,
                                   EnergyLossScaleMode energyLossScaleMode,
                                   BetaReferencePoint referencePoint,
                                   GlobalZetaTarget target)
    : _events(events),
      _energyLossScaleMode(energyLossScaleMode),
      _referencePoint(referencePoint),
      _target(target)
{
}

bool GlobalZetaFitter::PredictHitTimes(const GlobalZetaEvent &event, double zeta,
                                       std::array<double, 4> &hitTimes) const
{
    hitTimes.fill(0);
    std::array<double, 4> paths;
    std::copy(event.pathLength.begin(), event.pathLength.end(), paths.begin());
    if (_referencePoint == BetaReferencePoint::AMSCenter)
    {
        paths[0] -= paths[1];
        paths[3] -= paths[2];
    }
    else
    {
        for (int station = 3; station > 0; --station)
            paths[station] -= paths[station - 1];
        paths[0] = 0;
    }

    if (!(event.mcBeta > 0) || event.mcBeta >= 1 || !(event.mass > 0))
        return false;

    const double massSquared = event.mass * event.mass;
    double energy = event.mass / std::sqrt(1 - event.mcBeta * event.mcBeta);
    const auto scaledEnergyLoss = [&](int station)
    {
        const bool applyScale = _energyLossScaleMode == EnergyLossScaleMode::All ||
                                (_energyLossScaleMode == EnergyLossScaleMode::S1S2 && station < 2) ||
                                (_energyLossScaleMode == EnergyLossScaleMode::S2Only && station == 1);
        return event.energyDeposited[station] * (applyScale ? zeta : 1.0);
    };
    const auto inverseBetaAtEnergy = [&](double value)
    {
        if (!(value > event.mass) || !std::isfinite(value))
            return std::numeric_limits<double>::quiet_NaN();
        return 1.0 / std::sqrt(1.0 - massSquared / (value * value));
    };

    if (_referencePoint == BetaReferencePoint::AMSCenter)
    {
        const double centerInverseBeta = 1.0 / event.mcBeta;
        hitTimes[1] = paths[1] * kInverseSpeedOfLight * centerInverseBeta;
        hitTimes[2] = paths[2] * kInverseSpeedOfLight * centerInverseBeta;

        const double inverseBetaS1S2 = inverseBetaAtEnergy(energy + scaledEnergyLoss(1));
        const double inverseBetaS3S4 = inverseBetaAtEnergy(energy - scaledEnergyLoss(2));
        if (!std::isfinite(inverseBetaS1S2) || !std::isfinite(inverseBetaS3S4))
            return false;
        hitTimes[0] = hitTimes[1] + paths[0] * kInverseSpeedOfLight * inverseBetaS1S2;
        hitTimes[3] = hitTimes[2] + paths[3] * kInverseSpeedOfLight * inverseBetaS3S4;
        return std::all_of(hitTimes.begin(), hitTimes.end(),
                           [](double value) { return std::isfinite(value); });
    }

    for (int station = 1; station < 4; ++station)
    {
        energy -= scaledEnergyLoss(station - 1);
        const double inverseBeta = inverseBetaAtEnergy(energy);
        if (!std::isfinite(inverseBeta))
            return false;
        hitTimes[station] = hitTimes[station - 1] +
                            paths[station] * kInverseSpeedOfLight * inverseBeta;
    }
    return std::all_of(hitTimes.begin(), hitTimes.end(),
                       [](double value) { return std::isfinite(value); });
}

double GlobalZetaFitter::ProfiledEventChi2(const GlobalZetaEvent &event,
                                            double zeta) const
{
    std::array<double, 4> predicted;
    if (!PredictHitTimes(event, zeta, predicted))
        return std::numeric_limits<double>::infinity();

    double weightSum = 0;
    double weightedResidualSum = 0;
    for (int station = 0; station < 4; ++station)
    {
        const double sigma = 0.1544809;
        const double weight = 1.0 / (sigma * sigma);
        weightSum += weight;
        const double targetTime = _target == GlobalZetaTarget::MeasuredTime
                                      ? event.hitTime[station]
                                      : event.checkpointTruthTime[station];
        weightedResidualSum += weight * (predicted[station] - targetTime);
    }
    if (!(weightSum > 0))
        return std::numeric_limits<double>::infinity();

    const double profiledResidual = weightedResidualSum / weightSum;
    double chi2 = 0;
    for (int station = 0; station < 4; ++station)
    {
        const double targetTime = _target == GlobalZetaTarget::MeasuredTime
                                      ? event.hitTime[station]
                                      : event.checkpointTruthTime[station];
        const double residual = predicted[station] - targetTime - profiledResidual;
        const double sigma = 0.1544809;
        chi2 += residual * residual / (sigma * sigma);
    }
    return chi2;
}

double GlobalZetaFitter::Chi2(double zeta) const
{
    long double total = 0;
    for (const GlobalZetaEvent &event : _events)
    {
        const double eventChi2 = ProfiledEventChi2(event, zeta);
        if (!std::isfinite(eventChi2))
            return std::numeric_limits<double>::infinity();
        total += eventChi2;
    }
    return static_cast<double>(total);
}

bool GlobalZetaFitter::IsValidAt(const GlobalZetaEvent &event, double zeta) const
{
    std::array<double, 4> hitTimes;
    return PredictHitTimes(event, zeta, hitTimes);
}

GlobalZetaResult GlobalZetaFitter::Fit(double zetaMin, double zetaMax) const
{
    GlobalZetaResult result;
    result.entries = _events.size();
    if (_events.empty() || !(zetaMin < zetaMax))
        return result;

    const double goldenRatio = (std::sqrt(5.0) - 1.0) / 2.0;
    double lower = zetaMin;
    double upper = zetaMax;
    double left = upper - goldenRatio * (upper - lower);
    double right = lower + goldenRatio * (upper - lower);
    double leftChi2 = Chi2(left);
    double rightChi2 = Chi2(right);
    for (int iteration = 0; iteration < 60; ++iteration)
    {
        if (leftChi2 < rightChi2)
        {
            upper = right;
            right = left;
            rightChi2 = leftChi2;
            left = upper - goldenRatio * (upper - lower);
            leftChi2 = Chi2(left);
        }
        else
        {
            lower = left;
            left = right;
            leftChi2 = rightChi2;
            right = lower + goldenRatio * (upper - lower);
            rightChi2 = Chi2(right);
        }
    }

    result.zeta = 0.5 * (lower + upper);
    result.chi2 = Chi2(result.zeta);
    result.chi2PerEvent = result.chi2 / result.entries;
    result.valid = std::isfinite(result.chi2);
    return result;
}
