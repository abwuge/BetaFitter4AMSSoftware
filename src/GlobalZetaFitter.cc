#include "GlobalZetaFitter.hh"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{
constexpr double kInverseSpeedOfLight = 1.0 / 29.9792458; // ns/cm

using Point = std::array<double, 3>;

double distance(const Point &first, const Point &second)
{
    const double dx = first[0] - second[0];
    const double dy = first[1] - second[1];
    const double dz = first[2] - second[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
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

    std::array<Point, 4> tofPoints;
    std::array<Point, 9> trackerPoints;
    for (int coordinate = 0; coordinate < 3; ++coordinate)
    {
        for (int station = 0; station < 4; ++station)
            tofPoints[station][coordinate] = event.tofPosition[station][coordinate];
        for (int layer = 0; layer < 9; ++layer)
            trackerPoints[layer][coordinate] = event.trackerPosition[layer][coordinate];
    }

    if (!(tofPoints[1][2] > trackerPoints[1][2] &&
          trackerPoints[1][2] > trackerPoints[2][2] &&
          trackerPoints[2][2] > trackerPoints[3][2] &&
          trackerPoints[3][2] > trackerPoints[4][2] &&
          trackerPoints[4][2] > 0 && 0 > trackerPoints[5][2] &&
          trackerPoints[5][2] > trackerPoints[6][2] &&
          trackerPoints[6][2] > trackerPoints[7][2] &&
          trackerPoints[7][2] > tofPoints[2][2]))
        return false;

    const double centerFraction = trackerPoints[4][2] /
                                  (trackerPoints[4][2] - trackerPoints[5][2]);
    Point center;
    for (int coordinate = 0; coordinate < 3; ++coordinate)
        center[coordinate] = trackerPoints[4][coordinate] +
                             centerFraction * (trackerPoints[5][coordinate] -
                                               trackerPoints[4][coordinate]);

    const auto segmentTime = [&](double length, double segmentEnergy)
    {
        const double inverseBeta = inverseBetaAtEnergy(segmentEnergy);
        return std::isfinite(inverseBeta)
                   ? length * kInverseSpeedOfLight * inverseBeta
                   : std::numeric_limits<double>::quiet_NaN();
    };

    const auto trackerPathScale = [&](const Point &first, const Point &last,
                                      int firstLayer, int lastLayer,
                                      double targetLength)
    {
        double geometricLength = distance(first, trackerPoints[firstLayer]);
        for (int layer = firstLayer; layer != lastLayer; layer += firstLayer < lastLayer ? 1 : -1)
            geometricLength += distance(trackerPoints[layer],
                                        trackerPoints[layer + (firstLayer < lastLayer ? 1 : -1)]);
        geometricLength += distance(trackerPoints[lastLayer], last);
        return geometricLength > 0 ? targetLength / geometricLength : 0.0;
    };

    if (_referencePoint == BetaReferencePoint::AMSCenter)
    {
        const double upperScale = trackerPathScale(center, tofPoints[1], 4, 1,
                                                   std::abs(paths[1]));
        const double lowerScale = trackerPathScale(center, tofPoints[2], 5, 7,
                                                   std::abs(paths[2]));
        if (!(upperScale > 0) || !(lowerScale > 0))
            return false;

        double upperEnergy = energy;
        double upperTime = 0;
        Point current = center;
        for (int layer = 4; layer >= 1; --layer)
        {
            upperTime += segmentTime(distance(current, trackerPoints[layer]) * upperScale,
                                     upperEnergy);
            upperEnergy += event.trackerEnergyDeposited[layer];
            current = trackerPoints[layer];
        }
        upperTime += segmentTime(distance(current, tofPoints[1]) * upperScale, upperEnergy);
        hitTimes[1] = -upperTime;

        double lowerEnergy = energy;
        double lowerTime = 0;
        current = center;
        for (int layer = 5; layer <= 7; ++layer)
        {
            lowerTime += segmentTime(distance(current, trackerPoints[layer]) * lowerScale,
                                     lowerEnergy);
            lowerEnergy -= event.trackerEnergyDeposited[layer];
            current = trackerPoints[layer];
        }
        lowerTime += segmentTime(distance(current, tofPoints[2]) * lowerScale, lowerEnergy);
        hitTimes[2] = lowerTime;

        hitTimes[0] = hitTimes[1] + segmentTime(paths[0], upperEnergy + scaledEnergyLoss(1));
        hitTimes[3] = hitTimes[2] + segmentTime(paths[3], lowerEnergy - scaledEnergyLoss(2));
        return std::all_of(hitTimes.begin(), hitTimes.end(),
                           [](double value) { return std::isfinite(value); });
    }

    energy -= scaledEnergyLoss(0);
    hitTimes[1] = segmentTime(paths[1], energy);
    energy -= scaledEnergyLoss(1);

    const double middleScale = trackerPathScale(tofPoints[1], tofPoints[2], 1, 7,
                                                paths[2]);
    if (!(middleScale > 0))
        return false;
    Point current = tofPoints[1];
    for (int layer = 1; layer <= 7; ++layer)
    {
        hitTimes[2] += segmentTime(distance(current, trackerPoints[layer]) * middleScale,
                                   energy);
        energy -= event.trackerEnergyDeposited[layer];
        current = trackerPoints[layer];
    }
    hitTimes[2] += segmentTime(distance(current, tofPoints[2]) * middleScale, energy);
    hitTimes[2] += hitTimes[1];

    energy -= scaledEnergyLoss(2);
    hitTimes[3] = hitTimes[2] + segmentTime(paths[3], energy);
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
