#include "BetaNL.h"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TMath.h>

#include <limits>

namespace
{
using Point = std::array<double, 3>;

double distance(const Point &first, const Point &second)
{
    const double dx = first[0] - second[0];
    const double dy = first[1] - second[1];
    const double dz = first[2] - second[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}
}

BetaNLPars::BetaNLPars(
    const double beta,
    const double mass,
    const std::vector<double> energyDeposited,
    const std::vector<double> hitTime,
    const std::vector<double> hitTimeError,
    const std::vector<double> pathLength)
    : BetaNLPars(beta, mass)
{
    _energyDeposited = energyDeposited;
    _hitTime = hitTime;
    _hitTimeError = hitTimeError;
    _pathLength = pathLength;
    init();
}

BetaNLPars::BetaNLPars(
    const double beta,
    const double mass,
    const float energyDeposited[nTOF],
    const float hitTime[nTOF],
    const float hitTimeError[nTOF],
    const float pathLength[nTOF],
    const float tofPosition[nTOF][3],
    const float trackerEnergyDeposited[9],
    const float trackerPosition[9][3])
    : BetaNLPars(beta, mass, energyDeposited, hitTime, hitTimeError, pathLength)
{
    for (int station = 0; station < nTOF; ++station)
        for (int coordinate = 0; coordinate < 3; ++coordinate)
            _tofPosition[station][coordinate] = tofPosition[station][coordinate];
    for (int layer = 0; layer < 9; ++layer)
    {
        _trackerEnergyDeposited[layer] = trackerEnergyDeposited[layer] * 1e-3;
        for (int coordinate = 0; coordinate < 3; ++coordinate)
            _trackerPosition[layer][coordinate] = trackerPosition[layer][coordinate];
    }
    _useTrackerEnergyLoss = true;
}

BetaNLPars::BetaNLPars(
    const double beta,
    const double mass,
    const double energyDeposited[nTOF],
    const double hitTime[nTOF],
    const double hitTimeError[nTOF],
    const double pathLength[nTOF])
    : BetaNLPars(beta, mass)
{
    _energyDeposited.assign(energyDeposited, energyDeposited + nTOF);
    _hitTime.assign(hitTime, hitTime + nTOF);
    _hitTimeError.assign(hitTimeError, hitTimeError + nTOF);
    _pathLength.assign(pathLength, pathLength + nTOF);
    init();
}

BetaNLPars::BetaNLPars(
    const double beta,
    const double mass,
    const float energyDeposited[nTOF],
    const float hitTime[nTOF],
    const float hitTimeError[nTOF],
    const float pathLength[nTOF])
    : BetaNLPars(beta, mass)
{
    _energyDeposited.assign(energyDeposited, energyDeposited + nTOF);
    _hitTime.assign(hitTime, hitTime + nTOF);
    _hitTimeError.assign(hitTimeError, hitTimeError + nTOF);
    _pathLength.assign(pathLength, pathLength + nTOF);
    init();
}

void BetaNLPars::init()
{
    for (auto &edep : _energyDeposited)
        edep *= 1e-3;
}

double BetaNL::EnergyLossScale(double mcBeta)
{
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    ROOT::Math::Functor functor(
        [&](const double *params)
        { return scaleChi2(params, mcBeta); },
        2);
    minimizer->SetFunction(functor);

    const double scaleRange = 10;
    const double initialScale = 2;
    const double lowerScale = initialScale - scaleRange;
    const double upperScale = initialScale + scaleRange;
    minimizer->SetLimitedVariable(0, "scale", initialScale, 0.1 * scaleRange, lowerScale, upperScale);

    const double timeError = _pars->_hitTimeError[0];
    const double initialTimeOffset = _referencePoint == BetaReferencePoint::AMSCenter
                                         ? (_pars->_hitTime[0] + _pars->_hitTime[3]) / 2
                                         : _pars->_hitTime[0];
    const double lowerTimeOffset = initialTimeOffset - 5 * timeError;
    const double upperTimeOffset = initialTimeOffset + 5 * timeError;
    minimizer->SetLimitedVariable(1, "timeOffset", initialTimeOffset, 0.1 * timeError, lowerTimeOffset, upperTimeOffset);

    minimizer->Minimize();

    return _energyLossScale = minimizer->X()[0];
}

/**
 * @note When a particle comes to rest, the hit time remains 0.
 *       This design choice serves two purposes:
 *
 *       1. A zero hit time naturally penalizes the χ² (chi-square) statistic
 *          through standard error propagation in the optimization process.
 *
 *       2. Alternative approaches using large numerical values (e.g., MAX_FLOAT)
 *          would create discontinuous jumps in the χ² landscape, causing
 *          numerical instability in Hessian matrix calculations during
 *          minimization.
 *
 * The current implementation maintains better numerical stability while
 * preserving the physical interpretation of stationary particles in the
 * tracking algorithm.
 */
std::vector<double> BetaNL::propagate(double beta) const
{
    constexpr double inv_c = 1.0 / BetaNLPars::SPEED_OF_LIGHT;
    std::vector<double> hitTimes(BetaNLPars::nTOF, 0.0);

    auto paths = _pars->_pathLength;

    if (_referencePoint == BetaReferencePoint::AMSCenter)
    {
        paths[0] -= paths[1];
        paths[3] -= paths[2];
    }
    else
    {
        for (int i = BetaNLPars::nTOF - 1; i > 0; --i)
            paths[i] -= paths[i - 1];
        paths[0] = 0;
    }

    if (beta >= 1 - 1e-10)
    {
        if (_referencePoint == BetaReferencePoint::AMSCenter)
        {
            hitTimes[1] = paths[1] * inv_c / beta;
            hitTimes[2] = paths[2] * inv_c / beta;
            hitTimes[0] = hitTimes[1] + paths[0] * inv_c / beta;
            hitTimes[3] = hitTimes[2] + paths[3] * inv_c / beta;
            return hitTimes;
        }

        for (int i = 1; i < BetaNLPars::nTOF; ++i)
            hitTimes[i] = hitTimes[i - 1] + paths[i] * inv_c / beta;

        return hitTimes;
    }

    const double mass = _pars->_mass;
    const double massSquared = _pars->_massSquared;
    double energy = mass / TMath::Sqrt(1 - beta * beta);

    const auto &deps = _pars->_energyDeposited;

    const auto scaledEnergyLoss = [&](int station)
    {
        const bool applyScale = _energyLossScaleMode == EnergyLossScaleMode::All ||
                                (_energyLossScaleMode == EnergyLossScaleMode::S1S2 && station < 2) ||
                                (_energyLossScaleMode == EnergyLossScaleMode::S2Only && station == 1);
        return deps[station] * (applyScale ? _energyLossScale : 1.0);
    };
    const auto scaledTrackerEnergyLoss = [&](int layer)
    {
        return _pars->_trackerEnergyDeposited[layer] * _trackerEnergyLossScale;
    };

    if (_pars->_useTrackerEnergyLoss)
    {
        const auto inverseBetaAtEnergy = [&](double value)
        {
            if (!(value > mass) || !std::isfinite(value))
                return std::numeric_limits<double>::quiet_NaN();
            return 1.0 / std::sqrt(1.0 - massSquared / (value * value));
        };
        const auto segmentTime = [&](double length, double segmentEnergy)
        {
            return length * inv_c * inverseBetaAtEnergy(segmentEnergy);
        };
        const auto &tofPoints = _pars->_tofPosition;
        const auto &trackerPoints = _pars->_trackerPosition;
        const auto trackerPathScale = [&](const Point &first, const Point &last,
                                          int firstLayer, int lastLayer,
                                          double targetLength)
        {
            double geometricLength = distance(first, trackerPoints[firstLayer]);
            const int step = firstLayer < lastLayer ? 1 : -1;
            for (int layer = firstLayer; layer != lastLayer; layer += step)
                geometricLength += distance(trackerPoints[layer], trackerPoints[layer + step]);
            geometricLength += distance(trackerPoints[lastLayer], last);
            return geometricLength > 0 ? std::abs(targetLength) / geometricLength : 0.0;
        };

        if (_referencePoint == BetaReferencePoint::AMSCenter)
        {
            // Integrate from the AMS center to the upper and lower TOF sides.
            const double centerFraction = trackerPoints[4][2] /
                                          (trackerPoints[4][2] - trackerPoints[5][2]);
            Point center;
            for (int coordinate = 0; coordinate < 3; ++coordinate)
                center[coordinate] = trackerPoints[4][coordinate] +
                                     centerFraction * (trackerPoints[5][coordinate] -
                                                       trackerPoints[4][coordinate]);

            const double upperScale = trackerPathScale(center, tofPoints[1], 4, 1, paths[1]);
            const double lowerScale = trackerPathScale(center, tofPoints[2], 5, 7, paths[2]);
            if (!(upperScale > 0) || !(lowerScale > 0))
                return hitTimes;

            double upperEnergy = energy;
            double upperTime = 0;
            Point current = center;
            for (int layer = 4; layer >= 1; --layer)
            {
                upperTime += segmentTime(distance(current, trackerPoints[layer]) * upperScale,
                                         upperEnergy);
                upperEnergy += scaledTrackerEnergyLoss(layer);
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
                lowerEnergy -= scaledTrackerEnergyLoss(layer);
                current = trackerPoints[layer];
            }
            lowerTime += segmentTime(distance(current, tofPoints[2]) * lowerScale, lowerEnergy);
            hitTimes[2] = lowerTime;

            hitTimes[0] = hitTimes[1] + segmentTime(paths[0], upperEnergy + scaledEnergyLoss(1));
            hitTimes[3] = hitTimes[2] + segmentTime(paths[3], lowerEnergy - scaledEnergyLoss(2));
            return hitTimes;
        }

        // BeforeTOF propagation crosses the tracker once between TOF stations 2 and 3.
        energy -= scaledEnergyLoss(0);
        hitTimes[1] = segmentTime(paths[1], energy);
        energy -= scaledEnergyLoss(1);

        const double middleScale = trackerPathScale(tofPoints[1], tofPoints[2], 1, 7, paths[2]);
        if (!(middleScale > 0))
            return hitTimes;
        Point current = tofPoints[1];
        for (int layer = 1; layer <= 7; ++layer)
        {
            hitTimes[2] += segmentTime(distance(current, trackerPoints[layer]) * middleScale,
                                       energy);
            energy -= scaledTrackerEnergyLoss(layer);
            current = trackerPoints[layer];
        }
        hitTimes[2] += segmentTime(distance(current, tofPoints[2]) * middleScale, energy);
        hitTimes[2] += hitTimes[1];

        energy -= scaledEnergyLoss(2);
        hitTimes[3] = hitTimes[2] + segmentTime(paths[3], energy);
        return hitTimes;
    }

    if (_referencePoint == BetaReferencePoint::AMSCenter)
    {
        hitTimes[1] = paths[1] * inv_c / beta;
        hitTimes[2] = paths[2] * inv_c / beta;

        const double energyAtS1S2 = energy + scaledEnergyLoss(1);
        const double energyAtS3S4 = energy - scaledEnergyLoss(2);
        const double invBetaS1S2 = 1.0 / std::sqrt(1.0 - massSquared / (energyAtS1S2 * energyAtS1S2));
        const double invBetaS3S4 = 1.0 / std::sqrt(1.0 - massSquared / (energyAtS3S4 * energyAtS3S4));

        hitTimes[0] = hitTimes[1] + paths[0] * inv_c * invBetaS1S2;
        hitTimes[3] = hitTimes[2] + paths[3] * inv_c * invBetaS3S4;
        return hitTimes;
    }

    for (int i = 1; i < BetaNLPars::nTOF; ++i)
    {
        const int station = i - 1;
        energy -= scaledEnergyLoss(station);
        const double inv_beta = 1.0 / std::sqrt(1.0 - massSquared / (energy * energy));
        hitTimes[i] = hitTimes[i - 1] + paths[i] * inv_c * inv_beta;
    }

    return hitTimes;
}

/**
 * @param params[0] Inverse beta (1/beta)
 * @param params[1] Time offset
 */
double BetaNL::betaChi2(const double *params)
{
    const double invBeta = params[0];
    _timeOffset = params[1];

    const auto &hitTimeReconstructed = propagate(1 / invBeta);
    const auto &hitTimeMeasured = _pars->_hitTime.data();
    const auto &hitTimeError = _pars->_hitTimeError.data();

    double chi2 = 0;
    for (int i = 0; i < BetaNLPars::nTOF; ++i)
    {
        if (hitTimeMeasured[i] == -1)
            continue;
        const double dt = hitTimeReconstructed[i] - (hitTimeMeasured[i] - _timeOffset);
        const double sigma = hitTimeError[i];
        chi2 += (dt * dt) / (sigma * sigma);
    }

    return chi2;
}

/**
 * @param params[0] Energy loss scale factor
 * @param params[1] Time offset
 * @param mcBeta Monte Carlo beta value
 */
double BetaNL::scaleChi2(const double *params, const double mcBeta)
{
    _energyLossScale = params[0];
    _timeOffset = params[1];

    const auto &hitTimeReconstructed = propagate(mcBeta);
    const auto &hitTimeMeasured = _pars->_hitTime.data();
    const auto &hitTimeError = _pars->_hitTimeError.data();

    double chi2 = 0;
    for (int i = 0; i < BetaNLPars::nTOF; ++i)
    {
        if (hitTimeMeasured[i] == -1)
            continue;
        const double dt = hitTimeReconstructed[i] - (hitTimeMeasured[i] - _timeOffset);
        const double sigma = hitTimeError[i];
        chi2 += (dt * dt) / (sigma * sigma);
    }

    return chi2;
}

/**
 * @brief Performs β⁻¹ reconstruction using Minuit2 optimization framework.
 *
 * Mathematical formulation:
 *
 *   χ² = ∑[(t_reco - (t_tofMeasured - timeOffset))² / hitTimeError²]
 *
 * Where:
 *   - t_reco:      Reconstructed time from particle hypothesis
 *   - timeOffset:  Detector timing calibration constant (ns)
 *   - hitTimeError: Timing resolution (σ) of the detection system
 *
 * @note Critical design choices:
 *
 * 1. Variable selection:
 *    - Uses β⁻¹ (1/β) as minimization parameter instead of β because:
 *      a) Better Hessian matrix condition number in relativistic regime
 *      b) Maintains linearity in dE/dx relationships
 *
 * @see Minuit2 documentation: https://root.cern.ch/doc/master/Minuit2Page.html
 */
double BetaNL::reconstruct()
{
    if (_invBeta)
        return *_invBeta;

    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    ROOT::Math::Functor functor(this, &BetaNL::betaChi2, 2);
    minimizer->SetFunction(functor);

    const double lowerInvBeta = 0.6; // beta < 1.67
    const double upperInvBeta = 10;  // beta > 0.1
    const double initialInvBeta = TMath::Range(lowerInvBeta, upperInvBeta, 1 / _pars->_beta);
    minimizer->SetLimitedVariable(0, "invBeta", initialInvBeta, 1e-5, lowerInvBeta, upperInvBeta);

    const double timeError = _pars->_hitTimeError[0];
    const double initialTimeOffset = _referencePoint == BetaReferencePoint::AMSCenter
                                         ? (_pars->_hitTime[0] + _pars->_hitTime[3]) / 2
                                         : _pars->_hitTime[0];
    const double lowerTimeOffset = initialTimeOffset - 5 * timeError;
    const double upperTimeOffset = initialTimeOffset + 5 * timeError;
    minimizer->SetLimitedVariable(1, "timeOffset", initialTimeOffset, 0.1 * timeError, lowerTimeOffset, upperTimeOffset);

    minimizer->Minimize();

    _invBeta = std::make_shared<double>(minimizer->X()[0]);
    _timeOffset = minimizer->X()[1];

    return *_invBeta;
}
