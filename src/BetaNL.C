
//////////////////////////////////////////////////////////////////////////
///
///\file  BetaNL.C
///\brief Source file of BetaNLPars & BetaNL class
///
///\date  2025/3/15 HW First stable version
///$Date: 2023/12/08 15:32:50 $
///
///$Revision: 1.0 $
///
//////////////////////////////////////////////////////////////////////////

#include "BetaNL.h"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <TMath.h>

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
    initPathLength();
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
    initPathLength();
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
    initPathLength();
}

void BetaNLPars::initPathLength()
{
    _pathLength[0] -= _pathLength[1];
    _pathLength[3] -= _pathLength[2];
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
    const double initialScale = 1;
    const double lowerScale = initialScale - scaleRange;
    const double upperScale = initialScale + scaleRange;
    minimizer->SetLimitedVariable(0, "scale", initialScale, 0.1 * scaleRange, lowerScale, upperScale);

    const double timeError = _pars->_hitTimeError[0];
    const double initialTimeOffset = (_pars->_hitTime[0] + _pars->_hitTime[3]) / 2;
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

    const auto &paths = _pars->_pathLength;

    if (beta >= 1 - 1e-10)
    {
        hitTimes[1] = paths[1] * inv_c / beta;
        hitTimes[2] = paths[2] * inv_c / beta;

        hitTimes[0] = hitTimes[1] + paths[0] * inv_c / beta;
        hitTimes[3] = hitTimes[2] + paths[3] * inv_c / beta;

        return hitTimes;
    }

    const double mass = _pars->_mass;
    const double massSquared = _pars->_massSquared;
    const double energy = mass / TMath::Sqrt(1 - beta * beta);

    const auto &deps = _pars->_energyDeposited;

    hitTimes[1] = paths[1] * inv_c / beta;
    hitTimes[2] = paths[2] * inv_c / beta;

    const double energy1 = energy + deps[1] * _energyLossScale;
    const double energy2 = energy - deps[2] * _energyLossScale;

    const double inv_beta1 = 1.0 / std::sqrt(1.0 - massSquared / (energy1 * energy1));
    const double inv_beta2 = 1.0 / std::sqrt(1.0 - massSquared / (energy2 * energy2));

    hitTimes[0] = hitTimes[1] + paths[0] * inv_c * inv_beta1;
    hitTimes[3] = hitTimes[2] + paths[3] * inv_c * inv_beta2;

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
    for (size_t i = 0; i < BetaNLPars::nTOF; ++i)
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

    const auto &hitTimeReconstructed = propagate(mcBeta).data();
    const auto &hitTimeMeasured = _pars->_hitTime.data();
    const auto &hitTimeError = _pars->_hitTimeError.data();

    double chi2 = 0;
    for (size_t i = 0; i < BetaNLPars::nTOF; ++i)
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
    const double upperInvBeta = 10;   // beta > 0.1
    const double initialInvBeta = TMath::Range(lowerInvBeta, upperInvBeta, 1 / _pars->_beta);
    minimizer->SetLimitedVariable(0, "invBeta", initialInvBeta, 1e-5, lowerInvBeta, upperInvBeta);

    const double timeError = _pars->_hitTimeError[0];
    const double initialTimeOffset = (_pars->_hitTime[0] + _pars->_hitTime[3]) / 2;
    const double lowerTimeOffset = initialTimeOffset - 5 * timeError;
    const double upperTimeOffset = initialTimeOffset + 5 * timeError;
    minimizer->SetLimitedVariable(1, "timeOffset", initialTimeOffset, 0.1 * timeError, lowerTimeOffset, upperTimeOffset);

    minimizer->Minimize();

    _invBeta = std::make_shared<double>(minimizer->X()[0]);
    _timeOffset = minimizer->X()[1];

    return *_invBeta;
}