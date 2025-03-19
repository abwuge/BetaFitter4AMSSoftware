
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
    AMSPoint pos,
    AMSDir dir,
    double beta,
    double mass,
    int charge,
    std::vector<double> zTOF,
    std::vector<double> energyDeposited,
    std::vector<double> hitTime,
    std::vector<double> hitTimeError)
    : BetaNLPars(pos, dir, beta, mass, charge)
{
    _zTOF = zTOF;
    _energyDeposited = energyDeposited;
    _hitTime = hitTime;
    _hitTimeError = hitTimeError;
}

BetaNLPars::BetaNLPars(
    const AMSPoint pos,
    const AMSDir dir,
    const double beta,
    const double mass,
    const int charge,
    const double zTOF[nTOF],
    const double energyDeposited[nTOF],
    const double hitTime[nTOF],
    const double hitTimeError[nTOF])
    : BetaNLPars(pos, dir, beta, mass, charge)
{
    _zTOF.assign(zTOF, zTOF + nTOF);
    _energyDeposited.assign(energyDeposited, energyDeposited + nTOF);
    _hitTime.assign(hitTime, hitTime + nTOF);
    _hitTimeError.assign(hitTimeError, hitTimeError + nTOF);
}

BetaNLPars::BetaNLPars(
    const AMSPoint pos,
    const AMSDir dir,
    const double beta,
    const double mass,
    const int charge,
    const float zTOF[nTOF],
    const float energyDeposited[nTOF],
    const float hitTime[nTOF],
    const float hitTimeError[nTOF])
    : BetaNLPars(pos, dir, beta, mass, charge)
{
    _zTOF.assign(zTOF, zTOF + nTOF);
    _energyDeposited.assign(energyDeposited, energyDeposited + nTOF);
    _hitTime.assign(hitTime, hitTime + nTOF);
    _hitTimeError.assign(hitTimeError, hitTimeError + nTOF);
}

BetaNLPars::BetaNLPars(AMSPoint pos, AMSDir dir, double beta, double mass, int charge)
    : _pos(pos), _dir(dir), _beta(beta), _mass(mass), _charge(charge)
{
    if (charge == 0)
        std::cerr << "Info in <BetaNLPars::BetaNLPars>: Charge is zero!" << std::endl;
}

double BetaNL::BetaZ0()
{
    double beta = Beta();
    double energy = _pars->_mass / TMath::Sqrt(1 - beta * beta);

    energy -= (_pars->_energyDeposited[0] + _pars->_energyDeposited[1]) * _energyLossScale;

    double momentum = TMath::Sqrt(energy * energy - _pars->_mass * _pars->_mass);
    return momentum / energy;
}

double BetaNL::BetaS4()
{
    double beta = Beta();
    double energy = _pars->_mass / TMath::Sqrt(1 - beta * beta);

    energy -= (_pars->_energyDeposited[0] + _pars->_energyDeposited[1] + _pars->_energyDeposited[2]) * _energyLossScale;

    double momentum = TMath::Sqrt(energy * energy - _pars->_mass * _pars->_mass);
    return momentum / energy;
}

double BetaNL::EnergyLossScale(double mcBeta)
{
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    ROOT::Math::Functor functor(
        [&](const double *params)
        { return scaleChi2(params, mcBeta); },
        2);
    minimizer->SetFunction(functor);

    double scaleRange = 10;
    double initialScale = 2;
    double lowerScale = initialScale - scaleRange;
    double upperScale = initialScale + scaleRange;
    minimizer->SetLimitedVariable(0, "scale", initialScale, 0.1 * scaleRange, lowerScale, upperScale);

    double timeError = _pars->_hitTimeError[0];
    double initialTimeOffset = _pars->_hitTime[0];
    double lowerTimeOffset = initialTimeOffset - 5 * timeError;
    double upperTimeOffset = initialTimeOffset + 5 * timeError;
    minimizer->SetLimitedVariable(1, "timeOffset", initialTimeOffset, 0.1 * timeError, lowerTimeOffset, upperTimeOffset);

    minimizer->Minimize();

    return _energyLossScale = minimizer->X()[0];
}

/**
 * @brief Get a TrProp object with the particle parameters
 */
TrProp BetaNL::Propagator() const
{
    TrProp propagator = TrProp(_pars->Pos(), _pars->Dir());
    propagator.SetMassChrg(_pars->Mass(), _pars->Charge());

    return propagator;
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
std::vector<double> BetaNL::propagate(const double beta) const
{
    std::vector<double> hitTimes(BetaNLPars::nTOF, 0.0);
    TrProp propagator = Propagator();

    double mass = _pars->Mass();
    double mass2 = mass * mass;

    double invCharge;
    if (_pars->Charge() == 0 || TMath::Abs(beta - 1) < 1e-10)
        invCharge = 0;
    else
        invCharge = 1.0 / _pars->Charge();

    double energy;
    if (TMath::Abs(beta - 1) < 1e-10)
        energy = mass;
    else if (beta > 1)
        energy = mass / TMath::Sqrt(1 - (2 - beta) * (2 - beta));
    else
        energy = mass / TMath::Sqrt(1 - beta * beta);

    const double *const energyDeposited = _pars->_energyDeposited.data();
    const double *const zTOF = _pars->_zTOF.data();

    for (size_t i = 1; i < BetaNLPars::nTOF; ++i)
    {
        // Update after passing last TOF layer
        // ------------------------------------------
        // Update particle energy
        energy -= energyDeposited[i - 1] * _energyLossScale;
        if (energy <= mass) // Particle stopped
            break;

        // Update particle rigidity
        double momentum = TMath::Sqrt(energy * energy - mass2);
        double rigidity = momentum * invCharge;
        propagator.SetRigidity(rigidity);

        // Propagate to this TOF layer
        // ------------------------------------------
        // Calculate path length
        // TODO: In TrProp::Propagate -> TrProp::Interpolate -> TrProp::VCFitPar,
        // change `if (fabs(h) > steps && imat++ == 0) { ... }` to `if (fabs(h) > steps && imat++ == 0 && m55) { ... }`
        // may improve the performance and do NOT change the result.
        double length = propagator.Propagate(zTOF[i]);
        if (length < 0)
            break;

        // Calculate hit time
        hitTimes[i] = hitTimes[i - 1] + length / (BetaNLPars::SPEED_OF_LIGHT * momentum / energy);
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

    const double *const hitTimeReconstructed = propagate(1 / invBeta).data();
    const double *const hitTimeMeasured = _pars->_hitTime.data();
    const double *const hitTimeError = _pars->_hitTimeError.data();

    double chi2 = 0;
    for (size_t i = 0; i < BetaNLPars::nTOF; ++i)
    {
        if (hitTimeMeasured[i] == -1) // Skip missing hit times
            continue;
        double dt = hitTimeReconstructed[i] - (hitTimeMeasured[i] - _timeOffset);
        double sigma = hitTimeError[i];
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

    const double *const hitTimeReconstructed = propagate(mcBeta).data();
    const double *const hitTimeMeasured = _pars->_hitTime.data();
    const double *const hitTimeError = _pars->_hitTimeError.data();

    double chi2 = 0;
    for (size_t i = 0; i < BetaNLPars::nTOF; ++i)
    {
        if (hitTimeMeasured[i] == -1)
            continue;
        double dt = hitTimeReconstructed[i] - (hitTimeMeasured[i] - _timeOffset);
        double sigma = hitTimeError[i];
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

    /**
     * ALLOW BETA > 1
     * ----------------
     * if beta > 1, we will treat the particle as having the opposite charge
     * and consider its actual beta to be 1 - (beta - 1).
     */
    double lowerInvBeta = 0.6; // beta < 1.11
    double upperInvBeta = 3;   // beta > 0.33
    double initialInvBeta = TMath::Range(lowerInvBeta, upperInvBeta, 1 / _pars->_beta);
    minimizer->SetLimitedVariable(0, "invBeta", initialInvBeta, 1e-5, lowerInvBeta, upperInvBeta);

    double timeError = _pars->_hitTimeError[0];
    double initialTimeOffset = _pars->_hitTime[0];
    double lowerTimeOffset = initialTimeOffset - 5 * timeError;
    double upperTimeOffset = initialTimeOffset + 5 * timeError;
    minimizer->SetLimitedVariable(1, "timeOffset", initialTimeOffset, 0.1 * timeError, lowerTimeOffset, upperTimeOffset);

    minimizer->Minimize();

    _invBeta = std::make_shared<double>(minimizer->X()[0]);
    _timeOffset = minimizer->X()[1];

    return *_invBeta;
}