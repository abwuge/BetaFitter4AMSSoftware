#include "BetaNL.h"

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

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
    {
        std::cerr << "BetaNLPars: Charge cannot be zero. Setting charge to 1." << std::endl;
        _charge = 1;
    }
}

double BetaNL::EnergyLossScale(double mcBeta)
{
    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    ROOT::Math::Functor functor(
        [&](const double *scale)
        { return scaleChi2(scale, mcBeta); },
        1);
    minimizer->SetFunction(functor);

    double lowerLimit = 1;
    double upperLimit = 3;
    minimizer->SetLimitedVariable(0, "scale", 2, 1e-5, lowerLimit, upperLimit);

    minimizer->Minimize();

    return _energyLossScale = minimizer->X()[0];
}

TrProp BetaNL::Propagator(const double beta) const
{
    TrProp propagator = TrProp(_pars->Pos(), _pars->Dir());
    propagator.SetMassChrg(_pars->Mass(), _pars->Charge());

    return propagator;
}

std::vector<double> BetaNL::propagate(const double beta) const
{
    std::vector<double> hitTimes(BetaNLPars::nTOF, 0.0);
    TrProp propagator = Propagator(beta);

    double mass = _pars->Mass();
    double mass2 = mass * mass;
    double energy = mass / TMath::Sqrt(1 - beta * beta);
    double invCharge = 1.0 / _pars->_charge;
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
        double length = propagator.Propagate(zTOF[i]);
        if (length < 0)
            break;

        // Calculate hit time
        hitTimes[i] = hitTimes[i - 1] + length / (BetaNLPars::SPEED_OF_LIGHT * momentum / energy);
    }

    return hitTimes;
}

double BetaNL::betaChi2(const double *params) const
{
    const double invBeta = params[0];
    const double timeOffset = params[1];

    const double *const hitTimeReconstructed = propagate(1 / invBeta).data();
    const double *const hitTimeMeasured = _pars->_hitTime.data();
    const double *const hitTimeError = _pars->_hitTimeError.data();

    double chi2 = 0;
    for (size_t i = 0; i < BetaNLPars::nTOF; ++i)
    {
        if (hitTimeMeasured[i] == -1)
            continue;
        double dt = hitTimeReconstructed[i] - (hitTimeMeasured[i] - timeOffset);
        double sigma = hitTimeError[i];
        chi2 += (dt * dt) / (sigma * sigma);
    }

    return chi2;
}

double BetaNL::scaleChi2(const double *scale, const double mcBeta)
{
    _energyLossScale = scale[0];
    const double *const hitTimeReconstructed = propagate(mcBeta).data();
    const double *const hitTimeMeasured = _pars->_hitTime.data();
    const double *const hitTimeError = _pars->_hitTimeError.data();

    double chi2 = 0;
    for (size_t i = 0; i < BetaNLPars::nTOF; ++i)
    {
        if (hitTimeMeasured[i] == -1)
            continue;
        double dt = hitTimeReconstructed[i] - hitTimeMeasured[i];
        double sigma = hitTimeError[i];
        chi2 += (dt * dt) / (sigma * sigma);
    }

    return chi2;
}

double BetaNL::reconstruct()
{
    if (_invBeta)
        return *_invBeta;

    ROOT::Math::Minimizer *minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    ROOT::Math::Functor functor(this, &BetaNL::betaChi2, 2);
    minimizer->SetFunction(functor);

    double lowerInvBeta = 1 + 1e-10; // beta < 1
    double upperInvBeta = 3;         // beta > 0.33
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