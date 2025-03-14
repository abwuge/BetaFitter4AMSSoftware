#include "BetaFitter.hh"
#include <algorithm>
#include <TMath.h>

int BetaFitter::fitOption = 0;

double BetaFitter::calculateChi2(const double *params,
                                 ParticlePropagator &propagator,
                                 const ParticleData &data)
{
    // Set particle beta for propagation
    double beta = 1.0 / params[0];
    propagator.resetPropagator(beta);

    // Set energy loss scale if fitOption == 2
    if (fitOption == 2)
    {
        propagator.SetEnergyLossScale(params[1]);
    }

    // Calculate TOF hits and path lengths
    double trackerHitX[ParticleData::TRACKER_MAX_HITS]{};
    double trackerHitY[ParticleData::TRACKER_MAX_HITS]{};
    double TOFHitTime[ParticleData::TOF_MAX_HITS]{};
    if (!propagator.PropagateToTOF(trackerHitX, trackerHitY, TOFHitTime))
        return 1e10 * params[0];

    double minTime = *std::min_element(TOFHitTime, TOFHitTime + 4);
    for (int i = 0; i < 4; ++i)
        TOFHitTime[i] -= minTime;

    // Calculate chi-square
    double chi2 = 0;
    for (int i = 0; i < ParticleData::TOF_MAX_HITS; ++i)
    {
        if (data.TOF_hitTime[i] == -1)
            continue;
        double dt = TOFHitTime[i] - data.TOF_hitTime[i];
        double sigma = data.TOF_hitTimeError[i];
        chi2 += (dt * dt) / (sigma * sigma);
    }

    if (fitOption > 0)
    {
        // Only consider tracker hits inside the TOF region
        for (int i = 1; i < ParticleData::TRACKER_MAX_HITS - 1; ++i)
        {
            double sigma = data.TRACKER_hitError[i];

            double dx = trackerHitX[i] - data.TRACKER_hitX[i];
            double dy = trackerHitY[i] - data.TRACKER_hitY[i];

            chi2 += (dx * dx + dy * dy) / (sigma * sigma);
        }
    }

    return std::isinf(chi2) ? 1e10 * params[0] : chi2;
}

double BetaFitter::reconstructBeta(const ParticleData *particle,
                                   ParticlePropagator &propagator)
{
    double initialBetaRecip = 1. / particle->betaLinear;

    // Set up minimizer
    ROOT::Math::Minimizer *minimizer =
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    // Create chi-square function
    auto chi2Function = [&](const double *params)
    {
        return calculateChi2(params, propagator, *particle);
    };

    // Set number of parameters based on fitOption
    int nParams = (fitOption == 2) ? 2 : 1;
    ROOT::Math::Functor functor(chi2Function, nParams);
    minimizer->SetFunction(functor);

    // Set parameter limits and initial values
    double lowerLimit = 1.0 + 1e-10;  // beta < 1
    double upperLimit = 10.0; // beta > 0.1
    initialBetaRecip = TMath::Range(lowerLimit, upperLimit, initialBetaRecip);
    minimizer->SetLimitedVariable(0, "betaReciprocal", initialBetaRecip, 1e-5, lowerLimit, upperLimit);

    if (fitOption == 2)
        minimizer->SetLimitedVariable(1, "energyLossScale", 1.0, 1e-3, 1, 3);

    // Do minimization
    minimizer->Minimize();
    double betaReciprocal = minimizer->X()[0];

    if (fitOption == 2 && betaReciprocal - 1. / particle->mcBeta > 0.2)
        printf("scale: %f\n", minimizer->X()[1]);

    delete minimizer;
    return betaReciprocal;
}