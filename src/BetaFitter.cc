#include "BetaFitter.hh"
#include <algorithm>
#include <TMath.h>

double BetaFitter::calculateChi2(const double *invBeta,
                                 ParticlePropagator &propagator,
                                 const double measuredTimes[4],
                                 const double timeErrors[4])
{
    // Set particle beta for propagation
    double beta = 1.0 / invBeta[0];
    if (!propagator.resetPropagator(beta))
        return 1e10 * invBeta[0];

    // Calculate TOF hits and path lengths
    double hitX[4], hitY[4], hitTime[4], pathLength[4];
    if (!propagator.PropagateToTOF(hitX, hitY, hitTime, pathLength))
        return 1e10 * invBeta[0];

    double minTime = *std::min_element(hitTime, hitTime + 4);
    for (int i = 0; i < 4; ++i)
        hitTime[i] -= minTime;

    // Calculate chi-square
    double chi2 = 0;
    for (int i = 0; i < 4; ++i)
    {
        if (measuredTimes[i] == -1)
            continue;
        double dt = hitTime[i] - measuredTimes[i];
        double sigma = timeErrors[i];
        chi2 += (dt * dt) / (sigma * sigma);
    }

    return std::isinf(chi2) ? 1e10 * invBeta[0] : chi2;
}

double BetaFitter::reconstructBeta(const ParticleData *particle,
                                   ParticlePropagator &propagator,
                                   const double measuredTimesOri[4],
                                   const double timeErrors[4])
{
    double measuredTimes[4] = {-1, -1, -1, -1};
    double minTime = 1e10;
    for (int i = 0; i < 4; ++i)
        if (measuredTimesOri[i] != -1)
            minTime = std::min(minTime, measuredTimesOri[i]);

    for (int i = 0; i < 4; ++i)
        if (measuredTimesOri[i] != -1)
            measuredTimes[i] = measuredTimesOri[i] - minTime;

    double initialBetaRecip = 1 / particle->beta;

    // Set up minimizer
    ROOT::Math::Minimizer *minimizer =
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    // minimizer->SetMaxFunctionCalls(1000);
    // minimizer->SetMaxIterations(100);
    // minimizer->SetTolerance(1e-3);

    // Create chi-square function
    auto chi2Function = [&](const double *params)
    {
        return calculateChi2(params, propagator, measuredTimes, timeErrors);
    };
    ROOT::Math::Functor functor(chi2Function, 1);
    minimizer->SetFunction(functor);

    // Set parameter limits and initial value
    double lowerLimit = 1.0 + 1e-4;  // beta > 0
    double upperLimit = 10.0 - 1e-4; // reasonable upper limit
    initialBetaRecip = TMath::Range(initialBetaRecip, lowerLimit, upperLimit);
    minimizer->SetLimitedVariable(0, "betaReciprocal", initialBetaRecip, 1e-5, lowerLimit, upperLimit);

    // Do minimization
    minimizer->Minimize();
    double betaReciprocal = minimizer->X()[0];

    delete minimizer;
    return betaReciprocal;
}