#include "BetaFitter.hh"
#include <algorithm>
#include <TMath.h>

double BetaFitter::getInitialBeta(const double measuredTimes[4],
                                 const double pathLengths[4]) 
{
    double sumXY = 0, sumX = 0, sumY = 0, sumX2 = 0;
    int nPoints = 0;

    // Linear fit of path length vs time
    for (int i = 0; i < 4; ++i) {
        if (measuredTimes[i] <= 0 || pathLengths[i] <= 0) continue;
        
        sumXY += measuredTimes[i] * pathLengths[i];
        sumX += measuredTimes[i];
        sumY += pathLengths[i];
        sumX2 += measuredTimes[i] * measuredTimes[i];
        nPoints++;
    }

    if (nPoints < 2) return 1.0; // Default to speed of light if not enough points

    // Calculate slope = c/beta
    double slope = (nPoints * sumXY - sumX * sumY) / (nPoints * sumX2 - sumX * sumX);
    return slope / ParticlePropagator::SPEED_OF_LIGHT;
}

double BetaFitter::calculateChi2(const double* invBeta,
                               ParticlePropagator& propagator,
                               const double measuredTimes[4],
                               const double timeErrors[4])
{
    // Set particle beta for propagation
    double beta = 1.0 / invBeta[0];
    if (!propagator.resetPropagator(beta))
        return 1e10 * invBeta[0];

    // Calculate TOF hits and path lengths
    double hitX[4], hitY[4], hitTime[4], pathLength[4];
    if (!propagator.PropagateToTOF(hitX, hitY, hitTime, pathLength)) {
        return 1e10 * invBeta[0];
    }

    // Calculate chi-square
    double chi2 = 0;
    int nHits = 0;
    for (int i = 0; i < 4; ++i) {
        if (measuredTimes[i] <= 0 || timeErrors[i] <= 0) continue;
        
        double dt = hitTime[i] - measuredTimes[i];
        double sigma = timeErrors[i];
        chi2 += (dt * dt) / (sigma * sigma);
        nHits++;
    }

    return nHits > 0 ? chi2 : 1e10 * invBeta[0];
}

double BetaFitter::reconstructBeta(const ParticleR* particle,
                                 ParticlePropagator& propagator,
                                 const double measuredTimes[4],
                                 const double timeErrors[4])
{
    // Get initial beta estimate from linear fit
    double pathLengths[4];
    for (int i = 0; i < 4; ++i) {
        pathLengths[i] = particle->TOFTLength[i];
    }
    double initialBetaRecip = getInitialBeta(measuredTimes, pathLengths);

    // Set up minimizer
    ROOT::Math::Minimizer* minimizer = 
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
    
    minimizer->SetMaxFunctionCalls(1000);
    minimizer->SetMaxIterations(100);
    minimizer->SetTolerance(1e-3);

    // Create chi-square function
    auto chi2Function = [&](const double* params) {
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