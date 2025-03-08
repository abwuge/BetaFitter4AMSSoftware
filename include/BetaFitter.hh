#ifndef __BETAFITTER_HH__
#define __BETAFITTER_HH__

#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include "TrFit.h"
#include "ParticlePropagator.hh"
#include "ParticleData.hh"

class BetaFitter
{
public:
    // TOF time resolution (ns)
    static constexpr double TOF_TIME_RESOLUTION[4] = {0.1, 0.1, 0.1, 0.1};

    /**
     * Reconstruct beta using non-linear method
     * @param particle Input particle containing track information
     * @param propagator Particle propagator for track extrapolation
     * @param measuredTimes Array of measured TOF hit times
     * @param timeErrors Array of TOF time measurement errors
     * @return Reconstructed 1/beta value
     */
    static double reconstructBeta(const ParticleData *particle,
                                  ParticlePropagator &propagator);

private:
    /**
     * Calculate chi-square for a given 1/beta value
     * @param invBeta 1/beta value to test
     * @param propagator Particle propagator
     * @param measuredTimes Measured TOF hit times
     * @param timeErrors TOF time measurement errors
     * @return Chi-square value
     */
    static double calculateChi2(const double *invBeta,
                                ParticlePropagator &propagator,
                                const ParticleData &data);
};

#endif