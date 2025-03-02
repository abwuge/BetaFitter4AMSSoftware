#ifndef __ParticlePropagator__
#define __ParticlePropagator__

#include "TrFit.h"
#include "Betalhd.h"
#include "AMSPlane.h"

class ParticlePropagator : public TrProp
{
public:
    // Physical constants
    static constexpr double SPEED_OF_LIGHT = 29.9792458; // Speed of light in cm/ns
    // Z coordinates of TOF layers (cm)
    static constexpr double TOF_Z[4] = {65.2, 62.1, -62.1, -65.2};

private:
    AMSPoint _initPos;      // Initial position of the particle
    AMSDir _initDir;        // Initial direction of the particle
    int _chrgSign;          // Sign of the particle charge
    
    AMSPoint _hitPoints[4]; // Hit positions on TOF layers
    AMSDir _hitDirs[4];     // Track directions at hit points

    double _momentum; // Momentum of the particle
    double _energy;   // Energy of the particle
    double _mass;     // Mass of the particle
    double _charge;   // Charge of the particle
    double _beta;     // Beta (v/c) of the particle

public:
    /**
     * Initialize particle propagator with given parameters
     * @param pos Initial position of the particle
     * @param dir Initial direction of the particle
     * @param momentum Initial momentum of the particle (GeV/c)
     * @param mass Mass of the particle (GeV/c^2)
     * @param charge Charge of the particle (must not be 0)
     * @throws runtime_error if charge is 0
     */
    ParticlePropagator(const AMSPoint &pos,
                       const AMSDir &dir,
                       double momentum,
                       double mass,
                       int charge);

    /**
     * Reset propagator state with new parameters
     * @param beta New beta of the particle
     * @return true if reset successful, false if beta is invalid
     */
    bool resetPropagator(double beta);

    /**
     * Propagate particle to specified z position with energy loss
     * @param z_target Target z coordinate to propagate to
     * @return Beta value at target position, -1 if propagation failed
     */
    double PropagateToZ(double z_target);

    /**
     * Propagate particle through all TOF layers
     * @param hitX Array to store x coordinates at each TOF layer [output]
     * @param hitY Array to store y coordinates at each TOF layer [output]
     * @param hitTime Array to store hit times at each TOF layer [output]
     * @param pathLength Array to store cumulative path lengths to each TOF layer [output]
     * @return true if propagation successful through all layers
     */
    bool PropagateToTOF(double hitX[4], double hitY[4],
                        double hitTime[4], double pathLength[4]);

    /**
     * Get particle's current beta (v/c)
     * @return Beta value
     */
    double GetBeta() const;

    /**
     * Get particle's current gamma factor
     * @return Gamma value
     */
    double GetGamma() const;

    /**
     * Get particle's current energy
     * @return Energy in GeV
     */
    double GetEnergy() { return _energy; }

    /**
     * Get particle's current momentum
     * @return Momentum in GeV/c
     */
    double GetMomentum() const;

    /**
     * Get array of hit positions on TOF layers
     * @return Array of hit positions
     */
    const AMSPoint *GetHitPoints() const { return _hitPoints; }

    /**
     * Get array of particle directions at TOF hits
     * @return Array of directions
     */
    const AMSDir *GetHitDirs() const { return _hitDirs; }

private:
    /**
     * Update particle kinematics considering energy loss
     * @param start_point Starting position for energy loss calculation
     * @param direction Particle direction for energy loss calculation
     * @param z_target Target z position for energy loss calculation
     */
    void UpdateWithEnergyLoss(const AMSPoint &start_point,
                              const AMSDir &direction,
                              double z_target);
};

#endif