#ifndef __PARTICLEDATA_HH__
#define __PARTICLEDATA_HH__

#include <vector>

struct ParticleData
{
    // Particle properties
    float mass;     // Mass
    float charge;   // Charge (Always positive)
    float momentum; // Momentum (positive if charge > 0, negative if charge < 0)

    // Particle direction
    float Theta;
    float Phi;

    // Particle hit information
    std::vector<float> hitX;
    std::vector<float> hitY;
    std::vector<float> hitZ;
    std::vector<float> hitTime;
    std::vector<float> hitTimeError;

    // MC truth information
    float mcBeta;       // MC truth beta
    float mcMomentum;   // MC truth momentum
    float mcMass;       // MC truth mass
    int mcPdgId;        // MC particle PDG ID
    bool isMC;          // Flag to indicate if it's MC particle

    // Constructor to initialize vectors and particle properties with appropriate default values
    ParticleData() : mass(0.0f),
                     charge(0.0f),
                     momentum(0.0f),
                     Theta(0.0f),
                     Phi(0.0f),
                     hitX(4, 0.0f),
                     hitY(4, 0.0f),
                     hitZ(4, 0.0f),
                     hitTime(4, 0.0f),
                     hitTimeError(4, 0.0f),
                     mcBeta(0.0f),
                     mcMomentum(0.0f),
                     mcMass(0.0f),
                     mcPdgId(0),
                     isMC(false) {}
};

#endif // __PARTICLEDATA_HH__