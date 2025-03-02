#ifndef __PARTICLEDATA_HH__
#define __PARTICLEDATA_HH__

struct ParticleData
{
    // Particle properties
    float mass;     // Mass
    float charge;   // Charge (Always positive)
    float momentum; // Momentum (positive if charge > 0, negative if charge < 0)
    float beta;     // Reconstructed beta value

    // Particle direction
    float Theta;
    float Phi;

    // Particle hit information
    static const int MAX_HITS = 4; // Maximum number of hits
    float hitX[MAX_HITS];
    float hitY[MAX_HITS];
    float hitZ[MAX_HITS];
    float hitTime[MAX_HITS];
    float hitTimeError[MAX_HITS];

    // MC truth information
    float mcBeta;     // MC truth beta
    float mcMomentum; // MC truth momentum
    float mcMass;     // MC truth mass
    int mcPdgId;      // MC particle PDG ID
    bool isMC;        // Flag to indicate if it's MC particle

    // Constructor to initialize arrays and particle properties with appropriate default values
    ParticleData() : mass(0.0f),
                     charge(0.0f),
                     momentum(0.0f),
                     beta(0.0f),
                     Theta(0.0f),
                     Phi(0.0f),
                     mcBeta(0.0f),
                     mcMomentum(0.0f),
                     mcMass(0.0f),
                     mcPdgId(0),
                     isMC(false)
    {
        // Initialize arrays with zeros
        for (int i = 0; i < MAX_HITS; i++)
        {
            hitX[i] = 0.0f;
            hitY[i] = 0.0f;
            hitZ[i] = 0.0f;
            hitTime[i] = 0.0f;
            hitTimeError[i] = 0.0f;
        }
    }
};

#endif // __PARTICLEDATA_HH__