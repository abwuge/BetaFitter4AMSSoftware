#ifndef __PARTICLEDATA_HH__
#define __PARTICLEDATA_HH__

struct ParticleData
{
    // Particle properties
    float mass;       // Mass
    int charge;       // Charge
    float momentum;   // Momentum
    float betaLinear; // Reconstructed beta value

    // Particle hit information (TOF)
    static const int TOF_MAX_HITS = 4;
    float TOF_hitZ[TOF_MAX_HITS];
    float TOF_hitTime[TOF_MAX_HITS];
    float TOF_hitTimeError[TOF_MAX_HITS];

    // Particle hit information (Tracker)
    float TRACKER_dir[3]; // only load the first hit direction
    static const int TRACKER_MAX_HITS = 9;
    float TRACKER_hitX[TRACKER_MAX_HITS];
    float TRACKER_hitY[TRACKER_MAX_HITS];
    float TRACKER_hitZ[TRACKER_MAX_HITS];
    float TRACKER_hitError[TRACKER_MAX_HITS];

    // MC truth information
    bool isMC;        // Flag to indicate if it's MC particle
    int mcPdgId;      // MC particle PDG ID
    int mcCharge;     // MC truth charge
    float mcBeta;     // MC truth beta
    float mcMomentum; // MC truth momentum
    float mcMass;     // MC truth mass
    float mcCoo[3];   // MC truth initial position (cm)
    float mcDir[3];   // MC truth initial direction cosines

    // Constructor to initialize arrays and particle properties with appropriate default values
    ParticleData() : mass(0.0f),
                     charge(0),
                     momentum(0.0f),
                     betaLinear(0.0f),
                     TOF_hitTime{},
                     TOF_hitTimeError{},
                     TRACKER_dir{},
                     TRACKER_hitX{},
                     TRACKER_hitY{},
                     TRACKER_hitZ{},
                     TRACKER_hitError{},
                     isMC(false),
                     mcPdgId(0),
                     mcCharge(0),
                     mcBeta(0.0f),
                     mcMomentum(0.0f),
                     mcMass(0.0f),
                     mcCoo{},
                     mcDir{}
    {
    }
};

#endif // __PARTICLEDATA_HH__